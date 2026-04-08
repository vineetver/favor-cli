//! Lazy cohort and chromosome handles.
//!
//! `CohortHandle::manifest`/`samples`/`chromosomes` materialize on first
//! call via `OnceLock`. `ChromosomeView` mmap's `sparse_g.bin` and parses
//! `variants.parquet` + `membership.parquet` on first method call. Both
//! are `Send + Sync` so the gene-parallel scoring loop can share a view
//! across rayon worker threads.

use std::path::PathBuf;
use std::sync::OnceLock;

use crate::error::CohortError;
use crate::output::Output;
use crate::resource::Resources;
use crate::store::query::materialize::{DenseDosageBatch, DenseDosageScratch};
use crate::store::query::{intersect_sorted, union_sorted, Region, VariantSelector};
use crate::store::Store;
use crate::types::{Chromosome, Consequence};

use super::sparse_g::SparseG;
use super::types::{as_u32_slice, CarrierBatch, SortedVcfs, VariantMetadata, VariantRow, VariantVcf};
use super::variants::{CarrierList, VariantIndex};
use super::{
    build_or_load_store_with_probe, probe as cohort_probe, read_sample_names_at, CohortId,
    CohortManifest, GenoStoreResult, StoreProbe,
};

pub struct CohortHandle<'a> {
    store: &'a Store,
    id: CohortId,
    dir: PathBuf,
    manifest: OnceLock<CohortManifest>,
    samples: OnceLock<Vec<String>>,
    chroms: OnceLock<Vec<Chromosome>>,
}

impl<'a> CohortHandle<'a> {
    pub(crate) fn new(store: &'a Store, id: CohortId, dir: PathBuf) -> Self {
        Self {
            store,
            id,
            dir,
            manifest: OnceLock::new(),
            samples: OnceLock::new(),
            chroms: OnceLock::new(),
        }
    }

    pub fn id(&self) -> &CohortId {
        &self.id
    }

    pub fn dir(&self) -> &std::path::Path {
        &self.dir
    }

    pub fn store(&self) -> &Store {
        self.store
    }

    /// Materialise `manifest.json` on first call. O(1) afterwards.
    pub fn manifest(&self) -> Result<&CohortManifest, CohortError> {
        if let Some(m) = self.manifest.get() {
            return Ok(m);
        }
        let path = self.dir.join("manifest.json");
        let s = std::fs::read_to_string(&path).map_err(|e| {
            CohortError::Resource(format!("read {}: {e}", path.display()))
        })?;
        let m: CohortManifest = serde_json::from_str(&s).map_err(|e| {
            CohortError::Resource(format!("parse {}: {e}", path.display()))
        })?;
        Ok(self.manifest.get_or_init(|| m))
    }

    /// Sample list in column order. Index `i` is the sample_id stored in
    /// sparse_g.bin carrier entries.
    pub fn samples(&self) -> Result<&[String], CohortError> {
        if let Some(s) = self.samples.get() {
            return Ok(s);
        }
        let names = read_sample_names_at(&self.dir)?;
        Ok(self.samples.get_or_init(|| names))
    }

    /// Chromosomes in manifest order, parsed once into the typed enum.
    pub fn chromosomes(&self) -> Result<&[Chromosome], CohortError> {
        if let Some(c) = self.chroms.get() {
            return Ok(c);
        }
        let m = self.manifest()?;
        let chroms: Vec<Chromosome> = m
            .chromosomes
            .iter()
            .map(|ci| ci.name.parse::<Chromosome>())
            .collect::<Result<Vec<_>, _>>()
            .map_err(CohortError::Input)?;
        Ok(self.chroms.get_or_init(|| chroms))
    }

    /// Per-chromosome view. Lazy: opens nothing until a method on the
    /// returned view is called.
    pub fn chromosome(&self, chrom: &Chromosome) -> Result<ChromosomeView<'_>, CohortError> {
        let chrom_dir = self.dir.join(format!("chromosome={}", chrom.label()));
        if !chrom_dir.is_dir() {
            return Err(CohortError::DataMissing(format!(
                "cohort {} has no chromosome {}",
                self.id.as_str(),
                chrom.label()
            )));
        }
        Ok(ChromosomeView {
            cohort: self,
            chrom: *chrom,
            dir: chrom_dir,
            sparse_g: OnceLock::new(),
            index: OnceLock::new(),
        })
    }

    /// Probe + (re)build path for the cohort. Wraps the existing free
    /// functions so the pipeline never reaches into raw `cohort::probe`.
    pub fn probe(&self, genotypes: &std::path::Path, annotations: &std::path::Path) -> StoreProbe {
        cohort_probe(&self.dir, genotypes, annotations)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn build_or_load(
        &self,
        genotypes: &std::path::Path,
        annotations: &std::path::Path,
        geno_staging_dir: &std::path::Path,
        rebuild_store: bool,
        probe_result: StoreProbe,
        res: &Resources,
        out: &dyn Output,
    ) -> Result<GenoStoreResult, CohortError> {
        if let Some(parent) = self.dir.parent() {
            std::fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", parent.display()))
            })?;
        }
        build_or_load_store_with_probe(
            genotypes,
            annotations,
            &self.dir,
            geno_staging_dir,
            rebuild_store,
            probe_result,
            res,
            out,
        )
    }
}

/// Borrowed view of one (cohort, chromosome). Methods open
/// `sparse_g.bin`, `variants.parquet`, and `membership.parquet` on the
/// first access that needs them. Both inner caches are `OnceLock` so
/// rayon parallel iteration over genes is safe.
pub struct ChromosomeView<'a> {
    cohort: &'a CohortHandle<'a>,
    chrom: Chromosome,
    dir: PathBuf,
    sparse_g: OnceLock<SparseG>,
    index: OnceLock<VariantIndex>,
}

impl<'a> ChromosomeView<'a> {
    pub fn cohort(&self) -> &CohortHandle<'a> {
        self.cohort
    }

    pub fn chromosome(&self) -> Chromosome {
        self.chrom
    }

    pub fn dir(&self) -> &std::path::Path {
        &self.dir
    }

    /// Internal: mmap'd sparse genotype matrix.
    pub fn sparse_g(&self) -> Result<&SparseG, CohortError> {
        if let Some(sg) = self.sparse_g.get() {
            return Ok(sg);
        }
        let sg = SparseG::open(&self.dir)?;
        Ok(self.sparse_g.get_or_init(|| sg))
    }

    /// Internal: parsed variant metadata + gene membership.
    pub fn index(&self) -> Result<&VariantIndex, CohortError> {
        if let Some(idx) = self.index.get() {
            return Ok(idx);
        }
        let idx = VariantIndex::load(&self.dir)?;
        Ok(self.index.get_or_init(|| idx))
    }

    /// O(1) — read from header.
    pub fn n_variants(&self) -> Result<u32, CohortError> {
        Ok(self.sparse_g()?.n_variants())
    }

    /// O(1) — read from header.
    pub fn n_samples(&self) -> Result<u32, CohortError> {
        Ok(self.sparse_g()?.n_samples())
    }

    /// O(1) direct array index. Panics on out-of-range vcf — that's a
    /// programming error, not a runtime condition.
    pub fn metadata(&self, vcf: VariantVcf) -> Result<&VariantMetadata, CohortError> {
        Ok(self.index()?.get(vcf.get()))
    }

    /// Streaming iterator over every variant in this chromosome. O(n)
    /// total, O(1) heap.
    pub fn stream(&self) -> Result<impl Iterator<Item = VariantRow<'_>> + '_, CohortError> {
        let entries = self.index()?.all_entries();
        Ok(entries
            .iter()
            .enumerate()
            .map(|(i, meta)| VariantRow { vcf: VariantVcf(i as u32), meta }))
    }

    /// Sparse carriers for one variant. O(MAC_v).
    pub fn carriers(&self, vcf: VariantVcf) -> Result<CarrierList, CohortError> {
        Ok(self.sparse_g()?.load_variant(vcf.get()))
    }

    /// Batch carriers for many variants. Slice MUST be sorted by
    /// VariantVcf so the mmap reads land sequentially. Debug-asserts.
    pub fn carriers_batch(&self, vcfs: &[VariantVcf]) -> Result<CarrierBatch, CohortError> {
        debug_assert!(
            vcfs.windows(2).all(|w| w[0] <= w[1]),
            "carriers_batch requires sorted input"
        );
        let entries = self.sparse_g()?.load_variants(as_u32_slice(vcfs));
        Ok(CarrierBatch { entries })
    }

    /// Variant_vcfs belonging to a gene, sorted, dedup'd at build time.
    /// Empty slice if the gene is not in this chromosome.
    pub fn gene_variants(&self, gene: &str) -> Result<&[VariantVcf], CohortError> {
        Ok(super::types::from_u32_slice(
            self.index()?.gene_variant_vcfs(gene),
        ))
    }

    /// Resolve a vid string to its `VariantVcf` in this chromosome.
    pub fn resolve_vid(&self, vid: &str) -> Result<Option<VariantVcf>, CohortError> {
        Ok(self.index()?.resolve_vid(vid).map(VariantVcf))
    }

    /// Build a sorted+dedup'd `SortedVcfs` scoped to this view from an
    /// already-sorted slice.
    pub fn sorted_from(&self, vcfs: Vec<VariantVcf>) -> SortedVcfs {
        SortedVcfs::new(self.cohort.id().clone(), self.chrom, vcfs)
    }

    /// Resolve a `VariantSelector` against this chromosome's metadata
    /// and membership tables. Returns a sorted, dedup'd `SortedVcfs`
    /// scoped to (cohort, chromosome).
    pub fn select(&self, sel: &VariantSelector) -> Result<SortedVcfs, CohortError> {
        let vcfs = resolve_selector(self, sel)?;
        Ok(self.sorted_from(vcfs))
    }

    /// Dense dosage matrix [n_samples × len(vcfs)] in row-major order
    /// (`data[s * len + i]` = dosage of variant `vcfs[i]` in sample `s`).
    /// The engine fills caller-owned scratch — no allocation in the
    /// hot path. Reference samples (everyone not in a carrier list) are
    /// implicitly zero, matching the sparse → dense convention.
    pub fn dense_dosages_into<'s>(
        &self,
        vcfs: &[VariantVcf],
        scratch: &'s mut DenseDosageScratch,
    ) -> Result<DenseDosageBatch<'s>, CohortError> {
        let n_samples = self.n_samples()?;
        let len = vcfs.len();
        scratch.reserve(n_samples, len);
        let sg = self.sparse_g()?;
        let buf = scratch.live_mut();
        for byte in buf.iter_mut() {
            *byte = 0.0;
        }
        for (i, &vcf) in vcfs.iter().enumerate() {
            let cl = sg.load_variant(vcf.get());
            for entry in &cl.entries {
                let row = entry.sample_idx as usize;
                buf[row * len + i] = entry.dosage as f32;
            }
        }
        Ok(DenseDosageBatch {
            n_samples,
            len,
            data: buf,
        })
    }
}

/// Resolve a `VariantSelector` arm against `view`. Returns the matching
/// vcfs as a sorted, deduplicated `Vec<VariantVcf>`. The `view.select`
/// wrapper just stamps the (cohort, chromosome) on top.
fn resolve_selector(
    view: &ChromosomeView<'_>,
    sel: &VariantSelector,
) -> Result<Vec<VariantVcf>, CohortError> {
    let index = view.index()?;
    match sel {
        VariantSelector::All => Ok((0..index.len() as u32).map(VariantVcf).collect()),
        VariantSelector::Gene(name) => {
            Ok(super::types::from_u32_slice(index.gene_variant_vcfs(name)).to_vec())
        }
        VariantSelector::Region(r) => {
            if r.chromosome != view.chromosome() {
                return Err(CohortError::Input(format!(
                    "region chromosome {} does not match view chromosome {}",
                    r.chromosome.label(),
                    view.chromosome().label()
                )));
            }
            // Inclusive interval overlap on [position, end_position]:
            // a variant is in the region whenever its reference span
            // touches [start, end]. Linear today; the IntervalLookup
            // (`store.lookup::<IntervalLookup>(&cohort, region)`) makes
            // this O(log n + matches) once the index is built.
            let mut out = Vec::new();
            for (i, e) in index.all_entries().iter().enumerate() {
                if e.end_position >= r.start && e.position <= r.end {
                    out.push(VariantVcf(i as u32));
                }
            }
            Ok(out)
        }
        VariantSelector::Vcf(v) => Ok(vec![*v]),
        VariantSelector::Vcfs(vs) => {
            let mut out = vs.clone();
            out.sort_unstable();
            out.dedup();
            Ok(out)
        }
        VariantSelector::Vid(vid) => Ok(index
            .resolve_vid(vid)
            .map(|v| vec![VariantVcf(v)])
            .unwrap_or_default()),
        VariantSelector::Vids(vids) => {
            let mut out: Vec<VariantVcf> = vids
                .iter()
                .filter_map(|v| index.resolve_vid(v).map(VariantVcf))
                .collect();
            out.sort_unstable();
            out.dedup();
            Ok(out)
        }
        VariantSelector::List(_id) => Err(CohortError::DataMissing(
            "VariantSelector::List requires the variant-list subsystem (Phase 5)".into(),
        )),
        VariantSelector::MafBelow(maf) => {
            let mut out = Vec::new();
            for (i, e) in index.all_entries().iter().enumerate() {
                if e.maf < *maf {
                    out.push(VariantVcf(i as u32));
                }
            }
            Ok(out)
        }
        VariantSelector::Consequence(c) => {
            let mut out = Vec::new();
            for (i, e) in index.all_entries().iter().enumerate() {
                if e.consequence == *c {
                    out.push(VariantVcf(i as u32));
                }
            }
            Ok(out)
        }
        VariantSelector::And(l, r) => {
            let lv = resolve_selector(view, l)?;
            let rv = resolve_selector(view, r)?;
            Ok(intersect_sorted(&lv, &rv))
        }
        VariantSelector::Or(l, r) => {
            let lv = resolve_selector(view, l)?;
            let rv = resolve_selector(view, r)?;
            Ok(union_sorted(&lv, &rv))
        }
    }
}


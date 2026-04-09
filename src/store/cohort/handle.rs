//! Lazy cohort and chromosome handles.
//!
//! `CohortHandle::manifest`/`samples`/`chromosomes` materialize on first
//! call via `OnceLock`. `ChromosomeView` mmap's `sparse_g.bin` and parses
//! `variants.parquet` + `membership.parquet` on first method call. Both
//! are `Send + Sync` so the gene-parallel scoring loop can share a view
//! across rayon worker threads.

use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use crate::column;
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::store::Store;
use crate::types::Chromosome;

use super::sparse_g::SparseG;
use super::types::{as_u32_slice, from_u32_slice, CarrierBatch, VariantVcf};
use super::variants::VariantIndex;
use super::{
    build, describe_miss, probe as cohort_probe, read_sample_names_at, run_annotation_join,
    CohortId, CohortManifest, GenoStoreResult, StoreProbe,
};

pub struct CohortSources<'a> {
    pub genotypes: &'a Path,
    pub annotations: &'a Path,
}

pub struct BuildOpts<'a> {
    pub staging_dir: &'a Path,
    pub rebuild: bool,
    pub probe: StoreProbe,
}

pub struct CohortHandle<'a> {
    store: &'a Store,
    id: CohortId,
    dir: PathBuf,
}

impl<'a> CohortHandle<'a> {
    pub(crate) fn new(store: &'a Store, id: CohortId, dir: PathBuf) -> Self {
        Self { store, id, dir }
    }

    pub fn dir(&self) -> &std::path::Path {
        &self.dir
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
            store: self.store,
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

    /// Load a pre-built cohort by id without recomputing the source
    /// fingerprint. Used by `cohort staar --cohort <id>`: the operator
    /// already built the cohort via `cohort ingest`, and we trust the
    /// manifest. Errors with `DataMissing` if the manifest is absent,
    /// unparseable, on a wrong schema version, or missing chromosome
    /// artifacts.
    pub fn load(&self, engine: &Engine) -> Result<GenoStoreResult, CohortError> {
        let manifest_path = self.dir.join("manifest.json");
        if !manifest_path.exists() {
            return Err(CohortError::DataMissing(format!(
                "Cohort '{}' not found at {}. Run `cohort ingest <vcf> --annotations <set> \
                 --cohort-id {}` first.",
                self.id.as_str(),
                self.dir.display(),
                self.id.as_str()
            )));
        }
        let manifest_json = std::fs::read_to_string(&manifest_path).map_err(|e| {
            CohortError::DataMissing(format!(
                "Cannot read cohort manifest '{}': {e}",
                manifest_path.display()
            ))
        })?;
        let manifest: CohortManifest = serde_json::from_str(&manifest_json).map_err(|e| {
            CohortError::DataMissing(format!(
                "Cohort manifest '{}' is unparseable: {e}",
                manifest_path.display()
            ))
        })?;
        if manifest.version != 4 {
            return Err(CohortError::DataMissing(format!(
                "Cohort '{}' is on schema v{}, this build expects v4. Re-run \
                 `cohort ingest <vcf> --annotations <set> --cohort-id {} --rebuild`.",
                self.id.as_str(),
                manifest.version,
                self.id.as_str()
            )));
        }
        for ci in &manifest.chromosomes {
            for artifact in ["sparse_g.bin", "variants.parquet"] {
                let p = self.dir.join(format!("chromosome={}/{}", ci.name, artifact));
                if !p.exists() {
                    return Err(CohortError::DataMissing(format!(
                        "Cohort '{}' chromosome={} missing {}",
                        self.id.as_str(),
                        ci.name,
                        artifact
                    )));
                }
            }
        }
        let sample_names = read_sample_names_at(&self.dir)?;
        let df = DfEngine::new(engine.resources())?;
        Ok(GenoStoreResult {
            geno: crate::staar::genotype::GenotypeResult {
                sample_names,
                output_dir: self.dir.clone(),
            },
            manifest,
            engine: df,
            store_dir: self.dir.clone(),
        })
    }

    pub fn build_or_load(
        &self,
        sources: CohortSources<'_>,
        opts: BuildOpts<'_>,
        engine: &Engine,
        out: &dyn Output,
    ) -> Result<GenoStoreResult, CohortError> {
        if let Some(parent) = self.dir.parent() {
            std::fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", parent.display()))
            })?;
        }

        let CohortSources { genotypes, annotations } = sources;
        let BuildOpts { staging_dir, rebuild, probe: probe_result } = opts;
        let res = engine.resources();

        if let Some(manifest) = probe_result.manifest {
            out.status(&format!(
                "Using cached genotype store ({} variants x {} samples)",
                manifest.n_variants, manifest.n_samples,
            ));

            let sample_names = read_sample_names_at(&probe_result.store_dir)?;
            let df = DfEngine::new(res)?;

            return Ok(GenoStoreResult {
                geno: crate::staar::genotype::GenotypeResult {
                    sample_names,
                    output_dir: probe_result.store_dir.clone(),
                },
                manifest,
                engine: df,
                store_dir: probe_result.store_dir,
            });
        }

        let why = match (&probe_result.miss_reason, rebuild) {
            (_, true) => "  Rebuild requested by --rebuild-store".to_string(),
            (Some(r), _) => format!("  Cache miss: {}", describe_miss(r)),
            (None, _) => "  Cache miss: probe returned no reason".to_string(),
        };
        out.status(&why);
        out.status("Building genotype store...");

        out.status("  Extracting genotypes from VCF...");
        let geno = crate::staar::genotype::extract_genotypes(
            genotypes,
            staging_dir,
            res.memory_bytes,
            res.threads,
            out,
        )?;

        out.status("  Joining genotypes with annotations...");
        let df = run_annotation_join(annotations, &geno, res)?;

        let dup_count = df.query_scalar(&column::dedup_count_sql())?;
        if dup_count > 0 {
            out.warn(&format!(
                "  {dup_count} duplicate variants found — keeping first occurrence."
            ));
            df.execute(&column::dedup_sql())?;
            df.execute("DROP TABLE _rare_all")?;
            df.execute("ALTER TABLE _rare_dedup RENAME TO _rare_all")?;
        }

        let n_all = df.query_scalar("SELECT COUNT(*) FROM _rare_all")?;
        let n_genes = df.query_scalar(&column::gene_count_sql())?;
        out.status(&format!(
            "  {} annotated variants, {} genes",
            n_all, n_genes
        ));

        if n_all == 0 {
            return Err(CohortError::Analysis(format!(
                "No variants found after joining genotypes ({}) with annotations ({}). \
                 Check that both use the same genome build and allele normalization.",
                genotypes.display(),
                annotations.display(),
            )));
        }

        let manifest = build(&df, &geno, &self.dir, genotypes, annotations, out)?;

        Ok(GenoStoreResult {
            geno,
            manifest,
            engine: df,
            store_dir: self.dir.clone(),
        })
    }
}

/// Borrowed view of one (cohort, chromosome). Methods open
/// `sparse_g.bin`, `variants.parquet`, and `membership.parquet` on the
/// first access that needs them. Both inner caches are `OnceLock` so
/// rayon parallel iteration over genes is safe.
pub struct ChromosomeView<'a> {
    store: &'a Store,
    chrom: Chromosome,
    dir: PathBuf,
    sparse_g: OnceLock<SparseG>,
    index: OnceLock<VariantIndex>,
}

impl<'a> ChromosomeView<'a> {
    pub fn chromosome(&self) -> Chromosome {
        self.chrom
    }

    /// Internal: mmap'd sparse genotype matrix.
    pub fn sparse_g(&self) -> Result<&SparseG, CohortError> {
        if let Some(sg) = self.sparse_g.get() {
            return Ok(sg);
        }
        let sg = SparseG::open(self.store.backend(), &self.dir)?;
        Ok(self.sparse_g.get_or_init(|| sg))
    }

    /// Internal: parsed variant metadata + gene membership.
    pub fn index(&self) -> Result<&VariantIndex, CohortError> {
        if let Some(idx) = self.index.get() {
            return Ok(idx);
        }
        let idx = VariantIndex::load(self.store.backend(), &self.dir)?;
        Ok(self.index.get_or_init(|| idx))
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

    /// Variants assigned to `gene` on this chromosome, in dense
    /// `VariantVcf` order. Empty when the gene is not in the membership
    /// table.
    pub fn gene_variants(&self, gene: &str) -> Result<Vec<VariantVcf>, CohortError> {
        Ok(from_u32_slice(self.index()?.gene_variant_vcfs(gene)).to_vec())
    }
}

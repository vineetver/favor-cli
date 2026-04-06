//! Variant metadata index and carrier data types.
//!
//! VariantIndex loads aligned metadata vectors from variants.parquet and
//! gene membership from membership.parquet. All data is indexed by variant_vcf
//! — a dense, immutable index assigned at store build time.
//!
//! CarrierList and CarrierEntry are the canonical sparse genotype types
//! consumed by scoring kernels.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use crate::error::FavorError;
use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, Consequence, FunctionalAnnotation,
    RegionType, RegulatoryFlags,
};

// ═══════════════════════════════════════════════════════════════════════
// Carrier types — format-agnostic, consumed by sparse_score.rs
// ═══════════════════════════════════════════════════════════════════════

/// A single carrier entry: (sample_index, dosage).
#[derive(Clone, Copy, Debug)]
pub struct CarrierEntry {
    pub sample_idx: u32,
    pub dosage: u8,
}

/// Carrier list for one variant — all non-reference samples.
#[derive(Clone, Debug)]
pub struct CarrierList {
    pub entries: Vec<CarrierEntry>,
}

impl CarrierList {
    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn mac(&self) -> u32 {
        self.entries.iter().map(|e| e.dosage as u32).sum()
    }
}

// ═══════════════════════════════════════════════════════════════════════
// WeightVector — named, composable annotation weight channel
// ═══════════════════════════════════════════════════════════════════════

/// A named weight vector aligned to variant_vcf.
/// First-class citizen: can be selected, combined, and composed at runtime.
#[allow(dead_code)] // fields accessed via VariantIndex weight methods
pub struct WeightVector {
    pub name: &'static str,
    pub values: Vec<f64>,
}

// ═══════════════════════════════════════════════════════════════════════
// VariantIndex — aligned metadata + mask resolution
// ═══════════════════════════════════════════════════════════════════════

/// One variant's metadata as loaded from variants.parquet.
#[derive(Clone, Debug)]
pub struct VariantIndexEntry {
    pub position: u32,
    pub ref_allele: Box<str>,
    pub alt_allele: Box<str>,
    pub vid: Box<str>,
    pub maf: f64,
    pub region_type: RegionType,
    pub consequence: Consequence,
    pub cadd_phred: f64,
    pub revel: f64,
    pub regulatory: RegulatoryFlags,
    pub weights: AnnotationWeights,
}

impl VariantIndexEntry {
    /// Convert to the canonical AnnotatedVariant type.
    pub fn to_annotated_variant(&self, chrom: Chromosome) -> AnnotatedVariant {
        AnnotatedVariant {
            chromosome: chrom,
            position: self.position,
            ref_allele: self.ref_allele.clone(),
            alt_allele: self.alt_allele.clone(),
            maf: self.maf,
            gene_name: "".into(), // gene is a derived property, not stored per-variant
            annotation: FunctionalAnnotation {
                region_type: self.region_type,
                consequence: self.consequence,
                cadd_phred: self.cadd_phred,
                revel: self.revel,
                regulatory: self.regulatory,
                weights: self.weights,
            },
        }
    }

    /// Convert with an explicit gene name (for mask predicates that check gene_name).
    pub fn to_annotated_variant_with_gene(
        &self,
        chrom: Chromosome,
        gene_name: &str,
    ) -> AnnotatedVariant {
        let mut av = self.to_annotated_variant(chrom);
        av.gene_name = gene_name.into();
        av
    }
}

/// All variant metadata for one chromosome, indexed by variant_vcf.
/// Provides aligned vectors and on-demand mask compilation.
pub struct VariantIndex {
    /// Variant metadata indexed by variant_vcf. entries[i] = variant_vcf i.
    entries: Vec<VariantIndexEntry>,
    /// gene_name → sorted Vec<u32> of variant_vcfs belonging to this gene.
    /// Loaded from membership.parquet.
    gene_variants: HashMap<String, Vec<u32>>,
    /// vid string → variant_vcf.
    vid_to_vcf: HashMap<Box<str>, u32>,
    /// Named weight vectors, each of length n_variants, aligned to variant_vcf.
    #[allow(dead_code)]
    pub weight_channels: Vec<WeightVector>,
}

impl VariantIndex {
    /// Load from variants.parquet + membership.parquet in a chromosome directory.
    pub fn load(chrom_dir: &Path) -> Result<Self, FavorError> {
        let entries = load_variant_entries(chrom_dir)?;
        let gene_variants = load_membership(chrom_dir)?;

        let vid_to_vcf: HashMap<Box<str>, u32> = entries
            .iter()
            .enumerate()
            .map(|(i, e)| (e.vid.clone(), i as u32))
            .collect();

        // Build named weight channels from entries
        let weight_names: [&str; 11] = [
            "cadd_phred", "linsight", "fathmm_xf",
            "apc_epigenetics_active", "apc_epigenetics_repressed",
            "apc_epigenetics_transcription", "apc_conservation",
            "apc_protein_function", "apc_local_nucleotide_diversity",
            "apc_mutation_density", "apc_transcription_factor",
        ];
        let weight_channels: Vec<WeightVector> = (0..11)
            .map(|ch| WeightVector {
                name: weight_names[ch],
                values: entries.iter().map(|e| e.weights.0[ch]).collect(),
            })
            .collect();

        Ok(Self { entries, gene_variants, vid_to_vcf, weight_channels })
    }

    // ── Mask resolution (compiled on demand) ─────────────────────────────

    /// Gene → variant_vcf set (sorted). Empty slice if gene not found.
    pub fn gene_variant_vcfs(&self, gene: &str) -> &[u32] {
        self.gene_variants.get(gene).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// All gene names.
    pub fn gene_names(&self) -> impl Iterator<Item = &str> {
        self.gene_variants.keys().map(|s| s.as_str())
    }

    /// Number of genes.
    pub fn n_genes(&self) -> usize {
        self.gene_variants.len()
    }

    /// Compile a mask on demand: applies MAF + predicate to a gene's variants.
    /// Returns qualifying variant_vcfs as LOCAL indices into gene_vcfs.
    /// No pre-materialization — memory stays O(N) for base metadata.
    #[allow(dead_code)]
    pub fn compile_mask(
        &self,
        gene_vcfs: &[u32],
        maf_cutoff: f64,
        chrom: Chromosome,
        gene_name: &str,
        predicate: fn(&AnnotatedVariant) -> bool,
    ) -> Vec<usize> {
        gene_vcfs
            .iter()
            .enumerate()
            .filter(|(_, &v)| {
                let e = &self.entries[v as usize];
                e.maf < maf_cutoff
                    && predicate(&e.to_annotated_variant_with_gene(chrom, gene_name))
            })
            .map(|(i, _)| i)
            .collect()
    }

    // ── Direct access ────────────────────────────────────────────────────

    /// Get variant metadata by variant_vcf (O(1) direct index).
    #[inline]
    pub fn get(&self, variant_vcf: u32) -> &VariantIndexEntry {
        &self.entries[variant_vcf as usize]
    }

    /// Resolve vid string → variant_vcf.
    pub fn resolve_vid(&self, vid: &str) -> Option<u32> {
        self.vid_to_vcf.get(vid).copied()
    }

    /// Total number of variant entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// All entries (for sliding window / SCANG that operate across genes).
    pub fn all_entries(&self) -> &[VariantIndexEntry] {
        &self.entries
    }

    // ── Weight channel access ────────────────────────────────────────────

    /// Get a weight channel by index.
    #[allow(dead_code)]
    pub fn weight_channel(&self, index: usize) -> &WeightVector {
        &self.weight_channels[index]
    }

    /// Extract weight values for a subset of variant_vcfs.
    #[allow(dead_code)]
    pub fn weight_values(&self, channel: usize, variant_vcfs: &[u32]) -> Vec<f64> {
        let vals = &self.weight_channels[channel].values;
        variant_vcfs.iter().map(|&v| vals[v as usize]).collect()
    }
}

// ═══════════════════════════════════════════════════════════════════════
// Loaders
// ═══════════════════════════════════════════════════════════════════════

fn load_variant_entries(chrom_dir: &Path) -> Result<Vec<VariantIndexEntry>, FavorError> {
    use arrow::array::{BooleanArray, Float64Array, Int32Array, StringArray, UInt32Array};

    let pq_path = chrom_dir.join("variants.parquet");
    let file = File::open(&pq_path)
        .map_err(|e| FavorError::Resource(format!("Open {}: {e}", pq_path.display())))?;
    let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| FavorError::Resource(format!("Parquet open: {e}")))?
        .build()
        .map_err(|e| FavorError::Resource(format!("Parquet reader: {e}")))?;

    let mut entries = Vec::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| FavorError::Resource(format!("Read: {e}")))?;
        let n = batch.num_rows();

        // Schema: variant_vcf(0), position(1), ref(2), alt(3), vid(4), maf(5),
        //         region_type(6), consequence(7), cadd(8), revel(9),
        //         cage_prom(10), cage_enh(11), ccre_prom(12), ccre_enh(13),
        //         w_cadd(14)..w_apc_tf(24)
        let _vvcf_arr = batch.column(0).as_any().downcast_ref::<UInt32Array>().unwrap();
        let pos_arr = batch.column(1).as_any().downcast_ref::<Int32Array>().unwrap();
        let ref_arr = batch.column(2).as_any().downcast_ref::<StringArray>().unwrap();
        let alt_arr = batch.column(3).as_any().downcast_ref::<StringArray>().unwrap();
        let vid_arr = batch.column(4).as_any().downcast_ref::<StringArray>().unwrap();
        let maf_arr = batch.column(5).as_any().downcast_ref::<Float64Array>().unwrap();
        let rt_arr = batch.column(6).as_any().downcast_ref::<StringArray>().unwrap();
        let csq_arr = batch.column(7).as_any().downcast_ref::<StringArray>().unwrap();
        let cadd_arr = batch.column(8).as_any().downcast_ref::<Float64Array>().unwrap();
        let revel_arr = batch.column(9).as_any().downcast_ref::<Float64Array>().unwrap();
        let cp_arr = batch.column(10).as_any().downcast_ref::<BooleanArray>().unwrap();
        let ce_arr = batch.column(11).as_any().downcast_ref::<BooleanArray>().unwrap();
        let crp_arr = batch.column(12).as_any().downcast_ref::<BooleanArray>().unwrap();
        let cre_arr = batch.column(13).as_any().downcast_ref::<BooleanArray>().unwrap();

        let mut w_arrs: Vec<&Float64Array> = Vec::with_capacity(11);
        for i in 0..11 {
            w_arrs.push(batch.column(14 + i).as_any().downcast_ref::<Float64Array>().unwrap());
        }

        for i in 0..n {
            let mut weights = [0.0f64; 11];
            for (ch, wa) in w_arrs.iter().enumerate() {
                let w = wa.value(i);
                weights[ch] = if w.is_finite() { w } else { 0.0 };
            }

            entries.push(VariantIndexEntry {
                position: pos_arr.value(i) as u32,
                ref_allele: ref_arr.value(i).into(),
                alt_allele: alt_arr.value(i).into(),
                vid: vid_arr.value(i).into(),
                maf: maf_arr.value(i),
                region_type: RegionType::from_str_lossy(rt_arr.value(i)),
                consequence: Consequence::from_str_lossy(csq_arr.value(i)),
                cadd_phred: cadd_arr.value(i),
                revel: revel_arr.value(i),
                regulatory: RegulatoryFlags {
                    cage_promoter: cp_arr.value(i),
                    cage_enhancer: ce_arr.value(i),
                    ccre_promoter: crp_arr.value(i),
                    ccre_enhancer: cre_arr.value(i),
                },
                weights: AnnotationWeights(weights),
            });
        }
    }

    Ok(entries)
}

fn load_membership(chrom_dir: &Path) -> Result<HashMap<String, Vec<u32>>, FavorError> {
    use arrow::array::{StringArray, UInt32Array};

    let pq_path = chrom_dir.join("membership.parquet");
    let file = File::open(&pq_path)
        .map_err(|e| FavorError::Resource(format!("Open {}: {e}", pq_path.display())))?;
    let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| FavorError::Resource(format!("Parquet open: {e}")))?
        .build()
        .map_err(|e| FavorError::Resource(format!("Parquet reader: {e}")))?;

    let mut gene_variants: HashMap<String, Vec<u32>> = HashMap::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| FavorError::Resource(format!("Read: {e}")))?;
        let vvcf_arr = batch.column(0).as_any().downcast_ref::<UInt32Array>().unwrap();
        let gene_arr = batch.column(1).as_any().downcast_ref::<StringArray>().unwrap();

        for i in 0..batch.num_rows() {
            let gene = gene_arr.value(i).to_string();
            let vvcf = vvcf_arr.value(i);
            gene_variants.entry(gene).or_default().push(vvcf);
        }
    }

    // Ensure variant_vcfs are sorted per gene
    for variants in gene_variants.values_mut() {
        variants.sort_unstable();
        variants.dedup();
    }

    Ok(gene_variants)
}

// ═══════════════════════════════════════════════════════════════════════
// Tests — alignment invariants
// ═══════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::AnnotationWeights;

    /// Helper: build a minimal VariantIndex in memory for testing.
    fn make_test_index() -> VariantIndex {
        let entries = vec![
            VariantIndexEntry {
                position: 100,
                ref_allele: "A".into(),
                alt_allele: "T".into(),
                vid: "22-100-A-T".into(),
                maf: 0.001,
                region_type: RegionType::Exonic,
                consequence: Consequence::MissenseVariant,
                cadd_phred: 25.0,
                revel: 0.8,
                regulatory: RegulatoryFlags::default(),
                weights: AnnotationWeights([1.0; 11]),
            },
            VariantIndexEntry {
                position: 200,
                ref_allele: "C".into(),
                alt_allele: "G".into(),
                vid: "22-200-C-G".into(),
                maf: 0.005,
                region_type: RegionType::Exonic,
                consequence: Consequence::Stopgain,
                cadd_phred: 35.0,
                revel: 0.95,
                regulatory: RegulatoryFlags::default(),
                weights: AnnotationWeights([0.5; 11]),
            },
            VariantIndexEntry {
                position: 300,
                ref_allele: "G".into(),
                alt_allele: "A".into(),
                vid: "22-300-G-A".into(),
                maf: 0.002,
                region_type: RegionType::Intronic,
                consequence: Consequence::Unknown,
                cadd_phred: 5.0,
                revel: 0.1,
                regulatory: RegulatoryFlags::default(),
                weights: AnnotationWeights([0.2; 11]),
            },
        ];

        let mut gene_variants = HashMap::new();
        gene_variants.insert("GENE1".to_string(), vec![0, 1]);
        gene_variants.insert("GENE2".to_string(), vec![1, 2]);

        let vid_to_vcf: HashMap<Box<str>, u32> = entries
            .iter()
            .enumerate()
            .map(|(i, e)| (e.vid.clone(), i as u32))
            .collect();

        let weight_names: [&str; 11] = [
            "cadd_phred", "linsight", "fathmm_xf",
            "apc_epigenetics_active", "apc_epigenetics_repressed",
            "apc_epigenetics_transcription", "apc_conservation",
            "apc_protein_function", "apc_local_nucleotide_diversity",
            "apc_mutation_density", "apc_transcription_factor",
        ];
        let weight_channels: Vec<WeightVector> = (0..11)
            .map(|ch| WeightVector {
                name: weight_names[ch],
                values: entries.iter().map(|e| e.weights.0[ch]).collect(),
            })
            .collect();

        VariantIndex { entries, gene_variants, vid_to_vcf, weight_channels }
    }

    #[test]
    fn gene_variant_vcfs_returns_correct_set() {
        let vi = make_test_index();
        assert_eq!(vi.gene_variant_vcfs("GENE1"), &[0, 1]);
        assert_eq!(vi.gene_variant_vcfs("GENE2"), &[1, 2]);
        assert_eq!(vi.gene_variant_vcfs("NONEXISTENT"), &[] as &[u32]);
    }

    #[test]
    fn overlapping_genes_share_variant() {
        let vi = make_test_index();
        let g1 = vi.gene_variant_vcfs("GENE1");
        let g2 = vi.gene_variant_vcfs("GENE2");
        // variant_vcf 1 belongs to both genes
        assert!(g1.contains(&1));
        assert!(g2.contains(&1));
    }

    #[test]
    fn resolve_vid_roundtrip() {
        let vi = make_test_index();
        for (i, e) in vi.all_entries().iter().enumerate() {
            let resolved = vi.resolve_vid(&e.vid);
            assert_eq!(resolved, Some(i as u32), "vid {} resolved to {:?}", e.vid, resolved);
        }
        assert_eq!(vi.resolve_vid("nonexistent"), None);
    }

    #[test]
    fn get_by_variant_vcf() {
        let vi = make_test_index();
        assert_eq!(vi.get(0).position, 100);
        assert_eq!(vi.get(1).position, 200);
        assert_eq!(vi.get(2).position, 300);
    }

    #[test]
    fn weight_channel_alignment() {
        let vi = make_test_index();
        for ch in 0..11 {
            let wv = vi.weight_channel(ch);
            assert_eq!(wv.values.len(), vi.len());
            for (i, &val) in wv.values.iter().enumerate() {
                assert_eq!(val, vi.get(i as u32).weights.0[ch],
                    "weight channel {ch} misaligned at variant_vcf {i}");
            }
        }
    }

    #[test]
    fn weight_values_subset() {
        let vi = make_test_index();
        let subset = &[0u32, 2];
        let vals = vi.weight_values(0, subset);
        assert_eq!(vals.len(), 2);
        assert_eq!(vals[0], vi.get(0).weights.0[0]);
        assert_eq!(vals[1], vi.get(2).weights.0[0]);
    }

    #[test]
    fn compile_mask_filters_correctly() {
        let vi = make_test_index();
        let gene_vcfs = vi.gene_variant_vcfs("GENE1");
        let chrom = Chromosome::Autosome(22);

        // All pass (maf_cutoff = 1.0, always-true predicate)
        let all = vi.compile_mask(gene_vcfs, 1.0, chrom, "GENE1", |_| true);
        assert_eq!(all, vec![0, 1]);

        // MAF filter: only variant 0 (maf=0.001) passes cutoff 0.002
        let maf_filtered = vi.compile_mask(gene_vcfs, 0.002, chrom, "GENE1", |_| true);
        assert_eq!(maf_filtered, vec![0]);

        // Predicate filter: only missense
        let missense_only = vi.compile_mask(gene_vcfs, 1.0, chrom, "GENE1", |av| {
            av.annotation.consequence.is_missense()
        });
        assert_eq!(missense_only, vec![0]); // variant 0 is MissenseVariant

        // Combined: maf < 0.01 AND is_plof
        let plof = vi.compile_mask(gene_vcfs, 0.01, chrom, "GENE1", |av| {
            av.annotation.consequence.is_plof()
        });
        assert_eq!(plof, vec![1]); // variant 1 is Stopgain (pLoF)
    }

    #[test]
    fn mask_composition_equivalence() {
        let vi = make_test_index();
        let gene_vcfs = vi.gene_variant_vcfs("GENE1");
        let chrom = Chromosome::Autosome(22);

        // Path A: compile_mask with combined predicate
        let combined = vi.compile_mask(gene_vcfs, 0.01, chrom, "GENE1", |av| {
            av.annotation.consequence.is_missense() || av.annotation.consequence.is_plof()
        });

        // Path B: union of two separate masks
        let missense = vi.compile_mask(gene_vcfs, 0.01, chrom, "GENE1", |av| {
            av.annotation.consequence.is_missense()
        });
        let plof = vi.compile_mask(gene_vcfs, 0.01, chrom, "GENE1", |av| {
            av.annotation.consequence.is_plof()
        });
        let mut union: Vec<usize> = missense.iter().chain(plof.iter()).copied().collect();
        union.sort_unstable();
        union.dedup();

        assert_eq!(combined, union);
    }

    #[test]
    fn sortedness_invariant() {
        let vi = make_test_index();
        let entries = vi.all_entries();
        for i in 1..entries.len() {
            let prev = &entries[i - 1];
            let curr = &entries[i];
            assert!(
                (curr.position, &*curr.ref_allele, &*curr.alt_allele)
                    >= (prev.position, &*prev.ref_allele, &*prev.alt_allele),
                "entries not sorted at index {i}"
            );
        }
        for gene in vi.gene_names() {
            let vcfs = vi.gene_variant_vcfs(gene);
            for i in 1..vcfs.len() {
                assert!(vcfs[i] > vcfs[i - 1],
                    "gene {gene} membership not sorted at index {i}");
            }
        }
    }

    #[test]
    fn to_annotated_variant_preserves_fields() {
        let vi = make_test_index();
        let chrom = Chromosome::Autosome(22);
        let av = vi.get(0).to_annotated_variant_with_gene(chrom, "GENE1");
        assert_eq!(av.chromosome, chrom);
        assert_eq!(av.position, 100);
        assert_eq!(&*av.ref_allele, "A");
        assert_eq!(&*av.alt_allele, "T");
        assert_eq!(av.maf, 0.001);
        assert_eq!(&*av.gene_name, "GENE1");
        assert_eq!(av.annotation.consequence, Consequence::MissenseVariant);
    }

    /// Region query == manual scan over entries[].
    /// Validates "region = mask over variant_vcf" assumption.
    #[test]
    fn region_mask_equivalence() {
        let vi = make_test_index();

        // Region [150, 250] should match variant_vcf 1 (position 200)
        let start = 150u32;
        let end = 250u32;

        // Path A: manual scan
        let manual: Vec<u32> = vi.all_entries().iter().enumerate()
            .filter(|(_, e)| e.position >= start && e.position <= end)
            .map(|(i, _)| i as u32)
            .collect();

        // Path B: this is what any region query resolves to
        let region_mask: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| {
                let pos = vi.get(v).position;
                pos >= start && pos <= end
            })
            .collect();

        assert_eq!(manual, region_mask);
        assert_eq!(manual, vec![1u32]);

        // Wider region [0, 500] should match all
        let all: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| vi.get(v).position <= 500)
            .collect();
        assert_eq!(all, vec![0, 1, 2]);

        // Empty region [400, 500] should match none
        let empty: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| {
                let pos = vi.get(v).position;
                (400..=500).contains(&pos)
            })
            .collect();
        assert!(empty.is_empty());
    }

    /// Cross-layer: vid computed from metadata matches stored vid.
    #[test]
    fn vid_from_metadata_matches_stored() {
        let vi = make_test_index();
        for v in 0..vi.len() as u32 {
            let e = vi.get(v);
            let computed = crate::types::format_vid("22", e.position, &e.ref_allele, &e.alt_allele);
            assert_eq!(&*e.vid, &*computed,
                "variant_vcf {v}: stored vid '{}' != computed '{computed}'", e.vid);
        }
    }
}

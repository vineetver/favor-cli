//! Variant metadata index and carrier data types.

use std::collections::HashMap;
use std::path::Path;

use crate::error::CohortError;
use crate::store::backend::{Backend, BoxedBatchReader};
use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, Consequence, FunctionalAnnotation, RegionType,
    RegulatoryFlags,
};

#[derive(Clone, Copy, Debug)]
pub struct CarrierEntry {
    pub sample_idx: u32,
    pub dosage: u8,
}

#[derive(Clone, Debug)]
pub struct CarrierList {
    pub entries: Vec<CarrierEntry>,
}

impl CarrierList {
    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }
}

#[derive(Clone, Debug)]
pub struct VariantIndexEntry {
    pub position: u32,
    /// Stored explicitly so interval queries do not need to recompute it from `ref_allele`.
    #[allow(dead_code)] // populated and round-trip-tested; reader added in a later phase
    pub end_position: u32,
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

pub struct VariantIndex {
    entries: Vec<VariantIndexEntry>,
    gene_variants: HashMap<String, Vec<u32>>,
    vid_to_vcf: HashMap<Box<str>, u32>,
}

impl VariantIndex {
    pub fn load(backend: &dyn Backend, chrom_dir: &Path) -> Result<Self, CohortError> {
        let variants_reader = backend.open_parquet(&chrom_dir.join("variants.parquet"))?;
        let entries = load_variant_entries(variants_reader)?;
        let gene_variants = super::membership::load(backend, chrom_dir)?;

        let vid_to_vcf: HashMap<Box<str>, u32> = entries
            .iter()
            .enumerate()
            .map(|(i, e)| (e.vid.clone(), i as u32))
            .collect();

        Ok(Self {
            entries,
            gene_variants,
            vid_to_vcf,
        })
    }

    pub fn gene_variant_vcfs(&self, gene: &str) -> &[u32] {
        self.gene_variants
            .get(gene)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn gene_names(&self) -> impl Iterator<Item = &str> {
        self.gene_variants.keys().map(|s| s.as_str())
    }

    pub fn n_genes(&self) -> usize {
        self.gene_variants.len()
    }

    #[inline]
    pub fn get(&self, variant_vcf: u32) -> &VariantIndexEntry {
        &self.entries[variant_vcf as usize]
    }

    pub fn resolve_vid(&self, vid: &str) -> Option<u32> {
        self.vid_to_vcf.get(vid).copied()
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn all_entries(&self) -> &[VariantIndexEntry] {
        &self.entries
    }
}

/// Look up a typed Arrow array by column name so schema reorderings stay safe.
/// Shared by the variants.parquet reader here and the builder that writes it.
pub(super) fn col_by_name<'a, T: arrow::array::Array + 'static>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &str,
) -> Result<&'a T, CohortError> {
    let idx = batch
        .schema()
        .index_of(name)
        .map_err(|_| CohortError::Resource(format!("variants.parquet missing column {name}")))?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<T>()
        .ok_or_else(|| {
            CohortError::Resource(format!("variants.parquet column {name} has wrong type"))
        })
}

fn load_variant_entries(reader: BoxedBatchReader) -> Result<Vec<VariantIndexEntry>, CohortError> {
    use arrow::array::{BooleanArray, Float64Array, Int32Array, StringArray};

    use crate::column::{Col, STAAR_WEIGHTS};

    let mut entries = Vec::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| CohortError::Resource(format!("Read: {e}")))?;
        let n = batch.num_rows();

        let pos_arr = col_by_name::<Int32Array>(&batch, Col::Position.as_str())?;
        let end_arr = col_by_name::<Int32Array>(&batch, Col::EndPosition.as_str())?;
        let ref_arr = col_by_name::<StringArray>(&batch, Col::RefAllele.as_str())?;
        let alt_arr = col_by_name::<StringArray>(&batch, Col::AltAllele.as_str())?;
        let vid_arr = col_by_name::<StringArray>(&batch, Col::Vid.as_str())?;
        let maf_arr = col_by_name::<Float64Array>(&batch, Col::Maf.as_str())?;
        let rt_arr = col_by_name::<StringArray>(&batch, Col::RegionType.as_str())?;
        let csq_arr = col_by_name::<StringArray>(&batch, Col::Consequence.as_str())?;
        let cadd_arr = col_by_name::<Float64Array>(&batch, Col::CaddPhred.as_str())?;
        let revel_arr = col_by_name::<Float64Array>(&batch, Col::Revel.as_str())?;
        let cp_arr = col_by_name::<BooleanArray>(&batch, Col::IsCagePromoter.as_str())?;
        let ce_arr = col_by_name::<BooleanArray>(&batch, Col::IsCageEnhancer.as_str())?;
        let crp_arr = col_by_name::<BooleanArray>(&batch, Col::IsCcrePromoter.as_str())?;
        let cre_arr = col_by_name::<BooleanArray>(&batch, Col::IsCcreEnhancer.as_str())?;

        let w_arrs: Vec<&Float64Array> = STAAR_WEIGHTS
            .iter()
            .map(|c| col_by_name::<Float64Array>(&batch, c.as_str()))
            .collect::<Result<Vec<_>, _>>()?;

        for i in 0..n {
            let mut weights = [0.0f64; 11];
            for (ch, wa) in w_arrs.iter().enumerate() {
                let w = wa.value(i);
                weights[ch] = if w.is_finite() { w } else { 0.0 };
            }

            entries.push(VariantIndexEntry {
                position: pos_arr.value(i) as u32,
                end_position: end_arr.value(i) as u32,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::AnnotationWeights;

    /// Helper: build a minimal VariantIndex in memory for testing.
    fn make_test_index() -> VariantIndex {
        let entries = vec![
            VariantIndexEntry {
                position: 100,
                end_position: 100,
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
                end_position: 200,
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
                end_position: 300,
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

        VariantIndex {
            entries,
            gene_variants,
            vid_to_vcf,
        }
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
        assert!(g1.contains(&1));
        assert!(g2.contains(&1));
    }

    #[test]
    fn resolve_vid_roundtrip() {
        let vi = make_test_index();
        for (i, e) in vi.all_entries().iter().enumerate() {
            let resolved = vi.resolve_vid(&e.vid);
            assert_eq!(
                resolved,
                Some(i as u32),
                "vid {} resolved to {:?}",
                e.vid,
                resolved
            );
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
                assert!(
                    vcfs[i] > vcfs[i - 1],
                    "gene {gene} membership not sorted at index {i}"
                );
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

    #[test]
    fn region_mask_equivalence() {
        let vi = make_test_index();

        let start = 150u32;
        let end = 250u32;

        let manual: Vec<u32> = vi
            .all_entries()
            .iter()
            .enumerate()
            .filter(|(_, e)| e.position >= start && e.position <= end)
            .map(|(i, _)| i as u32)
            .collect();

        let region_mask: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| {
                let pos = vi.get(v).position;
                pos >= start && pos <= end
            })
            .collect();

        assert_eq!(manual, region_mask);
        assert_eq!(manual, vec![1u32]);

        let all: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| vi.get(v).position <= 500)
            .collect();
        assert_eq!(all, vec![0, 1, 2]);

        let empty: Vec<u32> = (0..vi.len() as u32)
            .filter(|&v| {
                let pos = vi.get(v).position;
                (400..=500).contains(&pos)
            })
            .collect();
        assert!(empty.is_empty());
    }

    #[test]
    fn vid_from_metadata_matches_stored() {
        let vi = make_test_index();
        for v in 0..vi.len() as u32 {
            let e = vi.get(v);
            let computed = crate::types::format_vid("22", e.position, &e.ref_allele, &e.alt_allele);
            assert_eq!(
                &*e.vid, &*computed,
                "variant_vcf {v}: stored vid '{}' != computed '{computed}'",
                e.vid
            );
        }
    }

    #[test]
    fn variants_parquet_schema_round_trip() {
        use std::sync::Arc;

        use arrow::array::{
            BooleanArray, Float64Array, Int32Array, RecordBatchIterator, StringArray, UInt32Array,
        };
        use arrow::datatypes::{DataType, Field, Schema};
        use arrow::record_batch::RecordBatch;

        use crate::column::{Col, STAAR_WEIGHTS};

        let mut fields = vec![
            Field::new(Col::VariantVcf.as_str(), DataType::UInt32, false),
            Field::new(Col::Position.as_str(), DataType::Int32, false),
            Field::new(Col::EndPosition.as_str(), DataType::Int32, false),
            Field::new(Col::RefAllele.as_str(), DataType::Utf8, false),
            Field::new(Col::AltAllele.as_str(), DataType::Utf8, false),
            Field::new(Col::Vid.as_str(), DataType::Utf8, false),
            Field::new(Col::Maf.as_str(), DataType::Float64, false),
            Field::new(Col::RegionType.as_str(), DataType::Utf8, false),
            Field::new(Col::Consequence.as_str(), DataType::Utf8, false),
            Field::new(Col::CaddPhred.as_str(), DataType::Float64, false),
            Field::new(Col::Revel.as_str(), DataType::Float64, false),
            Field::new(Col::IsCagePromoter.as_str(), DataType::Boolean, false),
            Field::new(Col::IsCageEnhancer.as_str(), DataType::Boolean, false),
            Field::new(Col::IsCcrePromoter.as_str(), DataType::Boolean, false),
            Field::new(Col::IsCcreEnhancer.as_str(), DataType::Boolean, false),
        ];
        for col in &STAAR_WEIGHTS {
            fields.push(Field::new(col.as_str(), DataType::Float64, false));
        }
        let schema = Arc::new(Schema::new(fields));

        let mut columns: Vec<arrow::array::ArrayRef> = vec![
            Arc::new(UInt32Array::from(vec![0u32, 1u32])),
            Arc::new(Int32Array::from(vec![100i32, 200i32])),
            Arc::new(Int32Array::from(vec![100i32, 203i32])),
            Arc::new(StringArray::from(vec!["A", "ACGT"])),
            Arc::new(StringArray::from(vec!["T", "A"])),
            Arc::new(StringArray::from(vec!["22-100-A-T", "22-200-ACGT-A"])),
            Arc::new(Float64Array::from(vec![0.001, 0.005])),
            Arc::new(StringArray::from(vec!["exonic", "intronic"])),
            Arc::new(StringArray::from(vec!["missense_variant", ""])),
            Arc::new(Float64Array::from(vec![25.0, 5.0])),
            Arc::new(Float64Array::from(vec![0.8, 0.1])),
            Arc::new(BooleanArray::from(vec![false, false])),
            Arc::new(BooleanArray::from(vec![false, true])),
            Arc::new(BooleanArray::from(vec![true, false])),
            Arc::new(BooleanArray::from(vec![false, false])),
        ];
        for ch in 0..11 {
            let v0 = 0.10 + 0.01 * ch as f64;
            let v1 = 0.50 + 0.01 * ch as f64;
            columns.push(Arc::new(Float64Array::from(vec![v0, v1])));
        }

        let batch = RecordBatch::try_new(schema.clone(), columns).unwrap();
        let reader: BoxedBatchReader =
            Box::new(RecordBatchIterator::new(vec![Ok(batch)].into_iter(), schema));
        let entries = load_variant_entries(reader).unwrap();

        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].position, 100);
        assert_eq!(entries[0].end_position, 100);
        assert_eq!(&*entries[0].ref_allele, "A");
        assert_eq!(&*entries[0].alt_allele, "T");
        assert_eq!(&*entries[0].vid, "22-100-A-T");
        assert_eq!(entries[0].maf, 0.001);
        assert_eq!(entries[0].region_type, RegionType::Exonic);
        assert_eq!(entries[0].consequence, Consequence::MissenseVariant);
        assert_eq!(entries[0].cadd_phred, 25.0);
        assert!(entries[0].regulatory.ccre_promoter);
        assert!(!entries[0].regulatory.cage_enhancer);
        for ch in 0..11 {
            assert!((entries[0].weights.0[ch] - (0.10 + 0.01 * ch as f64)).abs() < 1e-12);
        }

        assert_eq!(entries[1].position, 200);
        assert_eq!(entries[1].end_position, 203);
        assert_eq!(&*entries[1].ref_allele, "ACGT");
        assert!(entries[1].regulatory.cage_enhancer);
    }
}

use super::MaskType;
use crate::types::{AnnotatedVariant, Chromosome, Consequence};

type MaskDef = (MaskType, fn(&AnnotatedVariant) -> bool);

#[derive(Debug)]
pub struct MaskGroup {
    pub name: String,
    pub chromosome: Chromosome,
    pub start: u32,
    pub end: u32,
    pub variant_indices: Vec<usize>,
}

pub const CODING_MASKS: &[MaskDef] = &[
    (MaskType::PLof, is_plof),
    (MaskType::Missense, is_missense),
    (MaskType::DisruptiveMissense, is_disruptive_missense),
    (MaskType::PLofMissense, is_plof_or_missense),
    (MaskType::Synonymous, is_synonymous),
    (MaskType::Ptv, is_ptv),
    (MaskType::PtvDs, is_ptv_ds),
];

pub const NONCODING_MASKS: &[MaskDef] = &[
    (MaskType::Upstream, is_upstream),
    (MaskType::Downstream, is_downstream),
    (MaskType::Utr, is_utr),
    (MaskType::PromoterCage, is_promoter_cage),
    (MaskType::PromoterDhs, is_promoter_dhs),
    (MaskType::EnhancerCage, is_enhancer_cage),
    (MaskType::EnhancerDhs, is_enhancer_dhs),
    (MaskType::Ncrna, is_ncrna),
];

fn is_plof(v: &AnnotatedVariant) -> bool {
    v.annotation.consequence.is_plof() || v.annotation.region_type.contains_splicing()
}

fn is_missense(v: &AnnotatedVariant) -> bool {
    v.annotation.consequence.is_missense()
}

fn is_disruptive_missense(v: &AnnotatedVariant) -> bool {
    is_missense(v) && (v.annotation.cadd_phred > 20.0 || v.annotation.revel > 0.5)
}

fn is_plof_or_missense(v: &AnnotatedVariant) -> bool {
    is_plof(v) || is_missense(v)
}

fn is_synonymous(v: &AnnotatedVariant) -> bool {
    v.annotation.consequence.is_synonymous()
}

fn is_ptv(v: &AnnotatedVariant) -> bool {
    v.annotation.consequence.is_ptv()
}

fn is_ptv_ds(v: &AnnotatedVariant) -> bool {
    is_ptv(v) || (is_splice(v) && v.annotation.cadd_phred > 20.0)
}

fn is_splice(v: &AnnotatedVariant) -> bool {
    v.annotation.region_type.contains_splicing() || v.annotation.consequence.is_splice()
}

fn is_upstream(v: &AnnotatedVariant) -> bool {
    v.annotation.region_type.contains_upstream()
        || v.annotation.consequence == Consequence::UpstreamGeneVariant
}

fn is_downstream(v: &AnnotatedVariant) -> bool {
    v.annotation.region_type.contains_downstream()
        || v.annotation.consequence == Consequence::DownstreamGeneVariant
}

fn is_utr(v: &AnnotatedVariant) -> bool {
    v.annotation.region_type.contains_utr()
        || matches!(
            v.annotation.consequence,
            Consequence::FivePrimeUtrVariant | Consequence::ThreePrimeUtrVariant
        )
}

fn is_promoter_cage(v: &AnnotatedVariant) -> bool {
    v.annotation.regulatory.cage_promoter
}

fn is_promoter_dhs(v: &AnnotatedVariant) -> bool {
    v.annotation.regulatory.ccre_promoter
}

fn is_enhancer_cage(v: &AnnotatedVariant) -> bool {
    v.annotation.regulatory.cage_enhancer
}

fn is_enhancer_dhs(v: &AnnotatedVariant) -> bool {
    v.annotation.regulatory.ccre_enhancer
}

fn is_ncrna(v: &AnnotatedVariant) -> bool {
    v.annotation.region_type.contains_ncrna()
        || matches!(
            v.annotation.consequence,
            Consequence::NonCodingTranscriptExonVariant | Consequence::NonCodingTranscriptVariant
        )
}

pub fn build_masks_from_registry(
    variants: &[AnnotatedVariant],
    registry: &[MaskDef],
    min_variants: usize,
) -> Vec<(MaskType, Vec<MaskGroup>)> {
    let mut gene_variants: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();
    for (i, v) in variants.iter().enumerate() {
        if !v.gene_name.is_empty() {
            gene_variants
                .entry(v.gene_name.to_string())
                .or_default()
                .push(i);
        }
    }

    registry
        .iter()
        .map(|(mask_type, predicate)| {
            let groups: Vec<MaskGroup> = gene_variants
                .iter()
                .filter_map(|(gene, indices)| {
                    let matching: Vec<usize> = indices
                        .iter()
                        .copied()
                        .filter(|&i| predicate(&variants[i]))
                        .collect();
                    if matching.len() < min_variants {
                        return None;
                    }
                    let positions: Vec<u32> =
                        matching.iter().map(|&i| variants[i].position).collect();
                    Some(MaskGroup {
                        name: gene.clone(),
                        chromosome: variants[matching[0]].chromosome,
                        start: positions.iter().copied().min().unwrap_or(0),
                        end: positions.iter().copied().max().unwrap_or(0),
                        variant_indices: matching,
                    })
                })
                .collect();
            (mask_type.clone(), groups)
        })
        .collect()
}

pub fn build_coding_masks(
    variants: &[AnnotatedVariant],
    min_variants: usize,
) -> Vec<(MaskType, Vec<MaskGroup>)> {
    build_masks_from_registry(variants, CODING_MASKS, min_variants)
}

pub fn build_noncoding_masks(
    variants: &[AnnotatedVariant],
    min_variants: usize,
) -> Vec<(MaskType, Vec<MaskGroup>)> {
    build_masks_from_registry(variants, NONCODING_MASKS, min_variants)
}

// ---------------------------------------------------------------------------
// SCANG: data-adaptive variable-size scanning windows (formerly scang.rs)
// ---------------------------------------------------------------------------

/// SCANG parameters: data-adaptive variable-size scanning windows.
pub struct ScangParams {
    pub lmin: usize,
    pub lmax: usize,
    pub step: usize,
}

impl Default for ScangParams {
    fn default() -> Self {
        Self {
            lmin: 40,
            lmax: 300,
            step: 10,
        }
    }
}

/// Build SCANG window groups across multiple variant-count sizes.
///
/// Returns `(window_length_variants, groups)` for each tested size.
/// Windows slide by 1 variant position for maximum resolution.
pub fn build_scang_windows(
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    chromosome: Chromosome,
    params: &ScangParams,
) -> Vec<(u32, Vec<MaskGroup>)> {
    let n = chrom_indices.len();
    if n < params.lmin {
        return Vec::new();
    }

    let lmax = params.lmax.min(n);
    let mut result = Vec::new();
    let chrom_label = chromosome.label();

    let mut wsize = params.lmin;
    while wsize <= lmax {
        let mut groups = Vec::new();

        for start in 0..=(n - wsize) {
            let window_indices: Vec<usize> = chrom_indices[start..start + wsize].to_vec();
            let first_pos = variants[window_indices[0]].position;
            let last_pos = variants[*window_indices.last().unwrap()].position;

            groups.push(MaskGroup {
                name: format!("scang_L{}_{chrom_label}:{first_pos}-{last_pos}", wsize),
                chromosome,
                start: first_pos,
                end: last_pos,
                variant_indices: window_indices,
            });
        }

        result.push((wsize as u32, groups));
        wsize += params.step;
    }

    result
}

pub fn build_sliding_windows(
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    chromosome: Chromosome,
    window_size: u32,
    step_size: u32,
) -> Vec<MaskGroup> {
    if chrom_indices.is_empty() {
        return Vec::new();
    }

    let min_pos = chrom_indices
        .iter()
        .map(|&i| variants[i].position)
        .min()
        .unwrap_or(0);
    let max_pos = chrom_indices
        .iter()
        .map(|&i| variants[i].position)
        .max()
        .unwrap_or(0);

    let chrom_label = chromosome.label();
    let mut windows = Vec::new();
    let mut start = min_pos;

    while start <= max_pos {
        let end = start + window_size;
        let indices: Vec<usize> = chrom_indices
            .iter()
            .copied()
            .filter(|&i| variants[i].position >= start && variants[i].position < end)
            .collect();

        if indices.len() >= 2 {
            windows.push(MaskGroup {
                name: format!("{chrom_label}:{start}-{end}"),
                chromosome,
                start,
                end,
                variant_indices: indices,
            });
        }
        start += step_size;
    }
    windows
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn v(
    region_type: &str,
    consequence: &str,
    cadd: f64,
    revel: f64,
    cage_prom: bool,
    cage_enh: bool,
    ccre_prom: bool,
    ccre_enh: bool,
) -> AnnotatedVariant {
    use crate::types::{
        AnnotationWeights, Chromosome, FunctionalAnnotation, RegionType, RegulatoryFlags,
    };
    AnnotatedVariant {
        chromosome: Chromosome::Autosome(22),
        position: 1000,
        ref_allele: "A".into(),
        alt_allele: "T".into(),
        maf: 0.001,
        gene_name: "BRCA2".into(),
        annotation: FunctionalAnnotation {
            region_type: RegionType::from_str_lossy(region_type),
            consequence: Consequence::from_str_lossy(consequence),
            cadd_phred: cadd,
            revel,
            weights: AnnotationWeights([0.0; 11]),
            regulatory: RegulatoryFlags {
                cage_promoter: cage_prom,
                cage_enhancer: cage_enh,
                ccre_promoter: ccre_prom,
                ccre_enhancer: ccre_enh,
            },
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -- coding masks --

    #[test]
    fn stopgain_is_plof_and_ptv() {
        let var = v("exonic", "stopgain", 23.7, 0.0, false, false, false, false);
        assert!(is_plof(&var));
        assert!(is_ptv(&var));
        assert!(is_ptv_ds(&var));
        assert!(is_plof_or_missense(&var));
        assert!(!is_missense(&var));
        assert!(!is_synonymous(&var));
    }

    #[test]
    fn frameshift_is_plof_and_ptv() {
        let var = v(
            "exonic",
            "frameshift deletion",
            20.9,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_plof(&var));
        assert!(is_ptv(&var));
    }

    #[test]
    fn splice_region_type_is_plof_not_ptv() {
        let var = v("splicing", "", 25.1, 0.0, false, false, false, false);
        assert!(is_plof(&var));
        assert!(!is_ptv(&var));
        assert!(is_ptv_ds(&var));
        assert!(is_splice(&var));
    }

    #[test]
    fn splice_low_cadd_not_ptv_ds() {
        let var = v("splicing", "", 5.0, 0.0, false, false, false, false);
        assert!(is_plof(&var));
        assert!(!is_ptv(&var));
        assert!(!is_ptv_ds(&var));
    }

    #[test]
    fn missense_high_cadd_is_disruptive() {
        let var = v(
            "exonic",
            "nonsynonymous SNV",
            25.1,
            0.3,
            false,
            false,
            false,
            false,
        );
        assert!(is_missense(&var));
        assert!(is_disruptive_missense(&var));
        assert!(!is_plof(&var));
        assert!(!is_ptv(&var));
    }

    #[test]
    fn missense_low_scores_not_disruptive() {
        let var = v(
            "exonic",
            "nonsynonymous SNV",
            10.0,
            0.2,
            false,
            false,
            false,
            false,
        );
        assert!(is_missense(&var));
        assert!(!is_disruptive_missense(&var));
    }

    #[test]
    fn missense_high_revel_is_disruptive() {
        let var = v(
            "exonic",
            "nonsynonymous SNV",
            5.0,
            0.6,
            false,
            false,
            false,
            false,
        );
        assert!(is_disruptive_missense(&var));
    }

    #[test]
    fn synonymous_matches() {
        let var = v(
            "exonic",
            "synonymous SNV",
            6.4,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_synonymous(&var));
        assert!(!is_plof(&var));
        assert!(!is_missense(&var));
    }

    #[test]
    fn stoploss_is_plof_not_ptv() {
        let var = v("exonic", "stoploss", 15.0, 0.0, false, false, false, false);
        assert!(is_plof(&var));
        assert!(!is_ptv(&var));
    }

    // -- noncoding masks --

    #[test]
    fn upstream_matches() {
        let var = v("upstream", "", 5.6, 0.0, false, false, false, false);
        assert!(is_upstream(&var));
        assert!(!is_downstream(&var));
    }

    #[test]
    fn compound_upstream_downstream() {
        let var = v(
            "upstream;downstream",
            "",
            3.0,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_upstream(&var));
        assert!(is_downstream(&var));
    }

    #[test]
    fn utr3_matches() {
        let var = v("UTR3", "", 7.4, 0.0, false, false, false, false);
        assert!(is_utr(&var));
        assert!(!is_upstream(&var));
    }

    #[test]
    fn utr5_utr3_compound() {
        let var = v("UTR5;UTR3", "", 4.0, 0.0, false, false, false, false);
        assert!(is_utr(&var));
    }

    #[test]
    fn cage_promoter_flag() {
        let var = v("intergenic", "", 12.3, 0.0, true, false, false, false);
        assert!(is_promoter_cage(&var));
        assert!(!is_enhancer_cage(&var));
        assert!(!is_promoter_dhs(&var));
    }

    #[test]
    fn cage_enhancer_flag() {
        let var = v("intergenic", "", 8.1, 0.0, false, true, false, false);
        assert!(is_enhancer_cage(&var));
        assert!(!is_promoter_cage(&var));
    }

    #[test]
    fn ccre_pls_is_promoter_dhs() {
        let var = v("intergenic", "", 7.0, 0.0, false, false, true, false);
        assert!(is_promoter_dhs(&var));
        assert!(!is_enhancer_dhs(&var));
    }

    #[test]
    fn ccre_els_is_enhancer_dhs() {
        let var = v("intergenic", "", 8.7, 0.0, false, false, false, true);
        assert!(is_enhancer_dhs(&var));
        assert!(!is_promoter_dhs(&var));
    }

    #[test]
    fn ncrna_exonic() {
        let var = v("ncRNA_exonic", "", 9.9, 0.0, false, false, false, false);
        assert!(is_ncrna(&var));
        assert!(!is_plof(&var));
        assert!(!is_upstream(&var));
    }

    #[test]
    fn ncrna_compound_splicing() {
        let var = v(
            "ncRNA_exonic;splicing",
            "",
            15.0,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_ncrna(&var));
    }

    // -- negative cases: intronic should match nothing --

    #[test]
    fn intronic_matches_no_mask() {
        let var = v("intronic", "", 13.9, 0.0, false, false, false, false);
        assert!(!is_plof(&var));
        assert!(!is_missense(&var));
        assert!(!is_synonymous(&var));
        assert!(!is_ptv(&var));
        assert!(!is_ptv_ds(&var));
        assert!(!is_upstream(&var));
        assert!(!is_downstream(&var));
        assert!(!is_utr(&var));
        assert!(!is_promoter_cage(&var));
        assert!(!is_promoter_dhs(&var));
        assert!(!is_enhancer_cage(&var));
        assert!(!is_enhancer_dhs(&var));
        assert!(!is_ncrna(&var));
    }

    // -- mask group building --

    #[test]
    fn build_coding_masks_groups_by_gene() {
        let variants = vec![
            AnnotatedVariant {
                position: 100,
                gene_name: "TP53".into(),
                ..v("exonic", "stopgain", 30.0, 0.0, false, false, false, false)
            },
            AnnotatedVariant {
                position: 200,
                gene_name: "TP53".into(),
                ..v(
                    "exonic",
                    "frameshift deletion",
                    25.0,
                    0.0,
                    false,
                    false,
                    false,
                    false,
                )
            },
            AnnotatedVariant {
                position: 300,
                gene_name: "BRCA1".into(),
                ..v(
                    "exonic",
                    "nonsynonymous SNV",
                    20.0,
                    0.0,
                    false,
                    false,
                    false,
                    false,
                )
            },
        ];
        let masks = build_coding_masks(&variants, 2);
        let plof = masks.iter().find(|(mt, _)| *mt == MaskType::PLof);
        assert!(plof.is_some());
        let (_, groups) = plof.unwrap();
        assert_eq!(groups.len(), 1); // only TP53 has >=2 pLoF
        assert_eq!(groups[0].name, "TP53");
        assert_eq!(groups[0].variant_indices.len(), 2);
    }

    #[test]
    fn sliding_windows_non_overlapping() {
        let variants: Vec<AnnotatedVariant> = (0..5)
            .map(|i| AnnotatedVariant {
                position: i * 500,
                ..v("exonic", "stopgain", 10.0, 0.0, false, false, false, false)
            })
            .collect();
        let indices: Vec<usize> = (0..5).collect();
        let windows =
            build_sliding_windows(&variants, &indices, Chromosome::Autosome(22), 1000, 1000);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].variant_indices.len(), 2);
        assert_eq!(windows[1].variant_indices.len(), 2);
    }

    // -- VEP-style consequences (forward compat) --

    #[test]
    fn vep_missense_variant() {
        let var = v(
            "exonic",
            "missense_variant",
            25.0,
            0.7,
            false,
            false,
            false,
            false,
        );
        assert!(is_missense(&var));
        assert!(is_disruptive_missense(&var));
    }

    #[test]
    fn vep_upstream_gene_variant() {
        let var = v(
            "",
            "upstream_gene_variant",
            5.0,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_upstream(&var));
    }

    #[test]
    fn vep_splice_donor_high_cadd_is_ptv_ds() {
        let var = v(
            "",
            "splice_donor_variant",
            30.0,
            0.0,
            false,
            false,
            false,
            false,
        );
        assert!(is_plof(&var));
        assert!(is_splice(&var));
        assert!(is_ptv_ds(&var));
        assert!(!is_ptv(&var));
    }

    // -- SCANG tests --

    fn var_at(pos: u32) -> AnnotatedVariant {
        use crate::types::{
            AnnotationWeights, Chromosome, FunctionalAnnotation, RegionType, RegulatoryFlags,
        };
        AnnotatedVariant {
            chromosome: Chromosome::Autosome(1),
            position: pos,
            ref_allele: "A".into(),
            alt_allele: "T".into(),
            maf: 0.005,
            gene_name: "".into(),
            annotation: FunctionalAnnotation {
                region_type: RegionType::Unknown,
                consequence: Consequence::Unknown,
                cadd_phred: 10.0,
                revel: 0.0,
                weights: AnnotationWeights([0.5; 11]),
                regulatory: RegulatoryFlags::default(),
            },
        }
    }

    #[test]
    fn scang_windows_basic() {
        let variants: Vec<AnnotatedVariant> = (0..100).map(|i| var_at(i * 100)).collect();
        let indices: Vec<usize> = (0..100).collect();
        let params = ScangParams {
            lmin: 10,
            lmax: 30,
            step: 10,
        };
        let result = build_scang_windows(&variants, &indices, Chromosome::Autosome(1), &params);

        assert_eq!(result.len(), 3);
        assert_eq!(result[0].0, 10);
        assert_eq!(result[1].0, 20);
        assert_eq!(result[2].0, 30);

        assert_eq!(result[0].1.len(), 91);
        assert_eq!(result[1].1.len(), 81);
        assert_eq!(result[2].1.len(), 71);
    }

    #[test]
    fn scang_too_few_variants() {
        let variants: Vec<AnnotatedVariant> = (0..5).map(|i| var_at(i * 100)).collect();
        let indices: Vec<usize> = (0..5).collect();
        let params = ScangParams::default();
        let result = build_scang_windows(&variants, &indices, Chromosome::Autosome(1), &params);
        assert!(result.is_empty());
    }
}

//! Parity tests against R ground truth (scripts/r/generate_grm_ground_truth.R).

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use arrow::array::*;
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use faer::Mat;
use parquet::arrow::ArrowWriter;

use crate::column::{Col, STAAR_PHRED_CHANNELS};
use crate::error::CohortError;
use crate::output::{Output, Progress};
use crate::store::cohort::encoding::{SparseGHeader, SPARSE_G_HEADER_SIZE};
use crate::store::cohort::{ChromInfo, CohortId, CohortManifest};
use crate::store::config::StoreConfig;
use crate::store::Store;
use crate::types::Chromosome;

use super::pca;
use super::types::PcaScores;
use super::unrelated;

struct NullOutput;
impl Output for NullOutput {
    fn status(&self, _: &str) {}
    fn success(&self, _: &str) {}
    fn warn(&self, _: &str) {}
    fn error(&self, _: &CohortError) {}
    fn result_json(&self, _: &serde_json::Value) {}
    fn table(&self, _: &[&str], _: &[Vec<String>]) {}
    fn progress(&self, _: u64, _: &str) -> Progress { Progress::noop() }
}

#[derive(serde::Deserialize)]
struct Fixture {
    n_samples: usize,
    n_var_chr1: usize,
    n_var_chr2: usize,
    unrelated_mask: Vec<u8>,
    carriers_chr1: Vec<Vec<[u32; 2]>>,
    carriers_chr2: Vec<Vec<[u32; 2]>>,

    king_pairs_i: Vec<usize>,
    king_pairs_j: Vec<usize>,
    king_pairs_ibd: Vec<f64>,
    divergence: Vec<i32>,

    expected_mu_chr1: Vec<f64>,
    expected_mu_chr2: Vec<f64>,
    expected_inv_sd_chr1: Vec<f64>,
    expected_inv_sd_chr2: Vec<f64>,

    v_post_rows: usize,
    v_post_cols: usize,
    v_post: Vec<f64>,
    expected_post_chr1: Vec<f64>,
    expected_post_chr2: Vec<f64>,

    v_pre_chr1: Vec<f64>,
    v_pre_chr2: Vec<f64>,
    expected_pre_chr1: Vec<f64>,
    expected_pre_chr2: Vec<f64>,

    expected_pca_eigenvalues: Vec<f64>,
    expected_pca_scores_rows: usize,
    expected_pca_scores_cols: usize,
    expected_pca_scores: Vec<f64>,

    kinship_pair_i: Vec<usize>,
    kinship_pair_j: Vec<usize>,
    expected_kinship: Vec<f64>,

    expected_unrelated: Vec<usize>,
}

fn load_fixture() -> Fixture {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("src/staar/grm/testdata/ground_truth.json");
    let data = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&data).expect("parse ground_truth.json")
}

fn row_major_to_mat(flat: &[f64], rows: usize, cols: usize) -> Mat<f64> {
    Mat::from_fn(rows, cols, |r, c| flat[r * cols + c])
}

fn unrel_mask(fixture: &Fixture) -> Vec<bool> {
    fixture.unrelated_mask.iter().map(|&v| v != 0).collect()
}

fn build_test_cohort(fixture: &Fixture, dir: &Path) -> (Store, CohortId) {
    let store = Store::open(StoreConfig { root: dir.join("store") }).unwrap();
    let cid = CohortId::new("grm_test");
    let cohort_dir = dir.join("store/cohorts/grm_test");
    std::fs::create_dir_all(&cohort_dir).unwrap();

    write_samples(&cohort_dir, fixture.n_samples);
    write_manifest(&cohort_dir, fixture);

    for (name, n_var, carriers) in [
        ("1", fixture.n_var_chr1, &fixture.carriers_chr1),
        ("2", fixture.n_var_chr2, &fixture.carriers_chr2),
    ] {
        let chrom_dir = cohort_dir.join(format!("chromosome={name}"));
        std::fs::create_dir_all(&chrom_dir).unwrap();
        write_sparse_g(&chrom_dir, fixture.n_samples, carriers);
        write_stub_variants_parquet(&chrom_dir, n_var, name);
        write_stub_membership_parquet(&chrom_dir);
    }

    (store, cid)
}

fn write_samples(cohort_dir: &Path, n: usize) {
    let mut f = std::fs::File::create(cohort_dir.join("samples.txt")).unwrap();
    for i in 0..n {
        writeln!(f, "S{i:02}").unwrap();
    }
}

fn write_manifest(cohort_dir: &Path, fix: &Fixture) {
    let manifest = CohortManifest {
        version: 5,
        key: "test".into(),
        n_samples: fix.n_samples,
        n_variants: fix.n_var_chr1 + fix.n_var_chr2,
        chromosomes: vec![
            ChromInfo { name: "1".into(), n_variants: fix.n_var_chr1 },
            ChromInfo { name: "2".into(), n_variants: fix.n_var_chr2 },
        ],
        created_at: "0".into(),
        cohort_version: "test".into(),
    };
    let json = serde_json::to_string_pretty(&manifest).unwrap();
    std::fs::write(cohort_dir.join("manifest.json"), json).unwrap();
}

fn write_sparse_g(chrom_dir: &Path, n_samples: usize, carriers: &[Vec<[u32; 2]>]) {
    let n_variants = carriers.len();
    let mut body = Vec::new();
    let mut offsets = Vec::new();

    for cl in carriers {
        offsets.push(body.len() as u64);
        let n = cl.len() as u16;
        body.extend_from_slice(&n.to_le_bytes());
        for &[sample, dosage] in cl {
            body.extend_from_slice(&(sample as u16).to_le_bytes());
            body.push(dosage as u8);
        }
    }

    let total_carriers: u64 = carriers.iter().map(|c| c.len() as u64).sum();
    let offsets_start = (SPARSE_G_HEADER_SIZE + body.len()) as u64;

    let header = SparseGHeader::new(
        n_samples as u32,
        n_variants as u32,
        total_carriers,
        offsets_start,
    );

    let path = chrom_dir.join("sparse_g.bin");
    let mut f = std::fs::File::create(&path).unwrap();
    let mut buf = Vec::new();
    header.write_to(&mut buf).unwrap();
    f.write_all(&buf).unwrap();
    f.write_all(&body).unwrap();
    for off in &offsets {
        f.write_all(&off.to_le_bytes()).unwrap();
    }
}

fn write_stub_variants_parquet(chrom_dir: &Path, n_variants: usize, chrom: &str) {
    let n = n_variants;
    let positions: Vec<i32> = (1..=n as i32).collect();
    let vids: Vec<String> = (0..n).map(|i| format!("{chrom}-{}-A-T", i + 1)).collect();
    let vid_refs: Vec<&str> = vids.iter().map(|s| s.as_str()).collect();

    let mut fields = vec![
        Field::new(Col::Position.as_str(), DataType::Int32, false),
        Field::new(Col::EndPosition.as_str(), DataType::Int32, false),
        Field::new(Col::RefAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::AltAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::Vid.as_str(), DataType::Utf8, false),
        Field::new(Col::Maf.as_str(), DataType::Float64, false),
        Field::new(Col::RegionType.as_str(), DataType::Utf8, false),
        Field::new(Col::Consequence.as_str(), DataType::Utf8, false),
        Field::new(Col::Revel.as_str(), DataType::Float64, false),
        Field::new(Col::MetaSvmPred.as_str(), DataType::Utf8, false),
        Field::new(Col::GeneHancer.as_str(), DataType::Utf8, false),
        Field::new(Col::IsCagePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCageEnhancer.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcrePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcreEnhancer.as_str(), DataType::Boolean, false),
    ];
    for col in &STAAR_PHRED_CHANNELS {
        fields.push(Field::new(col.as_str(), DataType::Float64, true));
    }
    let schema = Arc::new(Schema::new(fields));

    let str_n = vec!["A"; n];
    let bools = vec![false; n];
    let phreds = vec![10.0f64; n];
    let mafs = vec![0.1f64; n];

    let mut columns: Vec<Arc<dyn Array>> = vec![
        Arc::new(Int32Array::from(positions.clone())),
        Arc::new(Int32Array::from(positions)),
        Arc::new(StringArray::from(str_n.clone())),
        Arc::new(StringArray::from(vec!["T"; n])),
        Arc::new(StringArray::from(vid_refs)),
        Arc::new(Float64Array::from(mafs)),
        Arc::new(StringArray::from(vec!["exonic"; n])),
        Arc::new(StringArray::from(vec!["missense"; n])),
        Arc::new(Float64Array::from(vec![0.5f64; n])),
        Arc::new(StringArray::from(str_n.clone())),
        Arc::new(StringArray::from(vec!["."; n])),
        Arc::new(BooleanArray::from(bools.clone())),
        Arc::new(BooleanArray::from(bools.clone())),
        Arc::new(BooleanArray::from(bools.clone())),
        Arc::new(BooleanArray::from(bools)),
    ];
    for _ in &STAAR_PHRED_CHANNELS {
        columns.push(Arc::new(Float64Array::from(phreds.clone())));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns).unwrap();
    let file = std::fs::File::create(chrom_dir.join("variants.parquet")).unwrap();
    let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
    writer.write(&batch).unwrap();
    writer.close().unwrap();
}

fn write_stub_membership_parquet(chrom_dir: &Path) {
    let schema = Arc::new(Schema::new(vec![
        Field::new("variant_vcf", DataType::UInt32, false),
        Field::new("gene", DataType::Utf8, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(UInt32Array::from(Vec::<u32>::new())),
            Arc::new(StringArray::from(Vec::<&str>::new())),
        ],
    )
    .unwrap();
    let file = std::fs::File::create(chrom_dir.join("membership.parquet")).unwrap();
    let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
    writer.write(&batch).unwrap();
    writer.close().unwrap();
}

// -- Tests --

#[test]
fn allele_freq_matches_r() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let mask = unrel_mask(&fix);

    for (chrom_n, expected_mu, expected_inv_sd) in [
        (1, &fix.expected_mu_chr1, &fix.expected_inv_sd_chr1),
        (2, &fix.expected_mu_chr2, &fix.expected_inv_sd_chr2),
    ] {
        let view = handle.chromosome(&Chromosome::Autosome(chrom_n)).unwrap();
        let stats = pca::allele_freq_chrom(&view, &mask).unwrap();
        for (i, (&got, &want)) in stats.mu.iter().zip(expected_mu.iter()).enumerate() {
            assert!(
                (got - want).abs() < 1e-12,
                "chr{chrom_n} mu[{i}]: got {got}, want {want}"
            );
        }
        for (i, (&got, &want)) in stats.inv_sd.iter().zip(expected_inv_sd.iter()).enumerate() {
            assert!(
                (got - want).abs() < 1e-8,
                "chr{chrom_n} inv_sd[{i}]: got {got}, want {want}"
            );
        }
    }
}

#[test]
fn postmultiply_matches_r() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let mask = unrel_mask(&fix);
    let v = row_major_to_mat(&fix.v_post, fix.v_post_rows, fix.v_post_cols);

    for (chrom_n, n_var, expected_flat) in [
        (1, fix.n_var_chr1, &fix.expected_post_chr1),
        (2, fix.n_var_chr2, &fix.expected_post_chr2),
    ] {
        let expected = row_major_to_mat(expected_flat, n_var, fix.v_post_cols);
        let view = handle.chromosome(&Chromosome::Autosome(chrom_n)).unwrap();
        let stats = pca::allele_freq_chrom(&view, &mask).unwrap();
        let got = pca::postmultiply_chrom(&view, &v, &mask, &stats).unwrap();
        for r in 0..n_var {
            for c in 0..fix.v_post_cols {
                let g = got[(r, c)];
                let w = expected[(r, c)];
                assert!(
                    (g - w).abs() < 1e-8,
                    "post chr{chrom_n}[{r},{c}]: got {g}, want {w}"
                );
            }
        }
    }
}

#[test]
fn premultiply_matches_r() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let mask = unrel_mask(&fix);

    for (chrom_n, n_var, v_flat, expected_flat) in [
        (1, fix.n_var_chr1, &fix.v_pre_chr1, &fix.expected_pre_chr1),
        (2, fix.n_var_chr2, &fix.v_pre_chr2, &fix.expected_pre_chr2),
    ] {
        let v = row_major_to_mat(v_flat, n_var, 2);
        let expected = row_major_to_mat(expected_flat, fix.n_samples, 2);
        let view = handle.chromosome(&Chromosome::Autosome(chrom_n)).unwrap();
        let stats = pca::allele_freq_chrom(&view, &mask).unwrap();
        let mut got = Mat::<f64>::zeros(fix.n_samples, 2);
        pca::premultiply_chrom(&view, &v, &mask, &stats, &mut got).unwrap();
        for r in 0..fix.n_samples {
            for c in 0..2 {
                let g = got[(r, c)];
                let w = expected[(r, c)];
                assert!(
                    (g - w).abs() < 1e-8,
                    "pre chr{chrom_n}[{r},{c}]: got {g}, want {w}"
                );
            }
        }
    }
}

#[test]
fn postpremultiply_round_trip() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let mask = unrel_mask(&fix);
    let v = row_major_to_mat(&fix.v_post, fix.v_post_rows, fix.v_post_cols);

    let view = handle.chromosome(&Chromosome::Autosome(1)).unwrap();
    let stats = pca::allele_freq_chrom(&view, &mask).unwrap();
    let gv = pca::postmultiply_chrom(&view, &v, &mask, &stats).unwrap();
    let mut gtgv = Mat::<f64>::zeros(fix.n_samples, fix.v_post_cols);
    pca::premultiply_chrom(&view, &gv, &mask, &stats, &mut gtgv).unwrap();

    let view2 = handle.chromosome(&Chromosome::Autosome(1)).unwrap();
    let stats2 = pca::allele_freq_chrom(&view2, &mask).unwrap();
    let gv2 = pca::postmultiply_chrom(&view2, &v, &mask, &stats2).unwrap();
    let mut gtgv2 = Mat::<f64>::zeros(fix.n_samples, fix.v_post_cols);
    pca::premultiply_chrom(&view2, &gv2, &mask, &stats2, &mut gtgv2).unwrap();

    for r in 0..fix.n_samples {
        for c in 0..fix.v_post_cols {
            let a = gtgv[(r, c)];
            let b = gtgv2[(r, c)];
            assert!((a - b).abs() < 1e-10, "round-trip [{r},{c}]: {a} vs {b}");
        }
    }
}

#[test]
fn pca_eigenvalues_within_5_percent() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let result = handle.load().unwrap();
    let mask = unrel_mask(&fix);
    let all_mask = vec![true; fix.n_samples];
    let n_pcs = fix.expected_pca_scores_cols;

    let pca_result = pca::randomized_pca(
        &handle,
        &result.manifest,
        &mask,
        &all_mask,
        n_pcs,
        10,
    )
    .unwrap();

    for (i, (&got, &want)) in pca_result
        .eigenvalues
        .iter()
        .zip(fix.expected_pca_eigenvalues.iter())
        .enumerate()
    {
        let rel_err = (got - want).abs() / want;
        assert!(
            rel_err < 0.05,
            "PCA eigenvalue[{i}]: got {got}, want {want} (rel err {rel_err:.4})"
        );
    }
}

#[test]
fn pca_score_correlations_above_099() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let result = handle.load().unwrap();
    let mask = unrel_mask(&fix);
    let all_mask = vec![true; fix.n_samples];
    let n_pcs = fix.expected_pca_scores_cols;

    let pca_result = pca::randomized_pca(
        &handle,
        &result.manifest,
        &mask,
        &all_mask,
        n_pcs,
        10,
    )
    .unwrap();

    let expected =
        row_major_to_mat(&fix.expected_pca_scores, fix.expected_pca_scores_rows, n_pcs);

    for pc in 0..n_pcs {
        let corr = column_correlation(&pca_result.scores, &expected, pc, fix.n_samples);
        assert!(
            corr.abs() > 0.99,
            "PCA score PC{}: correlation {corr:.4} (sign flip ok)",
            pc + 1
        );
    }
}

fn column_correlation(a: &Mat<f64>, b: &Mat<f64>, col: usize, n: usize) -> f64 {
    let mut sum_a = 0.0;
    let mut sum_b = 0.0;
    for i in 0..n {
        sum_a += a[(i, col)];
        sum_b += b[(i, col)];
    }
    let mean_a = sum_a / n as f64;
    let mean_b = sum_b / n as f64;
    let mut cov = 0.0;
    let mut var_a = 0.0;
    let mut var_b = 0.0;
    for i in 0..n {
        let da = a[(i, col)] - mean_a;
        let db = b[(i, col)] - mean_b;
        cov += da * db;
        var_a += da * da;
        var_b += db * db;
    }
    if var_a < 1e-20 || var_b < 1e-20 {
        return 0.0;
    }
    cov / (var_a * var_b).sqrt()
}

#[test]
fn kinship_matches_r() {
    let fix = load_fixture();
    let dir = tempfile::tempdir().unwrap();
    let (store, cid) = build_test_cohort(&fix, dir.path());
    let handle = store.cohort(&cid);
    let result = handle.load().unwrap();
    let mask = unrel_mask(&fix);
    let n_pcs = fix.expected_pca_scores_cols;

    let pca_scores_mat =
        row_major_to_mat(&fix.expected_pca_scores, fix.expected_pca_scores_rows, n_pcs);
    let pca_input = PcaScores {
        scores: pca_scores_mat,
        eigenvalues: fix.expected_pca_eigenvalues.clone(),
    };

    let pairs: Vec<(usize, usize, f64)> = fix
        .king_pairs_i
        .iter()
        .zip(fix.king_pairs_j.iter())
        .zip(fix.king_pairs_ibd.iter())
        .map(|((&i, &j), &ibd)| (i, j, ibd))
        .collect();

    let out = NullOutput;
    let grm = super::estimate::estimate_grm(
        &handle,
        &result.manifest,
        &pca_input,
        &mask,
        &pairs,
        fix.n_samples,
        n_pcs,
        5000,
        4,
        &out,
    )
    .unwrap();

    let expected_kin: std::collections::HashMap<(usize, usize), f64> = fix
        .kinship_pair_i
        .iter()
        .zip(fix.kinship_pair_j.iter())
        .zip(fix.expected_kinship.iter())
        .map(|((&i, &j), &k)| ((i, j), k))
        .collect();

    for &(i, j, kin) in &grm.triplets {
        if i > j {
            continue;
        }
        if let Some(&want) = expected_kin.get(&(i, j)) {
            assert!(
                (kin - want).abs() < 1e-4,
                "kinship({i},{j}): got {kin:.6}, want {want:.6}"
            );
        }
    }
}

#[test]
fn select_unrelated_matches_r() {
    let fix = load_fixture();
    let pairs: Vec<(usize, usize, f64)> = fix
        .king_pairs_i
        .iter()
        .zip(fix.king_pairs_j.iter())
        .zip(fix.king_pairs_ibd.iter())
        .map(|((&i, &j), &ibd)| (i, j, ibd))
        .collect();

    let result = unrelated::select_unrelated(&pairs, fix.n_samples, &fix.divergence);
    let mut got = result.sample_indices.clone();
    got.sort();
    let mut want = fix.expected_unrelated.clone();
    want.sort();
    assert_eq!(got, want, "unrelated set mismatch");
}

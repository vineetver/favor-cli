//! Bit-identical regression net for the STAAR sumstats kernel.
//!
//! Loads the same `staar_continuous` fixture `ground_truth_test.rs` uses,
//! runs `run_staar_from_sumstats`, and compares every output f64 as raw
//! bits against a committed golden. No tolerance — any drift fails the
//! build. Refactors that claim to be mathematically equivalent must not
//! change the bits.
//!
//! Regenerate the golden (after an intentional numerical change) with:
//!     cargo test -- regenerate_invariance_golden --ignored

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use faer::Mat;
    use serde::{Deserialize, Serialize};

    use crate::staar::score;

    #[derive(Deserialize)]
    struct GroundTruthRoot {
        staar_continuous: StaarCase,
    }

    #[derive(Deserialize)]
    struct StaarCase {
        n_samples: usize,
        #[serde(rename = "U")]
        u: Vec<f64>,
        #[serde(rename = "K")]
        k: Vec<Vec<f64>>,
        mafs: Vec<f64>,
        annotation_rank: Vec<Vec<f64>>,
        sigma2: f64,
    }

    /// Golden is structurally isomorphic to `score::StaarResult` but every
    /// f64 is serialized as `to_bits()` in 16-char lowercase hex. That way
    /// equality is bit-exact and a diff is human-readable.
    #[derive(Serialize, Deserialize, PartialEq, Eq, Debug)]
    struct Golden {
        burden_1_25: String,
        burden_1_1: String,
        skat_1_25: String,
        skat_1_1: String,
        acat_v_1_25: String,
        acat_v_1_1: String,
        per_annotation: Vec<[String; 6]>,
        staar_b_1_25: String,
        staar_b_1_1: String,
        staar_s_1_25: String,
        staar_s_1_1: String,
        staar_a_1_25: String,
        staar_a_1_1: String,
        acat_o: String,
        staar_o: String,
    }

    fn bits(x: f64) -> String {
        format!("{:016x}", x.to_bits())
    }

    fn golden_of(r: &score::StaarResult) -> Golden {
        Golden {
            burden_1_25: bits(r.burden_1_25),
            burden_1_1: bits(r.burden_1_1),
            skat_1_25: bits(r.skat_1_25),
            skat_1_1: bits(r.skat_1_1),
            acat_v_1_25: bits(r.acat_v_1_25),
            acat_v_1_1: bits(r.acat_v_1_1),
            per_annotation: r
                .per_annotation
                .iter()
                .map(|row| [
                    bits(row[0]), bits(row[1]), bits(row[2]),
                    bits(row[3]), bits(row[4]), bits(row[5]),
                ])
                .collect(),
            staar_b_1_25: bits(r.staar_b_1_25),
            staar_b_1_1: bits(r.staar_b_1_1),
            staar_s_1_25: bits(r.staar_s_1_25),
            staar_s_1_1: bits(r.staar_s_1_1),
            staar_a_1_25: bits(r.staar_a_1_25),
            staar_a_1_1: bits(r.staar_a_1_1),
            acat_o: bits(r.acat_o),
            staar_o: bits(r.staar_o),
        }
    }

    fn testdata_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("src")
            .join("staar")
            .join("testdata")
    }

    fn load_inputs() -> StaarCase {
        let path = testdata_dir().join("ground_truth.json");
        let data = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
        // R serializes NULL as `{}` in some sections; normalize so serde parses.
        let data = data.replace(": {}", ": null");
        let root: GroundTruthRoot = serde_json::from_str(&data)
            .expect("failed to parse ground_truth.json");
        root.staar_continuous
    }

    fn run_kernel(case: &StaarCase) -> score::StaarResult {
        // The kernel expects U / sigma2 and K / sigma2 (sum-stat scaling);
        // ground_truth.json stores raw U = G'r and K = G'PG.
        let inv_s2 = 1.0 / case.sigma2;
        let u_scaled: Vec<f64> = case.u.iter().map(|&v| v * inv_s2).collect();
        let k_scaled: Vec<Vec<f64>> = case
            .k
            .iter()
            .map(|row| row.iter().map(|&v| v * inv_s2).collect())
            .collect();

        let u = Mat::from_fn(u_scaled.len(), 1, |i, _| u_scaled[i]);
        let k = Mat::from_fn(k_scaled.len(), k_scaled[0].len(), |i, j| k_scaled[i][j]);

        score::run_staar_from_sumstats(&u, &k, &case.annotation_rank, &case.mafs, case.n_samples)
    }

    #[test]
    fn invariance_bits_match_golden() {
        let case = load_inputs();
        let actual = golden_of(&run_kernel(&case));

        let golden_path = testdata_dir().join("invariance_golden.json");
        let text = std::fs::read_to_string(&golden_path).unwrap_or_else(|e| {
            panic!(
                "cannot read {}: {e}\n\
                 regenerate with: cargo test -- regenerate_invariance_golden --ignored",
                golden_path.display()
            )
        });
        let expected: Golden = serde_json::from_str(&text)
            .expect("failed to parse invariance_golden.json");

        assert_eq!(
            actual, expected,
            "STAAR output drift — run `cargo test -- regenerate_invariance_golden --ignored` \
             if the drift is intentional; otherwise the change is a regression."
        );
    }

    #[test]
    #[ignore = "regenerates the committed invariance_golden.json"]
    fn regenerate_invariance_golden() {
        let case = load_inputs();
        let golden = golden_of(&run_kernel(&case));
        let json = serde_json::to_string_pretty(&golden).unwrap();

        let path = testdata_dir().join("invariance_golden.json");
        std::fs::write(&path, json).unwrap_or_else(|e| {
            panic!("cannot write {}: {e}", path.display())
        });
        eprintln!("wrote {}", path.display());
    }
}

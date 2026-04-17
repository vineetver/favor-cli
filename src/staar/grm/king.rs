//! KING .seg file parser and connected-component discovery.
//!
//! Mirrors FastSparseGRM R/getUnrels.R:removeHigherDegree + R/calcGRM.R:29-68.
//! Parses the IBD segment output from KING, filters pairs by relatedness
//! degree, and builds connected components via union-find. Each component
//! is a cluster of related individuals whose pairwise kinships will be
//! estimated in the GRM step.

use std::collections::HashMap;
use std::path::Path;

use crate::error::CohortError;

use super::types::{KingSegEntry, RelatedComponent, RelatednessType};

/// Parse a KING .seg file and filter pairs up to the given degree.
///
/// Handles both KING column layouts:
///   FID1 ID1 FID2 ID2 ... PropIBD ... InfType
///   FID  ID1 FID  ID2 ... PropIBD ... InfType
///
/// Sample IDs are formatted as `FID_IID` to match upstream convention.
pub fn parse_king_seg(path: &Path, max_degree: u8) -> Result<Vec<KingSegEntry>, CohortError> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;

    let mut lines = content.lines();
    let header = lines.next().ok_or_else(|| {
        CohortError::Input(format!("KING .seg file '{}' is empty", path.display()))
    })?;

    let cols: Vec<&str> = header.split_whitespace().collect();
    let col_idx = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));

    let has_fid1 = col_idx("FID1").is_some();
    let (fid1_col, id1_col, fid2_col, id2_col) = if has_fid1 {
        (
            col_idx("FID1").unwrap(),
            col_idx("ID1").unwrap(),
            col_idx("FID2").unwrap(),
            col_idx("ID2").unwrap(),
        )
    } else {
        let fid = col_idx("FID").ok_or_else(|| {
            CohortError::Input(format!(
                "KING .seg '{}': missing FID1/FID column in header: {header}",
                path.display()
            ))
        })?;
        let id1 = col_idx("ID1").ok_or_else(|| {
            CohortError::Input(format!(
                "KING .seg '{}': missing ID1 column",
                path.display()
            ))
        })?;
        let id2 = col_idx("ID2").ok_or_else(|| {
            CohortError::Input(format!(
                "KING .seg '{}': missing ID2 column",
                path.display()
            ))
        })?;
        (fid, id1, fid, id2)
    };

    let ibd_col = col_idx("PropIBD")
        .or_else(|| col_idx("Kinship"))
        .ok_or_else(|| {
            CohortError::Input(format!(
                "KING .seg '{}': missing PropIBD/Kinship column",
                path.display()
            ))
        })?;
    let inf_col = col_idx("InfType").ok_or_else(|| {
        CohortError::Input(format!(
            "KING .seg '{}': missing InfType column",
            path.display()
        ))
    })?;
    let n_cols = cols.len();

    let mut entries = Vec::new();
    for (lineno, line) in lines.enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < n_cols {
            continue;
        }
        let inf_type = RelatednessType::from_king_label(parts[inf_col]);
        if inf_type == RelatednessType::Unrelated || inf_type.degree() > max_degree {
            continue;
        }
        let prop_ibd: f64 = parts[ibd_col].parse().map_err(|e| {
            CohortError::Input(format!(
                "KING .seg {}:{}: bad PropIBD '{}': {e}",
                path.display(),
                lineno + 2,
                parts[ibd_col]
            ))
        })?;
        let id1 = format!("{}_{}", parts[fid1_col], parts[id1_col]);
        let id2 = format!("{}_{}", parts[fid2_col], parts[id2_col]);
        entries.push(KingSegEntry {
            id1,
            id2,
            prop_ibd,
            inf_type,
        });
    }

    Ok(entries)
}

/// Map sample identifiers (FID_IID) to cohort sample indices. Pairs
/// referencing samples not in the cohort are silently dropped.
#[allow(clippy::type_complexity)]
pub fn map_to_cohort_indices(
    entries: &[KingSegEntry],
    sample_ids: &[String],
) -> (Vec<(usize, usize, f64)>, HashMap<String, usize>) {
    let id_to_idx: HashMap<String, usize> = sample_ids
        .iter()
        .enumerate()
        .map(|(i, s)| (s.clone(), i))
        .collect();

    let mut pairs = Vec::with_capacity(entries.len());
    for e in entries {
        if let (Some(&i), Some(&j)) = (id_to_idx.get(&e.id1), id_to_idx.get(&e.id2)) {
            let (lo, hi) = if i < j { (i, j) } else { (j, i) };
            pairs.push((lo, hi, e.prop_ibd));
        }
    }
    pairs.sort_unstable_by_key(|a| (a.0, a.1));
    pairs.dedup_by_key(|p| (p.0, p.1));
    (pairs, id_to_idx)
}

/// Build connected components from related pairs using union-find.
pub fn build_components(
    pairs: &[(usize, usize, f64)],
    n_samples: usize,
) -> Vec<RelatedComponent> {
    let mut parent: Vec<usize> = (0..n_samples).collect();
    let mut rank = vec![0u8; n_samples];

    let find = |parent: &mut [usize], mut x: usize| -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    };

    for &(i, j, _) in pairs {
        let ri = find(&mut parent, i);
        let rj = find(&mut parent, j);
        if ri != rj {
            if rank[ri] < rank[rj] {
                parent[ri] = rj;
            } else if rank[ri] > rank[rj] {
                parent[rj] = ri;
            } else {
                parent[rj] = ri;
                rank[ri] += 1;
            }
        }
    }

    let mut comp_members: HashMap<usize, Vec<usize>> = HashMap::new();
    let involved: std::collections::HashSet<usize> =
        pairs.iter().flat_map(|&(i, j, _)| [i, j]).collect();
    for &s in &involved {
        let root = find(&mut parent, s);
        comp_members.entry(root).or_default().push(s);
    }

    let mut components: Vec<RelatedComponent> = comp_members
        .into_values()
        .map(|mut members| {
            members.sort_unstable();
            let member_set: HashMap<usize, usize> = members
                .iter()
                .enumerate()
                .map(|(local, &global)| (global, local))
                .collect();
            let comp_pairs: Vec<(usize, usize)> = pairs
                .iter()
                .filter_map(|&(i, j, _)| {
                    let li = member_set.get(&i)?;
                    let lj = member_set.get(&j)?;
                    Some((*li, *lj))
                })
                .collect();
            RelatedComponent {
                members,
                pairs: comp_pairs,
            }
        })
        .collect();
    components.sort_by_key(|c| std::cmp::Reverse(c.members.len()));
    components
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg_content(rows: &[(&str, &str, &str, &str, &str, &str)]) -> String {
        let mut s = String::from("FID1\tID1\tFID2\tID2\tPropIBD\tInfType\n");
        for (f1, i1, f2, i2, ibd, inf) in rows {
            s.push_str(&format!("{f1}\t{i1}\t{f2}\t{i2}\t{ibd}\t{inf}\n"));
        }
        s
    }

    #[test]
    fn parse_filters_by_degree() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.seg");
        std::fs::write(
            &path,
            seg_content(&[
                ("fam1", "s1", "fam1", "s2", "0.5", "PO"),
                ("fam2", "s3", "fam2", "s4", "0.25", "2nd"),
                ("fam3", "s5", "fam3", "s6", "0.01", "4th"),
                ("fam4", "s7", "fam4", "s8", "0.001", "UN"),
            ]),
        )
        .unwrap();

        let all = parse_king_seg(&path, 4).unwrap();
        assert_eq!(all.len(), 3);

        let degree2 = parse_king_seg(&path, 2).unwrap();
        assert_eq!(degree2.len(), 2);

        let degree1 = parse_king_seg(&path, 1).unwrap();
        assert_eq!(degree1.len(), 1);
        assert_eq!(degree1[0].inf_type, RelatednessType::First);
    }

    #[test]
    fn build_components_finds_clusters() {
        let pairs = vec![(0, 1, 0.5), (1, 2, 0.25), (3, 4, 0.5)];
        let components = build_components(&pairs, 10);
        assert_eq!(components.len(), 2);
        let sizes: Vec<usize> = components.iter().map(|c| c.members.len()).collect();
        assert!(sizes.contains(&3));
        assert!(sizes.contains(&2));
    }

    #[test]
    fn singletons_excluded_from_components() {
        let pairs = vec![(0, 1, 0.5)];
        let components = build_components(&pairs, 100);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].members.len(), 2);
    }

    #[test]
    fn map_to_cohort_deduplicates() {
        let entries = vec![
            KingSegEntry {
                id1: "f_s1".into(),
                id2: "f_s2".into(),
                prop_ibd: 0.5,
                inf_type: RelatednessType::First,
            },
            KingSegEntry {
                id1: "f_s2".into(),
                id2: "f_s1".into(),
                prop_ibd: 0.5,
                inf_type: RelatednessType::First,
            },
        ];
        let samples = vec!["f_s1".into(), "f_s2".into(), "f_s3".into()];
        let (pairs, _) = map_to_cohort_indices(&entries, &samples);
        assert_eq!(pairs.len(), 1);
    }
}

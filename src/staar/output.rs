//! Result writing (parquet) and HTML summary report generation.
//!
//! Every artifact here is written via `write_parquet_atomic` or
//! `store::write_atomic` so a crash mid-write leaves either nothing or
//! the prior version, never a half-finished parquet on disk.

use std::collections::HashMap;
use std::fs::{self, File};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float64Builder, Int32Builder, StringBuilder, UInt32Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use serde_json::json;

use crate::column::{Col, STAAR_WEIGHTS};
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::store::{fsync_parent, tmp_path, write_atomic};
use crate::types::AnnotatedVariant;

use super::model::NullModel;
use super::{GeneResult, MaskType, TraitType};

/// Run `build` against `out_path.tmp`, then rename onto `out_path` and
/// fsync the parent. Closure shape so each caller keeps its arrow
/// builders local.
fn write_parquet_atomic<F>(out_path: &Path, build: F) -> Result<(), CohortError>
where
    F: FnOnce(File) -> Result<(), CohortError>,
{
    let tmp = tmp_path(out_path);
    let file = File::create(&tmp)
        .map_err(|e| CohortError::Resource(format!("Create {}: {e}", tmp.display())))?;
    build(file)?;
    fs::rename(&tmp, out_path).map_err(|e| {
        CohortError::Resource(format!(
            "rename {} -> {}: {e}",
            tmp.display(),
            out_path.display()
        ))
    })?;
    fsync_parent(out_path);
    Ok(())
}

pub fn write_individual_results(
    pvals: &[(usize, f64)],
    variants: &[AnnotatedVariant],
    output_dir: &Path,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let out_path = output_dir.join("individual.parquet");
    let n = pvals.len();

    let mut sorted: Vec<(usize, f64)> = pvals.to_vec();
    sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut b_chrom = StringBuilder::with_capacity(n, n * 2);
    let mut b_pos = Int32Builder::with_capacity(n);
    let mut b_ref = StringBuilder::with_capacity(n, n * 4);
    let mut b_alt = StringBuilder::with_capacity(n, n * 4);
    let mut b_maf = Float64Builder::with_capacity(n);
    let mut b_gene = StringBuilder::with_capacity(n, n * 8);
    let mut b_region = StringBuilder::with_capacity(n, n * 8);
    let mut b_consequence = StringBuilder::with_capacity(n, n * 8);
    let mut b_cadd = Float64Builder::with_capacity(n);
    let mut b_pvalue = Float64Builder::with_capacity(n);

    for &(idx, pval) in &sorted {
        let v = &variants[idx];
        b_chrom.append_value(v.chromosome.label());
        b_pos.append_value(v.position as i32);
        b_ref.append_value(&v.ref_allele);
        b_alt.append_value(&v.alt_allele);
        b_maf.append_value(v.maf);
        b_gene.append_value(&v.gene_name);
        b_region.append_value(v.annotation.region_type.as_str());
        b_consequence.append_value(v.annotation.consequence.as_str());
        b_cadd.append_value(v.annotation.cadd_phred);
        b_pvalue.append_value(pval);
    }

    let schema = Arc::new(Schema::new(vec![
        Field::new(Col::Chromosome.as_str(), DataType::Utf8, false),
        Field::new(Col::Position.as_str(), DataType::Int32, false),
        Field::new(Col::RefAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::AltAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::Maf.as_str(), DataType::Float64, false),
        Field::new(Col::GeneName.as_str(), DataType::Utf8, false),
        Field::new(Col::RegionType.as_str(), DataType::Utf8, false),
        Field::new(Col::Consequence.as_str(), DataType::Utf8, false),
        Field::new(Col::CaddPhred.as_str(), DataType::Float64, false),
        Field::new("pvalue", DataType::Float64, false),
    ]));
    let columns: Vec<ArrayRef> = vec![
        Arc::new(b_chrom.finish()),
        Arc::new(b_pos.finish()),
        Arc::new(b_ref.finish()),
        Arc::new(b_alt.finish()),
        Arc::new(b_maf.finish()),
        Arc::new(b_gene.finish()),
        Arc::new(b_region.finish()),
        Arc::new(b_consequence.finish()),
        Arc::new(b_cadd.finish()),
        Arc::new(b_pvalue.finish()),
    ];
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))?;

    write_parquet_atomic(&out_path, |file| {
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .build();
        let mut writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| CohortError::Resource(format!("Parquet writer: {e}")))?;
        writer
            .write(&batch)
            .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
        writer
            .close()
            .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
        Ok(())
    })?;

    let n_sig = pvals.iter().filter(|(_, p)| *p < 5e-8).count();
    out.success(&format!(
        "  individual -> {} variants, {} genome-wide significant",
        pvals.len(),
        n_sig
    ));
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn write_results(
    all_mask_results: &[(MaskType, Vec<GeneResult>)],
    trait_names: &[String],
    maf_cutoff: f64,
    output_dir: &Path,
    null_model: &NullModel,
    trait_type: TraitType,
    n: usize,
    n_rare: i64,
    out: &dyn Output,
) -> Result<(), CohortError> {
    out.status("Writing results...");
    let mut significant_genes: Vec<serde_json::Value> = Vec::new();

    let channels: Vec<&str> = STAAR_WEIGHTS
        .iter()
        .map(|c| c.weight_display_name().expect("STAAR_WEIGHTS entries have display names"))
        .collect();
    let n_channels = channels.len();

    for (mask_type, results) in all_mask_results {
        if results.is_empty() {
            continue;
        }
        let out_path = output_dir.join(format!("{}.parquet", mask_type.file_stem()));

        let nan_pvals = results.iter().filter(|r| r.staar.staar_o.is_nan()).count();
        if nan_pvals > 0 {
            out.warn(&format!(
                "  {}: {nan_pvals} gene-mask combinations produced NaN p-values",
                mask_type.file_stem()
            ));
        }

        let mut sorted_results: Vec<&GeneResult> = results.iter().collect();
        sorted_results.sort_by(|a, b| {
            a.staar
                .staar_o
                .partial_cmp(&b.staar.staar_o)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let nr = sorted_results.len();
        let mut b_ensembl = StringBuilder::with_capacity(nr, nr * 16);
        let mut b_symbol = StringBuilder::with_capacity(nr, nr * 12);
        let mut b_chrom = StringBuilder::with_capacity(nr, nr * 2);
        let mut b_start = UInt32Builder::with_capacity(nr);
        let mut b_end = UInt32Builder::with_capacity(nr);
        let mut b_nvariants = UInt32Builder::with_capacity(nr);
        let mut b_cmac = UInt32Builder::with_capacity(nr);

        // All p-value columns are Float64. f32 truncated 5e-324 (the
        // smallest f64 denormal, which CCT and chisq survival use as a
        // lower bound) to literal 0.0, masking the very underflows the
        // floor exists to prevent.
        let n_p_cols = 6 + 6 * n_channels + 6; // base + annotation + omnibus
        let mut p_builders: Vec<Float64Builder> = (0..n_p_cols)
            .map(|_| Float64Builder::with_capacity(nr))
            .collect();
        let mut b_acat_o = Float64Builder::with_capacity(nr);
        let mut b_staar_o = Float64Builder::with_capacity(nr);

        for r in &sorted_results {
            let s = &r.staar;
            b_ensembl.append_value(&r.ensembl_id);
            b_symbol.append_value(&r.gene_symbol);
            b_chrom.append_value(r.chromosome.label());
            b_start.append_value(r.start);
            b_end.append_value(r.end);
            b_nvariants.append_value(r.n_variants);
            b_cmac.append_value(r.cumulative_mac);

            let mut pi = 0;
            for p in [
                s.burden_1_25,
                s.burden_1_1,
                s.skat_1_25,
                s.skat_1_1,
                s.acat_v_1_25,
                s.acat_v_1_1,
            ] {
                p_builders[pi].append_value(p);
                pi += 1;
            }
            for ann_p in &s.per_annotation {
                for &v in ann_p {
                    p_builders[pi].append_value(v);
                    pi += 1;
                }
            }
            while pi < 6 + 6 * n_channels {
                p_builders[pi].append_value(f64::NAN);
                pi += 1;
            }
            for p in [
                s.staar_b_1_25,
                s.staar_b_1_1,
                s.staar_s_1_25,
                s.staar_s_1_1,
                s.staar_a_1_25,
                s.staar_a_1_1,
            ] {
                p_builders[pi].append_value(p);
                pi += 1;
            }
            b_acat_o.append_value(s.acat_o);
            b_staar_o.append_value(s.staar_o);
        }

        let test_names = [
            "Burden(1,25)",
            "Burden(1,1)",
            "SKAT(1,25)",
            "SKAT(1,1)",
            "ACAT-V(1,25)",
            "ACAT-V(1,1)",
        ];
        let mut fields = vec![
            Field::new("ensembl_id", DataType::Utf8, false),
            Field::new("gene_symbol", DataType::Utf8, false),
            Field::new(Col::Chromosome.as_str(), DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("end", DataType::UInt32, false),
            Field::new("n_variants", DataType::UInt32, false),
            Field::new("cMAC", DataType::UInt32, false),
        ];
        for test in &test_names {
            fields.push(Field::new(*test, DataType::Float64, true));
        }
        for ch in &channels {
            for test in &test_names {
                fields.push(Field::new(format!("{test}-{ch}"), DataType::Float64, true));
            }
        }
        for name in [
            "STAAR-B(1,25)",
            "STAAR-B(1,1)",
            "STAAR-S(1,25)",
            "STAAR-S(1,1)",
            "STAAR-A(1,25)",
            "STAAR-A(1,1)",
        ] {
            fields.push(Field::new(name, DataType::Float64, true));
        }
        fields.push(Field::new("ACAT-O", DataType::Float64, true));
        fields.push(Field::new("STAAR-O", DataType::Float64, true));
        let schema = Arc::new(Schema::new(fields));

        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(b_ensembl.finish()),
            Arc::new(b_symbol.finish()),
            Arc::new(b_chrom.finish()),
            Arc::new(b_start.finish()),
            Arc::new(b_end.finish()),
            Arc::new(b_nvariants.finish()),
            Arc::new(b_cmac.finish()),
        ];
        for b in &mut p_builders {
            columns.push(Arc::new(b.finish()));
        }
        columns.push(Arc::new(b_acat_o.finish()));
        columns.push(Arc::new(b_staar_o.finish()));

        let batch = RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))?;
        write_parquet_atomic(&out_path, |file| {
            let props = WriterProperties::builder()
                .set_compression(Compression::ZSTD(Default::default()))
                .build();
            let mut writer = ArrowWriter::try_new(file, schema, Some(props))
                .map_err(|e| CohortError::Resource(format!("Parquet writer: {e}")))?;
            writer
                .write(&batch)
                .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
            writer
                .close()
                .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
            Ok(())
        })?;

        let n_sig = results.iter().filter(|r| r.staar.staar_o < 2.5e-6).count();
        for r in results.iter().filter(|r| r.staar.staar_o < 2.5e-6) {
            significant_genes.push(json!({
                "gene": r.ensembl_id, "mask": mask_type.file_stem(),
                "STAAR-O": r.staar.staar_o, "n_variants": r.n_variants,
            }));
        }
        out.success(&format!(
            "  {} -> {} genes, {} significant",
            mask_type.file_stem(),
            results.len(),
            n_sig
        ));
    }

    let meta = json!({
        "cohort_staar_version": 1,
        "traits": trait_names, "trait_type": format!("{:?}", trait_type),
        "n_samples": n, "n_rare_variants": n_rare, "maf_cutoff": maf_cutoff,
        "sigma2": null_model.sigma2, "significant_genes": significant_genes,
    });
    let meta_bytes = serde_json::to_string_pretty(&meta)
        .map_err(|e| CohortError::Resource(format!("Serialize staar.meta.json: {e}")))?;
    write_atomic(&output_dir.join("staar.meta.json"), meta_bytes.as_bytes())?;

    match generate_report(
        all_mask_results,
        trait_names,
        n,
        n_rare,
        output_dir,
        "STAAR Rare Variant Association",
    ) {
        Ok(()) => out.success(&format!(
            "  summary.html -> {}",
            output_dir.join("summary.html").display()
        )),
        Err(e) => out.warn(&format!("  Summary report failed: {e}")),
    }

    out.success(&format!("STAAR complete -> {}", output_dir.display()));
    out.result_json(&meta);
    Ok(())
}

/// Chromosome sizes (hg38, bp) for Manhattan plot layout.
const CHROMS: &[(&str, u64)] = &[
    ("1", 249_000_000),
    ("2", 243_000_000),
    ("3", 198_000_000),
    ("4", 190_000_000),
    ("5", 182_000_000),
    ("6", 171_000_000),
    ("7", 159_000_000),
    ("8", 145_000_000),
    ("9", 138_000_000),
    ("10", 134_000_000),
    ("11", 135_000_000),
    ("12", 133_000_000),
    ("13", 114_000_000),
    ("14", 107_000_000),
    ("15", 102_000_000),
    ("16", 90_000_000),
    ("17", 84_000_000),
    ("18", 80_000_000),
    ("19", 59_000_000),
    ("20", 65_000_000),
    ("21", 47_000_000),
    ("22", 51_000_000),
    ("X", 156_000_000),
    ("Y", 57_000_000),
];

const SIGNIFICANCE: f64 = 2.5e-6;

struct PlotGene {
    gene: String,
    mask: String,
    chromosome: String,
    genome_x: f64,
    neg_log_p: f64,
    n_variants: u32,
    cmac: u32,
    burden: f64,
    skat: f64,
    acat_v: f64,
    staar_o: f64,
    chrom_idx: usize,
}

pub fn generate_report(
    results: &[(MaskType, Vec<GeneResult>)],
    traits: &[String],
    n_samples: usize,
    n_rare: i64,
    output_dir: &Path,
    title: &str,
) -> Result<(), CohortError> {
    let genes = collect_plot_genes(results);
    let pvals = collect_pvalues(results);
    let n_genes: usize = results.iter().map(|(_, r)| r.len()).sum();
    let n_sig = genes.iter().filter(|g| g.staar_o < SIGNIFICANCE).count();
    let masks: Vec<String> = results
        .iter()
        .filter(|(_, r)| !r.is_empty())
        .map(|(m, _)| m.file_stem())
        .collect();

    // Manhattan data: JSON arrays
    let manhattan_json = manhattan_traces(&genes);
    // QQ data
    let qq_json = qq_trace(&pvals);
    // Volcano data
    let volcano_json = volcano_trace(&genes);
    // Table rows
    let table_rows = table_json(&genes);

    let html = format!(
        r####"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
       color:#1a1a2e; background:#fafafa; padding:1.5rem; max-width:1200px; margin:0 auto; }}
h1 {{ font-size:1.4rem; font-weight:600; }}
.meta {{ color:#666; font-size:0.82rem; margin-bottom:1rem; }}
.cards {{ display:flex; gap:1rem; margin-bottom:1.5rem; flex-wrap:wrap; }}
.card {{ background:#fff; border:1px solid #e0e0e0; border-radius:6px; padding:0.6rem 1rem; }}
.card .v {{ font-size:1.3rem; font-weight:700; color:#16213e; }}
.card .l {{ font-size:0.7rem; color:#888; text-transform:uppercase; letter-spacing:0.04em; }}
.sig {{ color:#c62828; }}
.plot-box {{ background:#fff; border:1px solid #e0e0e0; border-radius:6px; padding:0.5rem; margin-bottom:1rem; }}
.row {{ display:flex; gap:1rem; flex-wrap:wrap; }}
.row > .plot-box {{ flex:1; min-width:380px; }}
table {{ width:100%; border-collapse:collapse; font-size:0.78rem; }}
th {{ text-align:left; padding:0.4rem 0.6rem; border-bottom:2px solid #16213e;
     font-weight:600; font-size:0.68rem; text-transform:uppercase; letter-spacing:0.04em; color:#555; cursor:pointer; }}
th:hover {{ color:#16213e; }}
td {{ padding:0.35rem 0.6rem; border-bottom:1px solid #eee; }}
tr:hover td {{ background:#f5f5ff; }}
.p {{ font-family:monospace; }}
input {{ width:100%; padding:0.4rem; margin-bottom:0.5rem; border:1px solid #ddd; border-radius:4px; font-size:0.8rem; }}
footer {{ margin-top:1.5rem; color:#aaa; font-size:0.65rem; text-align:center; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="meta">Traits: {traits} | Masks: {masks}</p>
<div class="cards">
  <div class="card"><div class="v">{n_samples}</div><div class="l">Samples</div></div>
  <div class="card"><div class="v">{n_rare}</div><div class="l">Rare variants</div></div>
  <div class="card"><div class="v">{n_genes}</div><div class="l">Gene-sets tested</div></div>
  <div class="card"><div class="v sig">{n_sig}</div><div class="l">Significant</div></div>
</div>

<div class="plot-box"><div id="manhattan" style="height:350px;"></div></div>
<div class="row">
  <div class="plot-box"><div id="qq" style="height:380px;"></div></div>
  <div class="plot-box"><div id="volcano" style="height:380px;"></div></div>
</div>

<div class="plot-box" style="overflow-x:auto;">
<h3 style="font-size:0.82rem;margin:0.3rem 0;">Gene Results</h3>
<input type="text" id="filter" placeholder="Filter by gene name..." oninput="filterTable()">
<table id="results">
<thead><tr>
<th onclick="sortTable(0)">Gene</th><th onclick="sortTable(1)">Mask</th><th onclick="sortTable(2)">Chr</th>
<th onclick="sortTable(3)">Variants</th><th onclick="sortTable(4)">cMAC</th>
<th onclick="sortTable(5)">STAAR-O</th><th onclick="sortTable(6)">Burden</th>
<th onclick="sortTable(7)">SKAT</th><th onclick="sortTable(8)">ACAT-V</th>
</tr></thead>
<tbody id="tbody">{table_rows}</tbody>
</table>
</div>

<footer>Generated by COHORT CLI | Plotly.js interactive report</footer>

<script>
// Manhattan
(function() {{
  {manhattan_json}
  var layout = {{
    title: {{text:'Manhattan Plot — STAAR-O', font:{{size:13}}}},
    xaxis: {{title:'Genomic position', showgrid:false, tickvals:tickvals, ticktext:ticktext, tickangle:-45, tickfont:{{size:9}}}},
    yaxis: {{title:'-log10(p)'}},
    shapes: [{{type:'line',x0:0,x1:1,xref:'paper',y0:{sig_line},y1:{sig_line},line:{{color:'#c62828',width:1,dash:'dash'}}}},
             {{type:'line',x0:0,x1:1,xref:'paper',y0:{sug_line},y1:{sug_line},line:{{color:'#ff9800',width:0.8,dash:'dot'}}}}],
    hovermode:'closest', margin:{{l:50,r:20,t:40,b:60}}, plot_bgcolor:'#fff',
  }};
  Plotly.newPlot('manhattan', traces, layout, {{responsive:true}});
}})();

// QQ
(function() {{
  {qq_json}
  var layout = {{
    title: {{text:'QQ Plot', font:{{size:13}}}},
    xaxis: {{title:'Expected -log10(p)'}},
    yaxis: {{title:'Observed -log10(p)'}},
    hovermode:'closest', margin:{{l:50,r:20,t:40,b:50}}, plot_bgcolor:'#fff',
  }};
  Plotly.newPlot('qq', traces, layout, {{responsive:true}});
}})();

// Volcano
(function() {{
  {volcano_json}
  var layout = {{
    title: {{text:'Volcano: cMAC vs -log10(p)', font:{{size:13}}}},
    xaxis: {{title:'log10(cMAC)', type:'log'}},
    yaxis: {{title:'-log10(STAAR-O)'}},
    shapes: [{{type:'line',x0:0,x1:1,xref:'paper',y0:{sig_line},y1:{sig_line},line:{{color:'#c62828',width:1,dash:'dash'}}}}],
    hovermode:'closest', margin:{{l:50,r:20,t:40,b:50}}, plot_bgcolor:'#fff',
  }};
  Plotly.newPlot('volcano', traces, layout, {{responsive:true}});
}})();

// Table sort + filter
var sortDir = {{}};
function sortTable(col) {{
  var tb = document.getElementById('tbody');
  var rows = Array.from(tb.rows);
  var dir = sortDir[col] = !(sortDir[col] || false);
  rows.sort(function(a,b) {{
    var av = a.cells[col].getAttribute('data-v') || a.cells[col].textContent;
    var bv = b.cells[col].getAttribute('data-v') || b.cells[col].textContent;
    var an = parseFloat(av), bn = parseFloat(bv);
    if (!isNaN(an) && !isNaN(bn)) return dir ? an-bn : bn-an;
    return dir ? av.localeCompare(bv) : bv.localeCompare(av);
  }});
  rows.forEach(function(r) {{ tb.appendChild(r); }});
}}
function filterTable() {{
  var q = document.getElementById('filter').value.toLowerCase();
  var rows = document.getElementById('tbody').rows;
  for (var i=0; i<rows.length; i++) {{
    rows[i].style.display = rows[i].cells[0].textContent.toLowerCase().includes(q) ? '' : 'none';
  }}
}}
</script>
</body>
</html>"####,
        title = title,
        traits = traits.join(", "),
        masks = masks.join(", "),
        sig_line = -SIGNIFICANCE.log10(),
        sug_line = -(5e-4_f64).log10(),
    );

    write_atomic(&output_dir.join("summary.html"), html.as_bytes())?;
    Ok(())
}

fn collect_plot_genes(results: &[(MaskType, Vec<GeneResult>)]) -> Vec<PlotGene> {
    let offsets = chrom_offsets();
    let mut genes = Vec::new();
    for (mask, gene_results) in results {
        for g in gene_results {
            let raw = g.staar.staar_o;
            if !raw.is_finite() || raw > 1.0 {
                continue;
            }
            let p = if raw <= 0.0 { 1e-300 } else { raw };
            let chrom_label = g.chromosome.label();
            let (idx, offset) = match offsets.get(&chrom_label) {
                Some(v) => *v,
                None => continue,
            };
            let midpoint = (g.start as f64 + g.end as f64) / 2.0;
            genes.push(PlotGene {
                gene: g.gene_symbol.clone(),
                mask: mask.file_stem(),
                chromosome: chrom_label,
                genome_x: offset as f64 + midpoint,
                neg_log_p: -p.log10(),
                n_variants: g.n_variants,
                cmac: g.cumulative_mac,
                burden: g.staar.burden_1_25,
                skat: g.staar.skat_1_25,
                acat_v: g.staar.acat_v_1_25,
                staar_o: p,
                chrom_idx: idx,
            });
        }
    }
    genes
}

fn chrom_offsets() -> HashMap<String, (usize, u64)> {
    let mut map = HashMap::new();
    let mut offset = 0u64;
    for (i, &(name, size)) in CHROMS.iter().enumerate() {
        map.insert(name.to_string(), (i, offset));
        offset += size;
    }
    map
}

fn collect_pvalues(results: &[(MaskType, Vec<GeneResult>)]) -> Vec<f64> {
    results
        .iter()
        .flat_map(|(_, genes)| genes.iter().map(|g| g.staar.staar_o))
        .filter(|p| p.is_finite() && *p <= 1.0)
        .map(|p| if p <= 0.0 { 1e-300 } else { p })
        .collect()
}

fn manhattan_traces(genes: &[PlotGene]) -> String {
    // Group by chromosome for alternating colors
    let colors = ["#4361ee", "#7f8fa6"];
    let mut traces = Vec::new();

    // Deduplicate genes across masks: keep best p-value per gene
    let mut best: HashMap<String, &PlotGene> = HashMap::new();
    for g in genes {
        let key = format!("{}:{}", g.gene, g.chromosome);
        match best.get(&key) {
            Some(existing) if existing.staar_o <= g.staar_o => {}
            _ => {
                best.insert(key, g);
            }
        }
    }
    let mut unique: Vec<&&PlotGene> = best.values().collect();
    unique.sort_by(|a, b| {
        a.genome_x
            .partial_cmp(&b.genome_x)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Split into even/odd chromosomes for color alternation
    for (color_idx, color) in colors.iter().enumerate() {
        let subset: Vec<&&PlotGene> = unique
            .iter()
            .filter(|g| g.chrom_idx % 2 == color_idx)
            .copied()
            .collect();
        if subset.is_empty() {
            continue;
        }

        let xs: Vec<String> = subset
            .iter()
            .map(|g| format!("{:.0}", g.genome_x))
            .collect();
        let ys: Vec<String> = subset
            .iter()
            .map(|g| format!("{:.4}", g.neg_log_p))
            .collect();
        let texts: Vec<String> = subset
            .iter()
            .map(|g| {
                format!(
                    "{}<br>chr{}<br>p={:.2e}<br>mask: {}<br>variants: {}<br>cMAC: {}",
                    g.gene, g.chromosome, g.staar_o, g.mask, g.n_variants, g.cmac
                )
            })
            .collect();

        traces.push(format!(
            "{{x:[{xs}],y:[{ys}],text:[{texts}],mode:'markers',type:'scatter',\
             marker:{{color:'{color}',size:5,opacity:0.6}},\
             hoverinfo:'text',showlegend:false}}",
            xs = xs.join(","),
            ys = ys.join(","),
            texts = texts
                .iter()
                .map(|t| format!("'{}'", t.replace('\'', "\\'")))
                .collect::<Vec<_>>()
                .join(","),
            color = color,
        ));
    }

    // Chromosome tick positions + labels
    let offsets = chrom_offsets();
    let mut tickvals = Vec::new();
    let mut ticktext = Vec::new();
    for &(name, size) in CHROMS {
        if let Some(&(_, off)) = offsets.get(name) {
            tickvals.push(format!("{}", off + size / 2));
            ticktext.push(format!("'{name}'"));
        }
    }

    format!(
        "var traces=[{}];\nvar tickvals=[{}];\nvar ticktext=[{}];",
        traces.join(","),
        tickvals.join(","),
        ticktext.join(","),
    )
}

fn qq_trace(pvalues: &[f64]) -> String {
    let mut valid: Vec<f64> = pvalues
        .iter()
        .filter(|p| p.is_finite() && **p > 0.0 && **p <= 1.0)
        .copied()
        .collect();
    valid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    if valid.is_empty() {
        return "var traces=[];".into();
    }

    let n = valid.len();
    let observed: Vec<f64> = valid.iter().map(|p| -p.log10()).collect();
    let expected: Vec<f64> = (0..n)
        .map(|i| -((i as f64 + 0.5) / n as f64).log10())
        .collect();
    let axis_max = observed
        .last()
        .copied()
        .unwrap_or(1.0)
        .max(*expected.last().unwrap_or(&1.0))
        + 0.5;

    // Subsample for performance, keep tail
    let step = if n > 3000 { n / 3000 } else { 1 };
    let mut exp_s = Vec::new();
    let mut obs_s = Vec::new();
    for i in (0..n).step_by(step) {
        exp_s.push(format!("{:.4}", expected[i]));
        obs_s.push(format!("{:.4}", observed[i]));
    }
    // Always include tail
    for i in n.saturating_sub(100)..n {
        exp_s.push(format!("{:.4}", expected[i]));
        obs_s.push(format!("{:.4}", observed[i]));
    }

    format!(
        "var traces=[{{x:[{exp}],y:[{obs}],mode:'markers',type:'scatter',\
         marker:{{color:'#4361ee',size:3,opacity:0.5}},hoverinfo:'x+y',showlegend:false}},\
         {{x:[0,{ax:.1}],y:[0,{ax:.1}],mode:'lines',line:{{color:'#ccc',dash:'dash'}},showlegend:false}}];",
        exp = exp_s.join(","),
        obs = obs_s.join(","),
        ax = axis_max,
    )
}

fn volcano_trace(genes: &[PlotGene]) -> String {
    // Deduplicate: best p per gene
    let mut best: HashMap<String, &PlotGene> = HashMap::new();
    for g in genes {
        let key = format!("{}:{}", g.gene, g.chromosome);
        match best.get(&key) {
            Some(existing) if existing.staar_o <= g.staar_o => {}
            _ => {
                best.insert(key, g);
            }
        }
    }

    let xs: Vec<String> = best
        .values()
        .map(|g| format!("{}", g.cmac.max(1)))
        .collect();
    let ys: Vec<String> = best
        .values()
        .map(|g| format!("{:.4}", g.neg_log_p))
        .collect();
    let texts: Vec<String> = best
        .values()
        .map(|g| {
            format!(
                "{}<br>chr{}<br>p={:.2e}<br>cMAC: {}<br>variants: {}",
                g.gene, g.chromosome, g.staar_o, g.cmac, g.n_variants
            )
        })
        .collect();
    let colors: Vec<String> = best
        .values()
        .map(|g| {
            if g.staar_o < SIGNIFICANCE {
                "'#c62828'".into()
            } else {
                "'#4361ee'".into()
            }
        })
        .collect();

    format!(
        "var traces=[{{x:[{xs}],y:[{ys}],text:[{texts}],mode:'markers',type:'scatter',\
         marker:{{color:[{colors}],size:5,opacity:0.6}},hoverinfo:'text',showlegend:false}}];",
        xs = xs.join(","),
        ys = ys.join(","),
        texts = texts
            .iter()
            .map(|t| format!("'{}'", t.replace('\'', "\\'")))
            .collect::<Vec<_>>()
            .join(","),
        colors = colors.join(","),
    )
}

fn table_json(genes: &[PlotGene]) -> String {
    let mut sorted: Vec<&PlotGene> = genes.iter().collect();
    sorted.sort_by(|a, b| {
        a.staar_o
            .partial_cmp(&b.staar_o)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut html = String::new();
    for g in &sorted {
        let sig_class = if g.staar_o < SIGNIFICANCE {
            " class='sig'"
        } else {
            ""
        };
        html.push_str(&format!(
            "<tr><td><b>{gene}</b></td><td>{mask}</td><td>{chr}</td>\
             <td data-v='{nv}'>{nv}</td><td data-v='{cmac}'>{cmac}</td>\
             <td{sig} data-v='{so_raw}'><span class='p'>{so}</span></td>\
             <td data-v='{b_raw}'><span class='p'>{b}</span></td>\
             <td data-v='{sk_raw}'><span class='p'>{sk}</span></td>\
             <td data-v='{av_raw}'><span class='p'>{av}</span></td></tr>",
            gene = g.gene,
            mask = g.mask,
            chr = g.chromosome,
            nv = g.n_variants,
            cmac = g.cmac,
            sig = sig_class,
            so_raw = g.staar_o,
            so = fmt_p(g.staar_o),
            b_raw = g.burden,
            b = fmt_p(g.burden),
            sk_raw = g.skat,
            sk = fmt_p(g.skat),
            av_raw = g.acat_v,
            av = fmt_p(g.acat_v),
        ));
    }
    html
}

fn fmt_p(p: f64) -> String {
    if p.is_nan() {
        "NaN".into()
    } else if p < 1e-300 {
        "<1e-300".into()
    } else {
        format!("{:.2e}", p)
    }
}

#[cfg(test)]
mod tests {
    use super::super::score::StaarResult;
    use super::*;
    use crate::types::Chromosome;

    fn gene(chr: &str, start: u32, p: f64) -> GeneResult {
        GeneResult {
            ensembl_id: format!("ENSG{start}"),
            gene_symbol: format!("GENE{start}"),
            chromosome: chr.parse::<Chromosome>().unwrap_or(Chromosome::Autosome(1)),
            start,
            end: start + 10000,
            n_variants: 5,
            cumulative_mac: 10,
            staar: StaarResult {
                burden_1_25: p * 2.0,
                burden_1_1: p * 3.0,
                skat_1_25: p * 1.5,
                skat_1_1: p * 2.5,
                acat_v_1_25: p * 1.2,
                acat_v_1_1: p * 2.0,
                per_annotation: Vec::new(),
                staar_b_1_25: p,
                staar_b_1_1: p,
                staar_s_1_25: p,
                staar_s_1_1: p,
                staar_a_1_25: p,
                staar_a_1_1: p,
                acat_o: p,
                staar_o: p,
            },
        }
    }

    #[test]
    fn report_generates_plotly_html() {
        let results = vec![(
            MaskType::PLof,
            vec![
                gene("1", 50_000_000, 1e-8),
                gene("2", 100_000_000, 0.05),
                gene("22", 30_000_000, 1e-3),
            ],
        )];
        let dir = std::env::temp_dir().join("cohort_test_report");
        let _ = std::fs::create_dir_all(&dir);
        generate_report(&results, &["BMI".into()], 1000, 500, &dir, "Test").unwrap();
        let html = std::fs::read_to_string(dir.join("summary.html")).unwrap();
        assert!(html.contains("plotly"));
        assert!(html.contains("Plotly.newPlot"));
        assert!(html.contains("manhattan"));
        assert!(html.contains("qq"));
        assert!(html.contains("volcano"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn top_genes_filters_significant() {
        let results = vec![(
            MaskType::PLof,
            vec![gene("1", 100, 1e-8), gene("1", 200, 0.5)],
        )];
        let genes = collect_plot_genes(&results);
        let sig: Vec<_> = genes.iter().filter(|g| g.staar_o < SIGNIFICANCE).collect();
        assert_eq!(sig.len(), 1);
    }

    #[test]
    fn empty_results_no_panic() {
        let results: Vec<(MaskType, Vec<GeneResult>)> = Vec::new();
        let dir = std::env::temp_dir().join("cohort_test_empty");
        let _ = std::fs::create_dir_all(&dir);
        generate_report(&results, &["BMI".into()], 0, 0, &dir, "Empty").unwrap();
        let html = std::fs::read_to_string(dir.join("summary.html")).unwrap();
        assert!(html.contains("plotly"));
        let _ = std::fs::remove_dir_all(&dir);
    }
}

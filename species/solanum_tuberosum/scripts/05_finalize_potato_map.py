#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd

BASE = Path("species/solanum_tuberosum")

MERGED_FILE = BASE / "results" / "intermediate" / "potato_genetic_map.sharma2013_TableS4_merged_chr_agree.tsv"

FINAL_DIR = BASE / "results" / "final"
QC_DIR = BASE / "results" / "qc"

FINAL_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

def chr_to_number(x):
    x = str(x)
    m = re.search(r"(\d+)", x)
    if not m:
        return x
    return str(int(m.group(1)))

def spearman_without_scipy(x, y):
    """
    Spearman correlation = Pearson correlation between ranks.
    This avoids scipy dependency.
    """
    xr = pd.Series(x).rank(method="average")
    yr = pd.Series(y).rank(method="average")
    return xr.corr(yr, method="pearson")

df = pd.read_csv(MERGED_FILE, sep="\t")

# Оставляем только то, что нужно для финальной карты.
final = df[["chr", "pos", "cM"]].copy()
final["chr"] = final["chr"].map(chr_to_number)
final["pos"] = pd.to_numeric(final["pos"], errors="coerce").astype("Int64")
final["cM"] = pd.to_numeric(final["cM"], errors="coerce")

final = final.dropna(subset=["chr", "pos", "cM"])
final = final.drop_duplicates()

# Сортировка: chr как число, потом физическая позиция.
final["_chr_sort"] = final["chr"].astype(int)
final = final.sort_values(["_chr_sort", "pos", "cM"]).drop(columns=["_chr_sort"])

final_file = FINAL_DIR / "potato_genetic_map.tsv"
final.to_csv(final_file, sep="\t", index=False)

# Более подробная таблица для трассировки: какие marker_id вошли в финальную карту.
trace = df.copy()
trace["chr_final"] = trace["chr"].map(chr_to_number)
trace = trace[
    [
        "marker_id_original",
        "marker_id",
        "genetic_chr",
        "chr",
        "chr_final",
        "pos",
        "cM",
        "num_mismatch",
        "multimap",
        "source",
        "marker_type",
        "tableS4_seqid_v403",
        "tableS4_start_v403",
        "tableS4_end_v403",
    ]
].sort_values(["chr_final", "pos", "cM"])

trace_file = FINAL_DIR / "potato_genetic_map.with_markers.tsv"
trace.to_csv(trace_file, sep="\t", index=False)

# QC по хромосомам.
qc_rows = []
for chrom, sub in final.groupby("chr", sort=False):
    sub = sub.sort_values("pos")
    n = len(sub)

    if n >= 2:
        spearman = spearman_without_scipy(sub["pos"], sub["cM"])
        cM_diff = sub["cM"].diff().dropna()
        increasing_steps = int((cM_diff > 0).sum())
        decreasing_steps = int((cM_diff < 0).sum())
        equal_steps = int((cM_diff == 0).sum())
    else:
        spearman = float("nan")
        increasing_steps = decreasing_steps = equal_steps = 0

    qc_rows.append(
        {
            "chr": chrom,
            "n_markers": n,
            "pos_min": int(sub["pos"].min()),
            "pos_max": int(sub["pos"].max()),
            "cM_min": float(sub["cM"].min()),
            "cM_max": float(sub["cM"].max()),
            "spearman_pos_cM": spearman,
            "increasing_cM_steps": increasing_steps,
            "decreasing_cM_steps": decreasing_steps,
            "equal_cM_steps": equal_steps,
        }
    )

qc = pd.DataFrame(qc_rows)
qc_file = QC_DIR / "potato_final_map_qc_by_chr.tsv"
qc.to_csv(qc_file, sep="\t", index=False)

summary = []
summary.append("Potato final genetic map")
summary.append("=" * 35)
summary.append("")
summary.append("Target reference: DM_1-3_516_R44_potato.v6.1")
summary.append("Genetic source: Sharma et al. 2013 TableS4")
summary.append("Physical source: SpudDB SolCAP 69K SNP positions on DM v6.1")
summary.append("")
summary.append(f"Final chr-pos-cM rows: {len(final)}")
summary.append(f"Chromosomes: {', '.join(final['chr'].drop_duplicates().astype(str))}")
summary.append(f"Final map file: {final_file}")
summary.append(f"Trace file with marker IDs: {trace_file}")
summary.append(f"QC by chromosome: {qc_file}")
summary.append("")
summary.append("Markers per chromosome:")
summary.append(final.groupby("chr").size().to_string())
summary.append("")
summary.append("Mismatch distribution for included markers:")
summary.append(trace.groupby("num_mismatch").size().to_string())
summary.append("")

summary_file = QC_DIR / "potato_final_map_summary.txt"
summary_file.write_text("\n".join(summary) + "\n")

print("Done.")
print(f"Final map:              {final_file}")
print(f"Trace map with markers: {trace_file}")
print(f"QC by chromosome:       {qc_file}")
print(f"Summary:                {summary_file}")

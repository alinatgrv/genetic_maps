#!/usr/bin/env python3

# QC коллинеарности финальной карты рапса: для каждой хромосомы считаем
# Spearman-корреляцию между физической позицией (pos) и генетической (cM),
# а также число возрастающих/убывающих шагов cM при сортировке по позиции.
#
# Дополнительно: общая сводка карты и проверка дублирующихся (chr, pos).
#
# Запускать из корня репозитория.

from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SPECIES_DIR = SCRIPT_DIR.parent

IN_MAP = SPECIES_DIR / "results/final/brassica_napus_genetic_map.tsv"
QC_DIR = SPECIES_DIR / "results/qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

OUT_TSV = QC_DIR / "brassica_napus_final_map_collinearity_qc.tsv"
OUT_DUP = QC_DIR / "brassica_napus_final_map_duplicate_positions.tsv"
OUT_SUMMARY = QC_DIR / "brassica_napus_final_map_summary.txt"

CHR_ORDER = [f"A{i}" for i in range(1, 11)] + [f"C{i}" for i in range(1, 10)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_ORDER)}

df = pd.read_csv(IN_MAP, sep="\t")
df["chr"] = df["chr"].astype(str)
df["pos"] = pd.to_numeric(df["pos"])
df["cM"] = pd.to_numeric(df["cM"])

rows = []
for chrom, g in df.groupby("chr"):
    g = g.sort_values(["pos", "cM"]).copy()
    n = len(g)

    if n >= 3:
        spearman = g["pos"].rank(method="average").corr(
            g["cM"].rank(method="average"), method="pearson")
    else:
        spearman = float("nan")

    cm_diff = g["cM"].diff().dropna()
    increasing_steps = int((cm_diff > 0).sum())
    decreasing_steps = int((cm_diff < 0).sum())
    equal_steps = int((cm_diff == 0).sum())
    informative = increasing_steps + decreasing_steps
    increasing_fraction = increasing_steps / informative if informative else float("nan")

    if pd.isna(spearman):
        orientation = "NA"
    elif spearman > 0:
        orientation = "increasing"
    elif spearman < 0:
        orientation = "decreasing"
    else:
        orientation = "flat_or_mixed"

    rows.append({
        "chr": chrom,
        "n_markers": n,
        "pos_min": int(g["pos"].min()),
        "pos_max": int(g["pos"].max()),
        "physical_span_Mb": round((g["pos"].max() - g["pos"].min()) / 1_000_000, 3),
        "cM_min": round(g["cM"].min(), 3),
        "cM_max": round(g["cM"].max(), 3),
        "genetic_span_cM": round(g["cM"].max() - g["cM"].min(), 3),
        "spearman_pos_cM": round(spearman, 4) if pd.notna(spearman) else "NA",
        "abs_spearman": round(abs(spearman), 4) if pd.notna(spearman) else "NA",
        "orientation": orientation,
        "increasing_cM_steps": increasing_steps,
        "decreasing_cM_steps": decreasing_steps,
        "equal_cM_steps": equal_steps,
        "monotonic_fraction": round(max(increasing_fraction, 1 - increasing_fraction), 3) if informative else "NA",
    })

qc = pd.DataFrame(rows)
qc["chr_idx"] = qc["chr"].map(CHR_INDEX)
qc = qc.sort_values("chr_idx").drop(columns=["chr_idx"])
qc.to_csv(OUT_TSV, sep="\t", index=False)

# Дублирующиеся (chr, pos).
dup = df[df.duplicated(subset=["chr", "pos"], keep=False)].copy()
dup["chr_idx"] = dup["chr"].map(CHR_INDEX)
dup = dup.sort_values(["chr_idx", "pos"]).drop(columns=["chr_idx"])
dup.to_csv(OUT_DUP, sep="\t", index=False)

# Общая сводка. abs_spearman усредняем по хромосомам, где он определён.
abs_vals = [r["abs_spearman"] for r in rows if r["abs_spearman"] != "NA"]
mean_abs = sum(abs_vals) / len(abs_vals) if abs_vals else float("nan")
median_abs = sorted(abs_vals)[len(abs_vals) // 2] if abs_vals else float("nan")

with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
    out.write("Brassica napus final genetic map — QC summary\n")
    out.write("=" * 47 + "\n\n")
    out.write(f"Input map: {IN_MAP}\n")
    out.write(f"Total markers: {len(df)}\n")
    out.write(f"Chromosomes covered: {df['chr'].nunique()} / 19\n")
    out.write(f"Duplicate (chr,pos) rows: {len(dup)}\n\n")
    out.write(f"Mean |Spearman(pos,cM)| over chromosomes:   {mean_abs:.4f}\n")
    out.write(f"Median |Spearman(pos,cM)| over chromosomes: {median_abs:.4f}\n\n")
    out.write("Per-chromosome |Spearman| (n markers):\n")
    for r in rows:
        order_idx = CHR_INDEX[r["chr"]]
        out.write(f"  {r['chr']:>3}\t{r['abs_spearman']}\t(n={r['n_markers']})\n")

print("Done.")
print(f"Collinearity QC: {OUT_TSV}")
print(f"Duplicate positions: {OUT_DUP} ({len(dup)} rows)")
print(f"Summary: {OUT_SUMMARY}")
print()
print(qc.to_string(index=False))
print()
print(f"Mean |Spearman| = {mean_abs:.4f} ; Median |Spearman| = {median_abs:.4f}")

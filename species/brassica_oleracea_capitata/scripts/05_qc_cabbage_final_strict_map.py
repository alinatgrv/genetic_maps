import pandas as pd
from pathlib import Path

input_file = Path("results/cabbage/cabbage_genetic_map.strict.clean.unique.tsv")
out_file = Path("results/cabbage/cabbage_genetic_map.strict.clean.unique.qc.tsv")

df = pd.read_csv(input_file, sep="\t")

rows = []

for chrom, g in df.groupby("chr"):
    g = g.sort_values("pos").copy()

    n = len(g)
    pos_min = g["pos"].min()
    pos_max = g["pos"].max()
    cm_min = g["cM"].min()
    cm_max = g["cM"].max()

    # Correlation between physical and genetic position.
    # Positive: cM mostly increases with pos.
    # Negative: cM mostly decreases with pos.
    if n >= 3:
        spearman = g["pos"].rank().corr(g["cM"].rank())
    else:
        spearman = float("nan")

    # Count local direction changes.
    cm_diff = g["cM"].diff().dropna()
    increasing_steps = (cm_diff > 0).sum()
    decreasing_steps = (cm_diff < 0).sum()
    equal_steps = (cm_diff == 0).sum()

    rows.append({
        "chr": chrom,
        "n_markers": n,
        "pos_min": pos_min,
        "pos_max": pos_max,
        "cM_min": cm_min,
        "cM_max": cm_max,
        "spearman_pos_cM": spearman,
        "increasing_cM_steps": int(increasing_steps),
        "decreasing_cM_steps": int(decreasing_steps),
        "equal_cM_steps": int(equal_steps),
    })

qc = pd.DataFrame(rows)
qc["chr_num"] = qc["chr"].astype(int)
qc = qc.sort_values("chr_num").drop(columns=["chr_num"])

qc.to_csv(out_file, sep="\t", index=False)

print("QC table written to:", out_file)
print(qc.to_string(index=False))

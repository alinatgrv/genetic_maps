from pathlib import Path
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
SPECIES_DIR = SCRIPT_DIR.parent

IN_MAP = SPECIES_DIR / "results/final/sunflower_genetic_map.bwa_exact_unique.tsv"
QC_DIR = SPECIES_DIR / "results/qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

OUT_TSV = QC_DIR / "sunflower_final_map_collinearity_qc.tsv"

df = pd.read_csv(IN_MAP, sep="\t")

df["chr"] = df["chr"].astype(str)
df["pos"] = pd.to_numeric(df["pos"])
df["cM"] = pd.to_numeric(df["cM"])

rows = []

for chrom, g in df.groupby("chr"):
    g = g.sort_values(["pos", "cM"]).copy()
    n = len(g)

    pos_min = g["pos"].min()
    pos_max = g["pos"].max()
    cm_min = g["cM"].min()
    cm_max = g["cM"].max()

    if n >= 3:
        spearman = g["pos"].rank(method="average").corr(
            g["cM"].rank(method="average"),
            method="pearson"
        )
    else:
        spearman = float("nan")

    cm_diff = g["cM"].diff().dropna()

    increasing_steps = int((cm_diff > 0).sum())
    decreasing_steps = int((cm_diff < 0).sum())
    equal_steps = int((cm_diff == 0).sum())

    informative_steps = increasing_steps + decreasing_steps

    if informative_steps > 0:
        increasing_fraction = increasing_steps / informative_steps
        decreasing_fraction = decreasing_steps / informative_steps
    else:
        increasing_fraction = float("nan")
        decreasing_fraction = float("nan")

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
        "pos_min": int(pos_min),
        "pos_max": int(pos_max),
        "cM_min": cm_min,
        "cM_max": cm_max,
        "spearman_pos_cM": spearman,
        "orientation": orientation,
        "increasing_cM_steps": increasing_steps,
        "decreasing_cM_steps": decreasing_steps,
        "equal_cM_steps": equal_steps,
        "increasing_fraction_without_ties": increasing_fraction,
        "decreasing_fraction_without_ties": decreasing_fraction,
    })

qc = pd.DataFrame(rows)
qc["chr_num"] = qc["chr"].astype(int)
qc = qc.sort_values("chr_num").drop(columns=["chr_num"])

qc.to_csv(OUT_TSV, sep="\t", index=False)

print("QC table written to:", OUT_TSV)
print()
print(qc.to_string(index=False))

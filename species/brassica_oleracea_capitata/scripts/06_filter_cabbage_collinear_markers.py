import pandas as pd
from pathlib import Path

input_file = Path("results/cabbage/cabbage_genetic_map.strict.clean_with_markers.tsv")

collinear_with_markers_out = Path("results/cabbage/cabbage_genetic_map.strict.clean.collinear_with_markers.tsv")
collinear_final_out = Path("results/cabbage/cabbage_genetic_map.strict.clean.collinear.tsv")
outliers_out = Path("results/cabbage/cabbage_genetic_map.strict.clean.collinearity_outliers.tsv")
summary_out = Path("results/cabbage/cabbage_genetic_map.strict.clean.collinearity.summary.txt")

df = pd.read_csv(input_file, sep="\t")

# Remove duplicate physical positions for this QC step.
# They were already excluded from the unique final map.
df = df[~df.duplicated(subset=["chr", "pos"], keep=False)].copy()

def longest_nondecreasing_subsequence_indices(values):
    """
    Return indices of one longest non-decreasing subsequence.
    Dynamic programming is fine here because each chromosome has few markers.
    """
    n = len(values)
    if n == 0:
        return []

    dp = [1] * n
    prev = [-1] * n

    for i in range(n):
        for j in range(i):
            if values[j] <= values[i] and dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                prev[i] = j

    best_i = max(range(n), key=lambda i: dp[i])

    seq = []
    while best_i != -1:
        seq.append(best_i)
        best_i = prev[best_i]

    return list(reversed(seq))

kept_parts = []
outlier_parts = []
summary_rows = []

for chrom, g in df.groupby("chr"):
    g = g.sort_values(["pos", "cM", "marker_id"]).reset_index(drop=True)

    # Determine physical-genetic orientation.
    if len(g) >= 3:
        spearman = g["pos"].rank().corr(g["cM"].rank())
    else:
        spearman = 0

    orientation = "increasing" if spearman >= 0 else "decreasing"

    # For decreasing orientation, multiply cM by -1 so that we can still search
    # for a non-decreasing sequence.
    if orientation == "increasing":
        values = g["cM"].tolist()
    else:
        values = (-g["cM"]).tolist()

    keep_local_idx = set(longest_nondecreasing_subsequence_indices(values))

    g["orientation"] = orientation
    g["spearman_pos_cM_before"] = spearman
    g["collinear_keep"] = [i in keep_local_idx for i in range(len(g))]

    kept = g[g["collinear_keep"]].copy()
    outliers = g[~g["collinear_keep"]].copy()

    kept_parts.append(kept)
    outlier_parts.append(outliers)

    if len(kept) >= 3:
        spearman_after = kept["pos"].rank().corr(kept["cM"].rank())
    else:
        spearman_after = float("nan")

    summary_rows.append({
        "chr": chrom,
        "orientation": orientation,
        "n_before": len(g),
        "n_kept": len(kept),
        "n_outliers": len(outliers),
        "spearman_before": spearman,
        "spearman_after": spearman_after,
    })

kept_all = pd.concat(kept_parts, ignore_index=True)
outliers_all = pd.concat(outlier_parts, ignore_index=True)
summary = pd.DataFrame(summary_rows)

# Sort outputs
kept_all["chr_order"] = kept_all["chr"].astype(int)
outliers_all["chr_order"] = outliers_all["chr"].astype(int)
summary["chr_order"] = summary["chr"].astype(int)

kept_all = kept_all.sort_values(["chr_order", "pos", "cM", "marker_id"]).drop(columns=["chr_order"])
outliers_all = outliers_all.sort_values(["chr_order", "pos", "cM", "marker_id"]).drop(columns=["chr_order"])
summary = summary.sort_values("chr_order").drop(columns=["chr_order"])

# Final task-format table
final = kept_all[["chr", "pos", "cM"]].copy()
final = final.sort_values(["chr", "pos"])

kept_all.to_csv(collinear_with_markers_out, sep="\t", index=False)
final.to_csv(collinear_final_out, sep="\t", index=False)
outliers_all.to_csv(outliers_out, sep="\t", index=False)
summary.to_csv(summary_out, sep="\t", index=False)

print(f"Collinear map with marker info: {collinear_with_markers_out}")
print(f"Final collinear chr-pos-cM map: {collinear_final_out}")
print(f"Outliers: {outliers_out}")
print(f"Summary: {summary_out}")

print("\nSummary:")
print(summary.to_string(index=False))

print("\nTotal markers before collinearity filter:", len(df))
print("Total markers kept:", len(kept_all))
print("Total outliers:", len(outliers_all))

print("\nFirst rows of final collinear map:")
print(final.head(10).to_string(index=False))

#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


BASE = Path("species/solanum_tuberosum")
FINAL = BASE / "results/final"
QC = BASE / "results/qc"

IN_MAP = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.tsv"
IN_META = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.with_metadata.tsv"

OUT_MAP = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.max_collinear.tsv"
OUT_META = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.max_collinear.with_metadata.tsv"

OUT_REMOVED = QC / "potato_max_collinear_removed_positions.tsv"
OUT_QC = QC / "potato_max_collinear_qc_by_chr.tsv"
OUT_SUMMARY = QC / "potato_max_collinear_summary.txt"


def normalize_chr_series(s):
    return (
        s.astype(str)
        .str.strip()
        .str.replace("chr", "", regex=False)
        .str.replace("Chr", "", regex=False)
        .str.replace("CHR", "", regex=False)
        .astype(int)
    )


def spearman_without_scipy(x, y):
    x = pd.Series(x).rank(method="average")
    y = pd.Series(y).rank(method="average")
    if len(x) < 2:
        return np.nan
    return x.corr(y, method="pearson")


def longest_non_decreasing_indices(values):
    """
    O(n^2) longest non-decreasing subsequence.
    n per chromosome is small here, so this is simple and transparent.
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

    best = max(range(n), key=lambda i: dp[i])

    indices = []
    while best != -1:
        indices.append(best)
        best = prev[best]

    return list(reversed(indices))


def collinearity_qc(df):
    rows = []

    for chr_id, g in df.groupby("chr", sort=True):
        g = g.sort_values(["pos", "cM"]).reset_index(drop=True)

        dcM = g["cM"].diff()
        dpos = g["pos"].diff()
        valid = (dpos > 0) & dcM.notna()

        inc = int(((dcM > 0) & valid).sum())
        dec = int(((dcM < 0) & valid).sum())
        eq = int(((dcM == 0) & valid).sum())
        non_tie = inc + dec

        rho = spearman_without_scipy(g["pos"], g["cM"])

        if pd.isna(rho):
            orientation = "not_enough_markers"
        elif rho >= 0.7:
            orientation = "increasing"
        elif rho <= -0.7:
            orientation = "decreasing"
        else:
            orientation = "mixed_or_weak"

        rows.append(
            {
                "chr": chr_id,
                "n_positions": len(g),
                "pos_min": int(g["pos"].min()),
                "pos_max": int(g["pos"].max()),
                "cM_min": float(g["cM"].min()),
                "cM_max": float(g["cM"].max()),
                "spearman_pos_cM": rho,
                "orientation": orientation,
                "increasing_cM_steps": inc,
                "decreasing_cM_steps": dec,
                "equal_cM_steps": eq,
                "increasing_fraction_without_ties": inc / non_tie if non_tie else np.nan,
                "decreasing_fraction_without_ties": dec / non_tie if non_tie else np.nan,
            }
        )

    return pd.DataFrame(rows)


def main():
    if not IN_MAP.exists():
        raise SystemExit(f"ERROR: input map not found: {IN_MAP}")

    df = pd.read_csv(IN_MAP, sep="\t")
    df["chr"] = normalize_chr_series(df["chr"])
    df["pos"] = df["pos"].astype(int)
    df["cM"] = df["cM"].astype(float)

    kept_parts = []
    removed_parts = []
    chr_decisions = []

    for chr_id, g in df.groupby("chr", sort=True):
        g = g.sort_values(["pos", "cM"]).reset_index(drop=True)

        rho = spearman_without_scipy(g["pos"], g["cM"])

        # Orientation decision:
        # chr12 is expected to be reversed. More generally, use Spearman sign.
        if rho < 0:
            orientation = "decreasing"
            values_for_lis = (-g["cM"]).tolist()
        else:
            orientation = "increasing"
            values_for_lis = g["cM"].tolist()

        keep_local_idx = set(longest_non_decreasing_indices(values_for_lis))

        kept = g.loc[sorted(keep_local_idx)].copy()
        removed = g.loc[[i for i in range(len(g)) if i not in keep_local_idx]].copy()

        kept["selected_orientation"] = orientation
        removed["selected_orientation"] = orientation
        removed["reason_removed"] = "breaks_monotonic_pos_cM_chain"

        kept_parts.append(kept)
        removed_parts.append(removed)

        chr_decisions.append(
            {
                "chr": chr_id,
                "input_positions": len(g),
                "kept_positions": len(kept),
                "removed_positions": len(removed),
                "input_spearman": rho,
                "selected_orientation": orientation,
            }
        )

    out = pd.concat(kept_parts, ignore_index=True)
    out = out.sort_values(["chr", "pos", "cM"])
    out[["chr", "pos", "cM"]].to_csv(OUT_MAP, sep="\t", index=False)

    removed = pd.concat(removed_parts, ignore_index=True)
    removed = removed.sort_values(["chr", "pos", "cM"])
    removed.to_csv(OUT_REMOVED, sep="\t", index=False)

    # Metadata version, if metadata table exists.
    if IN_META.exists():
        meta = pd.read_csv(IN_META, sep="\t")
        meta["chr"] = normalize_chr_series(meta["chr"])
        meta["pos"] = meta["pos"].astype(int)
        meta["cM"] = meta["cM"].astype(float)

        selected_keys = out[["chr", "pos", "cM", "selected_orientation"]].copy()

        out_meta = meta.merge(
            selected_keys,
            on=["chr", "pos", "cM"],
            how="inner",
        )
        out_meta = out_meta.sort_values(["chr", "pos", "cM"])
        out_meta.to_csv(OUT_META, sep="\t", index=False)

    qc = collinearity_qc(out)
    qc.to_csv(OUT_QC, sep="\t", index=False)

    decisions = pd.DataFrame(chr_decisions)

    with open(OUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("Potato maximum-collinearity filtered map\n")
        f.write("========================================\n\n")
        f.write(f"Input map: {IN_MAP}\n")
        f.write(f"Input rows: {len(df)}\n")
        f.write(f"Kept rows: {len(out)}\n")
        f.write(f"Removed rows: {len(removed)}\n\n")

        f.write("Method:\n")
        f.write("- Within each chromosome, rows were sorted by physical position.\n")
        f.write("- For chromosomes with positive Spearman(pos,cM), the longest non-decreasing cM subsequence was retained.\n")
        f.write("- For chromosomes with negative Spearman(pos,cM), the longest non-increasing cM subsequence was retained.\n")
        f.write("- cM values were not modified; only non-collinear positions were filtered out.\n\n")

        f.write("Chromosome-level filtering decisions:\n")
        f.write(decisions.to_string(index=False))
        f.write("\n\n")

        f.write("Collinearity QC after filtering:\n")
        f.write(qc.to_string(index=False))
        f.write("\n\n")

        f.write("Output files:\n")
        f.write(f"- max-collinear map: {OUT_MAP}\n")
        f.write(f"- max-collinear metadata: {OUT_META}\n")
        f.write(f"- removed positions: {OUT_REMOVED}\n")
        f.write(f"- QC by chromosome: {OUT_QC}\n")

    print("Done.")
    print(f"Max-collinear map:      {OUT_MAP}")
    print(f"Max-collinear metadata: {OUT_META}")
    print(f"Removed positions:      {OUT_REMOVED}")
    print(f"QC by chr:              {OUT_QC}")
    print(f"Summary:                {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

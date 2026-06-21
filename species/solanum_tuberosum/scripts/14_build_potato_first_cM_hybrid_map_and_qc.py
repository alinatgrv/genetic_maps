#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


BASE = Path("species/solanum_tuberosum")
FINAL = BASE / "results/final"
QC = BASE / "results/qc"

QC.mkdir(parents=True, exist_ok=True)
FINAL.mkdir(parents=True, exist_ok=True)

IN_WITH_MARKERS = FINAL / "potato_genetic_map.hybrid_high_confidence.with_markers.tsv"

OUT_MAP = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.tsv"
OUT_META = FINAL / "potato_genetic_map.hybrid_high_confidence.first_cM.with_metadata.tsv"

OUT_QC_BY_CHR = QC / "potato_hybrid_high_confidence_first_cM_collinearity_qc_by_chr.tsv"
OUT_DUPLICATES = QC / "potato_hybrid_high_confidence_first_cM_duplicate_positions.tsv"
OUT_SUMMARY = QC / "potato_hybrid_high_confidence_first_cM_summary.txt"


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


def join_unique(x, max_items=20):
    vals = [str(v) for v in sorted(set(x.dropna().astype(str))) if v != ""]
    if len(vals) > max_items:
        return ",".join(vals[:max_items]) + f",...(+{len(vals)-max_items})"
    return ",".join(vals)


def join_values_in_order(x, max_items=30):
    vals = [str(v) for v in x.dropna().tolist()]
    if len(vals) > max_items:
        return ",".join(vals[:max_items]) + f",...(+{len(vals)-max_items})"
    return ",".join(vals)


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
    if not IN_WITH_MARKERS.exists():
        raise SystemExit(f"ERROR: input file not found: {IN_WITH_MARKERS}")

    df = pd.read_csv(IN_WITH_MARKERS, sep="\t")

    required = {"chr", "pos", "cM"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: missing columns in {IN_WITH_MARKERS}: {missing}")

    df = df.copy()
    df["input_order"] = range(len(df))

    df["chr"] = normalize_chr_series(df["chr"])
    df["pos"] = df["pos"].astype(int)
    df["cM"] = df["cM"].astype(float)

    for col in ["marker_id", "marker_type", "source_mode", "evidence"]:
        if col not in df.columns:
            df[col] = ""

    # IMPORTANT:
    # For duplicated chr:pos, take the first row as it appears in
    # potato_genetic_map.hybrid_high_confidence.with_markers.tsv.
    first_rows = (
        df.sort_values("input_order")
        .groupby(["chr", "pos"], as_index=False, sort=True)
        .first()
    )

    final_map = first_rows[["chr", "pos", "cM"]].sort_values(["chr", "pos", "cM"])
    final_map.to_csv(OUT_MAP, sep="\t", index=False)

    duplicate_info = (
        df.groupby(["chr", "pos"])
        .agg(
            n_rows=("cM", "size"),
            n_unique_cM=("cM", "nunique"),
            first_cM=("cM", "first"),
            first_marker_id=("marker_id", "first"),
            first_marker_type=("marker_type", "first"),
            first_source_mode=("source_mode", "first"),
            cM_min=("cM", "min"),
            cM_median=("cM", "median"),
            cM_max=("cM", "max"),
            cM_range=("cM", lambda x: float(x.max() - x.min())),
            cM_values_in_input_order=("cM", join_values_in_order),
            marker_ids=("marker_id", join_unique),
            marker_types=("marker_type", join_unique),
            source_modes=("source_mode", join_unique),
        )
        .reset_index()
    )

    duplicate_problem = duplicate_info[
        (duplicate_info["n_rows"] > 1) & (duplicate_info["n_unique_cM"] > 1)
    ].copy()
    duplicate_problem.to_csv(OUT_DUPLICATES, sep="\t", index=False)

    metadata = first_rows.merge(
        duplicate_info,
        on=["chr", "pos"],
        how="left",
        suffixes=("", "_group"),
    )

    keep_cols = [
        "chr",
        "pos",
        "cM",
        "marker_id",
        "marker_type",
        "source_mode",
        "evidence",
        "n_rows",
        "n_unique_cM",
        "first_cM",
        "first_marker_id",
        "first_marker_type",
        "first_source_mode",
        "cM_min",
        "cM_median",
        "cM_max",
        "cM_range",
        "cM_values_in_input_order",
        "marker_ids",
        "marker_types",
        "source_modes",
    ]

    metadata[keep_cols].sort_values(["chr", "pos", "cM"]).to_csv(
        OUT_META, sep="\t", index=False
    )

    qc_by_chr = collinearity_qc(final_map)
    qc_by_chr.to_csv(OUT_QC_BY_CHR, sep="\t", index=False)

    flagged = qc_by_chr[
        (qc_by_chr["spearman_pos_cM"].abs() < 0.7)
        | (qc_by_chr["decreasing_fraction_without_ties"] > 0.35)
    ].copy()

    with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
        out.write("Potato first-cM hybrid high-confidence map\n")
        out.write("==========================================\n\n")
        out.write(f"Input with-marker rows: {len(df)}\n")
        out.write(f"Unique chr-pos rows after first-cM collapse: {len(final_map)}\n")
        out.write(f"Positions with multiple cM values before collapsing: {len(duplicate_problem)}\n\n")

        out.write("Collapse rule:\n")
        out.write("- For each unique chr:pos, cM was set to the first cM value in the input with-marker table.\n")
        out.write("- No median or averaging was applied.\n")
        out.write("- Full marker/cM information was preserved in the metadata table.\n\n")

        out.write("Collinearity QC by chromosome:\n")
        out.write(qc_by_chr.to_string(index=False))
        out.write("\n\n")

        out.write("Flagged chromosomes after first-cM collapse:\n")
        if flagged.empty:
            out.write("None by current thresholds.\n")
        else:
            out.write(flagged.to_string(index=False))
            out.write("\n")

        out.write("\nOutput files:\n")
        out.write(f"- first-cM map: {OUT_MAP}\n")
        out.write(f"- first-cM metadata: {OUT_META}\n")
        out.write(f"- duplicate positions before collapsing: {OUT_DUPLICATES}\n")
        out.write(f"- collinearity QC: {OUT_QC_BY_CHR}\n")

    print("Done.")
    print(f"First-cM map:       {OUT_MAP}")
    print(f"First-cM metadata:  {OUT_META}")
    print(f"Duplicate positions: {OUT_DUPLICATES}")
    print(f"QC by chr:          {OUT_QC_BY_CHR}")
    print(f"Summary:            {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

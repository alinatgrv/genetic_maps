#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


BASE = Path("species/solanum_tuberosum")
FINAL = BASE / "results/final"
QC = BASE / "results/qc"
QC.mkdir(parents=True, exist_ok=True)

MAPS = {
    "strict_direct_solcap": {
        "map": FINAL / "potato_genetic_map.tsv",
        "with_markers": FINAL / "potato_genetic_map.with_markers.tsv",
    },
    "liftover_extended_cigar": {
        "map": FINAL / "potato_genetic_map.liftover_extended_cigar.tsv",
        "with_markers": FINAL / "potato_genetic_map.liftover_extended_cigar.with_markers.tsv",
    },
    "hybrid_high_confidence": {
        "map": FINAL / "potato_genetic_map.hybrid_high_confidence.tsv",
        "with_markers": FINAL / "potato_genetic_map.hybrid_high_confidence.with_markers.tsv",
    },
}

OUT_BY_CHR = QC / "potato_collinearity_qc_by_map_chr.tsv"
OUT_SUMMARY = QC / "potato_collinearity_qc_summary.txt"
OUT_DECREASING = QC / "potato_collinearity_decreasing_steps.tsv"
OUT_DUPLICATES = QC / "potato_collinearity_duplicate_positions.tsv"


def spearman_without_scipy(x, y):
    x = pd.Series(x).rank(method="average")
    y = pd.Series(y).rank(method="average")
    if len(x) < 2:
        return np.nan
    return x.corr(y, method="pearson")


def normalize_chr_series(s):
    """Convert chr labels like 1, '1', 'chr01', 'chr1' to integer chromosome numbers."""
    return (
        s.astype(str)
        .str.strip()
        .str.replace("chr", "", regex=False)
        .str.replace("Chr", "", regex=False)
        .str.replace("CHR", "", regex=False)
        .astype(int)
    )


def load_with_optional_marker_info(map_path, marker_path):
    df = pd.read_csv(map_path, sep="\t")

    required = {"chr", "pos", "cM"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"ERROR: {map_path} missing columns: {missing}")

    df = df.copy()
    df["chr"] = normalize_chr_series(df["chr"])
    df["pos"] = df["pos"].astype(int)
    df["cM"] = df["cM"].astype(float)

    if marker_path.exists():
        m = pd.read_csv(marker_path, sep="\t")
        if {"chr", "pos", "cM"}.issubset(m.columns):
            m = m.copy()
            m["chr"] = normalize_chr_series(m["chr"])
            m["pos"] = m["pos"].astype(int)
            m["cM"] = m["cM"].astype(float)

            keep_cols = ["chr", "pos", "cM"]
            for col in ["marker_id", "marker_type", "source_mode", "evidence", "paf_identity", "mapq"]:
                if col in m.columns:
                    keep_cols.append(col)

            m = m[keep_cols].drop_duplicates()
            df = df.merge(m, on=["chr", "pos", "cM"], how="left")

    if "marker_id" not in df.columns:
        df["marker_id"] = ""
    if "marker_type" not in df.columns:
        df["marker_type"] = ""
    if "source_mode" not in df.columns:
        df["source_mode"] = ""

    return df


def analyze_map(map_name, df):
    rows = []
    decreasing_rows = []
    duplicate_rows = []

    for chr_id, g in df.groupby("chr", sort=True):
        # Some positions can have several cM values. Keep all rows for step analysis,
        # sorted by physical position and then cM.
        g = g.sort_values(["pos", "cM", "marker_id"]).reset_index(drop=True)

        n = len(g)
        pos_min = int(g["pos"].min())
        pos_max = int(g["pos"].max())
        cm_min = float(g["cM"].min())
        cm_max = float(g["cM"].max())

        rho = spearman_without_scipy(g["pos"], g["cM"]) if n >= 2 else np.nan

        dcM = g["cM"].diff()
        dpos = g["pos"].diff()

        valid_steps = (dpos > 0) & dcM.notna()

        inc = int(((dcM > 0) & valid_steps).sum())
        dec = int(((dcM < 0) & valid_steps).sum())
        eq = int(((dcM == 0) & valid_steps).sum())

        non_tie = inc + dec
        inc_fraction = inc / non_tie if non_tie else np.nan
        dec_fraction = dec / non_tie if non_tie else np.nan

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
                "map_name": map_name,
                "chr": chr_id,
                "n_rows": n,
                "n_unique_pos": g["pos"].nunique(),
                "n_unique_cM": g["cM"].nunique(),
                "pos_min": pos_min,
                "pos_max": pos_max,
                "cM_min": cm_min,
                "cM_max": cm_max,
                "spearman_pos_cM": rho,
                "orientation": orientation,
                "increasing_cM_steps": inc,
                "decreasing_cM_steps": dec,
                "equal_cM_steps": eq,
                "increasing_fraction_without_ties": inc_fraction,
                "decreasing_fraction_without_ties": dec_fraction,
            }
        )

        # Save local inversions/decreasing adjacent steps for inspection.
        for i in range(1, len(g)):
            if g.loc[i, "pos"] > g.loc[i - 1, "pos"] and g.loc[i, "cM"] < g.loc[i - 1, "cM"]:
                decreasing_rows.append(
                    {
                        "map_name": map_name,
                        "chr": chr_id,
                        "prev_pos": int(g.loc[i - 1, "pos"]),
                        "prev_cM": float(g.loc[i - 1, "cM"]),
                        "prev_marker_id": g.loc[i - 1, "marker_id"],
                        "prev_marker_type": g.loc[i - 1, "marker_type"],
                        "prev_source_mode": g.loc[i - 1, "source_mode"],
                        "curr_pos": int(g.loc[i, "pos"]),
                        "curr_cM": float(g.loc[i, "cM"]),
                        "curr_marker_id": g.loc[i, "marker_id"],
                        "curr_marker_type": g.loc[i, "marker_type"],
                        "curr_source_mode": g.loc[i, "source_mode"],
                        "delta_pos": int(g.loc[i, "pos"] - g.loc[i - 1, "pos"]),
                        "delta_cM": float(g.loc[i, "cM"] - g.loc[i - 1, "cM"]),
                    }
                )

        # Duplicate physical positions with more than one cM value.
        dup = (
            g.groupby("pos")
            .agg(
                n_rows=("pos", "size"),
                n_unique_cM=("cM", "nunique"),
                cM_values=("cM", lambda x: ",".join(map(str, sorted(set(x))))),
                marker_ids=("marker_id", lambda x: ",".join(map(str, sorted(set(x.dropna()))))),
                marker_types=("marker_type", lambda x: ",".join(map(str, sorted(set(x.dropna()))))),
            )
            .reset_index()
        )
        dup = dup[(dup["n_rows"] > 1) & (dup["n_unique_cM"] > 1)].copy()
        for _, r in dup.iterrows():
            duplicate_rows.append(
                {
                    "map_name": map_name,
                    "chr": chr_id,
                    "pos": int(r["pos"]),
                    "n_rows": int(r["n_rows"]),
                    "n_unique_cM": int(r["n_unique_cM"]),
                    "cM_values": r["cM_values"],
                    "marker_ids": r["marker_ids"],
                    "marker_types": r["marker_types"],
                }
            )

    return rows, decreasing_rows, duplicate_rows


def main():
    all_chr_rows = []
    all_decreasing = []
    all_duplicates = []

    for map_name, paths in MAPS.items():
        map_path = paths["map"]
        marker_path = paths["with_markers"]

        if not map_path.exists():
            print(f"WARNING: missing map file, skipping: {map_path}")
            continue

        df = load_with_optional_marker_info(map_path, marker_path)
        rows, decreasing, duplicates = analyze_map(map_name, df)

        all_chr_rows.extend(rows)
        all_decreasing.extend(decreasing)
        all_duplicates.extend(duplicates)

    by_chr = pd.DataFrame(all_chr_rows)
    by_chr.to_csv(OUT_BY_CHR, sep="\t", index=False)

    decreasing = pd.DataFrame(all_decreasing)
    if decreasing.empty:
        decreasing = pd.DataFrame(
            columns=[
                "map_name", "chr", "prev_pos", "prev_cM", "prev_marker_id",
                "curr_pos", "curr_cM", "curr_marker_id", "delta_pos", "delta_cM"
            ]
        )
    decreasing.to_csv(OUT_DECREASING, sep="\t", index=False)

    duplicates = pd.DataFrame(all_duplicates)
    if duplicates.empty:
        duplicates = pd.DataFrame(
            columns=["map_name", "chr", "pos", "n_rows", "n_unique_cM", "cM_values", "marker_ids", "marker_types"]
        )
    duplicates.to_csv(OUT_DUPLICATES, sep="\t", index=False)

    with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
        out.write("Potato genetic map collinearity QC\n")
        out.write("==================================\n\n")

        out.write("Input maps:\n")
        for map_name, paths in MAPS.items():
            out.write(f"- {map_name}: {paths['map']}\n")
        out.write("\n")

        out.write("Per-map overview:\n")
        overview = (
            by_chr.groupby("map_name")
            .agg(
                n_chromosomes=("chr", "nunique"),
                total_rows=("n_rows", "sum"),
                min_spearman=("spearman_pos_cM", "min"),
                median_spearman=("spearman_pos_cM", "median"),
                max_spearman=("spearman_pos_cM", "max"),
                total_decreasing_steps=("decreasing_cM_steps", "sum"),
                total_increasing_steps=("increasing_cM_steps", "sum"),
            )
            .reset_index()
        )
        out.write(overview.to_string(index=False))
        out.write("\n\n")

        out.write("Potentially problematic chromosomes:\n")
        flagged = by_chr[
            (by_chr["spearman_pos_cM"].abs() < 0.7)
            | (by_chr["decreasing_fraction_without_ties"] > 0.35)
        ].copy()

        if flagged.empty:
            out.write("None by current thresholds.\n")
        else:
            out.write(flagged.to_string(index=False))
            out.write("\n")

        out.write("\nDuplicate physical positions with multiple cM values:\n")
        if duplicates.empty:
            out.write("None.\n")
        else:
            dup_overview = (
                duplicates.groupby("map_name")
                .agg(n_duplicate_positions=("pos", "size"))
                .reset_index()
            )
            out.write(dup_overview.to_string(index=False))
            out.write("\n")

        out.write("\nOutput files:\n")
        out.write(f"- by chromosome: {OUT_BY_CHR}\n")
        out.write(f"- decreasing adjacent steps: {OUT_DECREASING}\n")
        out.write(f"- duplicate positions: {OUT_DUPLICATES}\n")

    print("Done.")
    print(f"By-chromosome QC:       {OUT_BY_CHR}")
    print(f"Decreasing steps:       {OUT_DECREASING}")
    print(f"Duplicate positions:    {OUT_DUPLICATES}")
    print(f"Summary:                {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

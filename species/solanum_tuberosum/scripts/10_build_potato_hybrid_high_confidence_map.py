#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


BASE = Path("species/solanum_tuberosum")

STRICT_DIRECT = BASE / "results/final/potato_genetic_map.with_markers.tsv"
CIGAR_LIFTED = BASE / "results/intermediate/liftover_v403_to_v61/potato_TableS4_liftover_v403_to_v61.cigar.all_markers.tsv"

OUT_FINAL = BASE / "results/final"
OUT_QC = BASE / "results/qc"

OUT_FINAL.mkdir(parents=True, exist_ok=True)
OUT_QC.mkdir(parents=True, exist_ok=True)

OUT_WITH_MARKERS = OUT_FINAL / "potato_genetic_map.hybrid_high_confidence.with_markers.tsv"
OUT_MAP = OUT_FINAL / "potato_genetic_map.hybrid_high_confidence.tsv"
OUT_QC_BY_CHR = OUT_QC / "potato_hybrid_high_confidence_qc_by_chr.tsv"
OUT_SUMMARY = OUT_QC / "potato_hybrid_high_confidence_summary.txt"


def normalize_marker_id(x):
    if pd.isna(x):
        return ""
    return str(x).strip().replace("solcap_stsnp_", "solcap_snp_").lower()


def find_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def load_strict_direct():
    df = pd.read_csv(STRICT_DIRECT, sep="\t")

    chr_col = find_col(df, ["chr", "target_chr", "physical_chr", "ref_chr", "chromosome"])
    pos_col = find_col(df, ["pos", "target_pos", "physical_pos", "SNP_POS", "snp_pos"])
    cm_col = find_col(df, ["cM", "cm", "genetic_position", "genetic_cm"])
    marker_col = find_col(df, ["marker_id", "snp_id", "marker", "Name", "name", "ID"])

    missing = []
    for label, col in [("chr", chr_col), ("pos", pos_col), ("cM", cm_col), ("marker_id", marker_col)]:
        if col is None:
            missing.append(label)

    if missing:
        raise SystemExit(
            f"ERROR: could not identify columns {missing} in {STRICT_DIRECT}. "
            f"Columns: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["chr"] = df[chr_col].astype(str).str.replace("chr", "", regex=False).astype(int)
    out["pos"] = df[pos_col].astype(int)
    out["cM"] = df[cm_col].astype(float)
    out["marker_id"] = df[marker_col].astype(str)
    out["marker_id_norm"] = out["marker_id"].map(normalize_marker_id)
    out["marker_type"] = "SNP"
    out["source_mode"] = "direct_spuddb_dm_v6_1"
    out["evidence"] = "SolCAP SNP exact non-multimap position from SpudDB DM v6.1"
    out["paf_identity"] = pd.NA
    out["mapq"] = pd.NA
    out["old_seqid"] = pd.NA
    out["old_midpoint"] = pd.NA

    return out


def load_filtered_liftover(direct_marker_ids):
    df = pd.read_csv(CIGAR_LIFTED, sep="\t")

    # High-confidence CIGAR liftover filters.
    filt = (
        (df["liftover_status"] == "lifted_chr_agree")
        & (df["cigar_status"] == "cigar_aligned_base")
        & (df["target_type"] == "chromosome")
        & (df["mapq"].astype(float) >= 60)
        & (df["paf_identity"].astype(float) >= 0.95)
        & (df["n_liftover_candidates"].astype(int) == 1)
    )

    filtered = df[filt].copy()

    # Do not use liftover positions for SolCAP SNPs if direct SpudDB positions are available.
    # Direct positions are more reliable and avoid known SolCAP liftover outliers.
    filtered["marker_id_norm"] = filtered["marker_id_norm"].map(normalize_marker_id)
    filtered = filtered[~filtered["marker_id_norm"].isin(direct_marker_ids)].copy()

    out = pd.DataFrame()
    out["chr"] = filtered["target_chr"].astype(int)
    out["pos"] = filtered["target_pos"].astype(int)
    out["cM"] = filtered["cM"].astype(float)
    out["marker_id"] = filtered["marker_id"].astype(str)
    out["marker_id_norm"] = filtered["marker_id_norm"].astype(str)
    out["marker_type"] = filtered["marker_type"].astype(str)
    out["source_mode"] = "cigar_liftover_v403_to_v61"
    out["evidence"] = "TableS4 marker lifted from PGSC/DM v4.03 to DM v6.1 using minimap2 PAF cg:Z CIGAR"
    out["paf_identity"] = filtered["paf_identity"].astype(float)
    out["mapq"] = filtered["mapq"].astype(int)
    out["old_seqid"] = filtered["old_seqid"].astype(str)
    out["old_midpoint"] = filtered["old_midpoint"].astype(int)

    return out, df, filtered


def main():
    if not STRICT_DIRECT.exists():
        raise SystemExit(f"ERROR: strict direct map not found: {STRICT_DIRECT}")

    if not CIGAR_LIFTED.exists():
        raise SystemExit(f"ERROR: CIGAR liftover file not found: {CIGAR_LIFTED}")

    direct = load_strict_direct()
    direct_marker_ids = set(direct["marker_id_norm"])

    liftover, all_lifted, filtered_liftover_raw = load_filtered_liftover(direct_marker_ids)

    combined = pd.concat([direct, liftover], ignore_index=True)
    combined = combined.sort_values(["chr", "pos", "cM", "marker_id", "source_mode"])

    combined.to_csv(OUT_WITH_MARKERS, sep="\t", index=False)

    final_map = combined[["chr", "pos", "cM"]].drop_duplicates()
    final_map = final_map.sort_values(["chr", "pos", "cM"])
    final_map.to_csv(OUT_MAP, sep="\t", index=False)

    qc_by_chr = (
        final_map.groupby("chr")
        .agg(
            n_markers=("pos", "size"),
            pos_min=("pos", "min"),
            pos_max=("pos", "max"),
            cM_min=("cM", "min"),
            cM_max=("cM", "max"),
        )
        .reset_index()
    )
    qc_by_chr.to_csv(OUT_QC_BY_CHR, sep="\t", index=False)

    with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
        out.write("Potato hybrid high-confidence genetic map\n")
        out.write("=========================================\n\n")
        out.write("Target reference: DM_1-3_516_R44_potato.v6.1\n\n")

        out.write("Inputs:\n")
        out.write(f"- Strict direct SolCAP map: {STRICT_DIRECT}\n")
        out.write(f"- CIGAR liftover TableS4 map: {CIGAR_LIFTED}\n\n")

        out.write("Construction rules:\n")
        out.write("- SolCAP SNP markers were taken from direct SpudDB DM v6.1 positions.\n")
        out.write("- Non-SolCAP markers were taken from CIGAR liftover only if:\n")
        out.write("  liftover_status == lifted_chr_agree\n")
        out.write("  cigar_status == cigar_aligned_base\n")
        out.write("  target_type == chromosome\n")
        out.write("  MAPQ >= 60\n")
        out.write("  paf_identity >= 0.95\n")
        out.write("  n_liftover_candidates == 1\n\n")

        out.write("Input counts:\n")
        out.write(f"Strict direct SolCAP rows: {len(direct)}\n")
        out.write(f"All CIGAR lifted TableS4 rows: {len(all_lifted)}\n")
        out.write(f"Filtered non-SolCAP CIGAR liftover rows added: {len(liftover)}\n")
        out.write(f"Combined with-marker rows: {len(combined)}\n")
        out.write(f"Final unique chr-pos-cM rows: {len(final_map)}\n\n")

        out.write("Rows by source_mode:\n")
        out.write(combined["source_mode"].value_counts().to_string())
        out.write("\n\n")

        out.write("Rows by marker_type:\n")
        out.write(combined["marker_type"].value_counts().to_string())
        out.write("\n\n")

        out.write("CIGAR liftover status counts before filtering:\n")
        out.write(all_lifted["liftover_status"].value_counts(dropna=False).to_string())
        out.write("\n\n")

        out.write("CIGAR marker type counts after high-confidence filtering before removing direct SNPs:\n")
        out.write(filtered_liftover_raw["marker_type"].value_counts(dropna=False).to_string())
        out.write("\n\n")

        out.write("QC by chromosome:\n")
        out.write(qc_by_chr.to_string(index=False))
        out.write("\n")

    print("Done.")
    print(f"Hybrid high-confidence map:      {OUT_MAP}")
    print(f"Hybrid high-confidence with IDs: {OUT_WITH_MARKERS}")
    print(f"QC by chr:                       {OUT_QC_BY_CHR}")
    print(f"Summary:                         {OUT_SUMMARY}")


if __name__ == "__main__":
    main()

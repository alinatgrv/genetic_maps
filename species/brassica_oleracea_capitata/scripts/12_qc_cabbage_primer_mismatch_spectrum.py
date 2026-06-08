from pathlib import Path
from collections import Counter, defaultdict

import pandas as pd


SPECIES_DIR = Path("species/brassica_oleracea_capitata")
ARCHIVE_CABBAGE_DIR = Path("archive_old_structure/results/cabbage")

METADATA_FILE = SPECIES_DIR / "data/metadata/cabbage_markers_metadata.tsv"
BLAST_FILE = ARCHIVE_CABBAGE_DIR / "cabbage_primers_vs_ref.blast.tsv"

QC_DIR = SPECIES_DIR / "results/qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

PRIMER_SPECTRUM_OUT = QC_DIR / "cabbage_primer_mismatch_spectrum.tsv"
PRIMER_DISTRIBUTION_OUT = QC_DIR / "cabbage_primer_min_mismatch_distribution.tsv"
MARKER_SPECTRUM_OUT = QC_DIR / "cabbage_marker_min_mismatch_spectrum.tsv"
MARKER_THRESHOLD_OUT = QC_DIR / "cabbage_marker_rescue_potential_by_threshold.tsv"
MARKER_THRESHOLD_TYPE_OUT = QC_DIR / "cabbage_marker_rescue_potential_by_threshold_and_type.tsv"
RESCUE_CANDIDATES_OUT = QC_DIR / "cabbage_marker_rescue_candidates_mm2_mm3.tsv"
SUMMARY_OUT = QC_DIR / "cabbage_primer_mismatch_spectrum.summary.md"


MM_EXACT_MAX = 5
THRESHOLDS = [0, 1, 2, 3, 4, 5]


def require_file(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Не найден файл: {path}")


def get_col(df, candidates, required=True):
    for c in candidates:
        if c in df.columns:
            return c
    if required:
        raise ValueError(f"Не найдена колонка из вариантов: {candidates}")
    return None


def empty_stat():
    stat = {
        "raw_hits": 0,
        "best_any_alignment_length": 0,
        "best_any_coverage": 0.0,
        "best_any_pident": None,
        "best_any_mismatch": None,
        "best_any_gapopen": None,
        "best_any_sseqid": None,
        "best_any_sstart": None,
        "best_any_send": None,
        "best_any_bitscore": None,

        "full_len_no_gap_hits": 0,
        "min_mismatch_full_len_no_gap": None,
        "hits_at_min_mismatch": 0,
        "best_full_len_no_gap_sseqid": None,
        "best_full_len_no_gap_sstart": None,
        "best_full_len_no_gap_send": None,
        "best_full_len_no_gap_strand": None,
        "best_full_len_no_gap_bitscore": None,
    }

    for mm in range(MM_EXACT_MAX + 1):
        stat[f"full_len_no_gap_mm{mm}_hits"] = 0

    stat["full_len_no_gap_mm_gt5_hits"] = 0

    for t in THRESHOLDS:
        stat[f"full_len_no_gap_mm_le{t}_hits"] = 0

    return stat


def parse_primer_id(qseqid):
    if "__" not in qseqid:
        return qseqid, None
    marker_id, primer_side = qseqid.rsplit("__", 1)
    return marker_id, primer_side


def better_any_hit(stat, aln_len, coverage, pident, bitscore):
    """
    Для объяснения проблемных праймеров:
    лучший любой BLAST-hit выбираем по coverage, затем pident, затем bitscore.
    """
    if coverage > stat["best_any_coverage"]:
        return True
    if coverage == stat["best_any_coverage"]:
        old_pident = stat["best_any_pident"]
        if old_pident is None or pident > old_pident:
            return True
        if pident == old_pident:
            old_bitscore = stat["best_any_bitscore"]
            if old_bitscore is None or bitscore > old_bitscore:
                return True
    return False


def better_full_len_hit(stat, mismatch, bitscore):
    """
    Для full-length no-gap hit:
    сначала меньше mismatch, затем выше bitscore.
    """
    old_mm = stat["min_mismatch_full_len_no_gap"]
    if old_mm is None:
        return True
    if mismatch < old_mm:
        return True
    if mismatch == old_mm:
        old_bitscore = stat["best_full_len_no_gap_bitscore"]
        if old_bitscore is None or bitscore > old_bitscore:
            return True
    return False


def scan_blast(expected_primers):
    stats = {q: empty_stat() for q in expected_primers}

    with open(BLAST_FILE, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                continue

            qseqid = fields[0]
            sseqid = fields[1]

            try:
                pident = float(fields[2])
                aln_len = int(float(fields[3]))
                mismatch = int(float(fields[4]))
                gapopen = int(float(fields[5]))
                sstart = int(float(fields[8]))
                send = int(float(fields[9]))
                bitscore = float(fields[11])
                qlen = int(float(fields[12]))
            except ValueError:
                continue

            if qseqid not in stats:
                stats[qseqid] = empty_stat()

            stat = stats[qseqid]
            stat["raw_hits"] += 1

            coverage = aln_len / qlen if qlen else 0

            if better_any_hit(stat, aln_len, coverage, pident, bitscore):
                stat["best_any_alignment_length"] = aln_len
                stat["best_any_coverage"] = coverage
                stat["best_any_pident"] = pident
                stat["best_any_mismatch"] = mismatch
                stat["best_any_gapopen"] = gapopen
                stat["best_any_sseqid"] = sseqid
                stat["best_any_sstart"] = sstart
                stat["best_any_send"] = send
                stat["best_any_bitscore"] = bitscore

            full_len_no_gap = (aln_len == qlen) and (gapopen == 0)

            if full_len_no_gap:
                stat["full_len_no_gap_hits"] += 1

                if mismatch <= MM_EXACT_MAX:
                    stat[f"full_len_no_gap_mm{mismatch}_hits"] += 1
                else:
                    stat["full_len_no_gap_mm_gt5_hits"] += 1

                for t in THRESHOLDS:
                    if mismatch <= t:
                        stat[f"full_len_no_gap_mm_le{t}_hits"] += 1

                old_min = stat["min_mismatch_full_len_no_gap"]
                if old_min is None or mismatch < old_min:
                    stat["min_mismatch_full_len_no_gap"] = mismatch
                    stat["hits_at_min_mismatch"] = 1
                elif mismatch == old_min:
                    stat["hits_at_min_mismatch"] += 1

                if better_full_len_hit(stat, mismatch, bitscore):
                    stat["best_full_len_no_gap_sseqid"] = sseqid
                    stat["best_full_len_no_gap_sstart"] = sstart
                    stat["best_full_len_no_gap_send"] = send
                    stat["best_full_len_no_gap_strand"] = "+" if sstart <= send else "-"
                    stat["best_full_len_no_gap_bitscore"] = bitscore

    return stats


def min_mm_bin(x):
    if pd.isna(x):
        return "no_full_length_no_gap"
    x = int(x)
    if x <= 5:
        return str(x)
    return ">5"


def main():
    require_file(METADATA_FILE)
    require_file(BLAST_FILE)

    metadata = pd.read_csv(METADATA_FILE, sep="\t")

    marker_col = get_col(metadata, ["marker_id", "Marker Name", "marker_name"])
    type_col = get_col(metadata, ["marker_type", "Marker Type", "type"], required=False)
    lg_col = get_col(metadata, ["linkage_group", "LG"], required=False)
    cm_col = get_col(metadata, ["cM", "Position", "position"], required=False)

    metadata = metadata.copy()
    metadata["marker_id_for_qc"] = metadata[marker_col].astype(str)

    expected_primers = set()
    for marker_id in metadata["marker_id_for_qc"]:
        expected_primers.add(f"{marker_id}__F")
        expected_primers.add(f"{marker_id}__R")

    print("Scanning BLAST file for primer mismatch spectrum...")
    stats = scan_blast(expected_primers)

    primer_rows = []

    meta_by_marker = metadata.set_index("marker_id_for_qc", drop=False)

    for qseqid in sorted(expected_primers):
        marker_id, primer_side = parse_primer_id(qseqid)
        stat = stats.get(qseqid, empty_stat())

        marker_type = meta_by_marker.loc[marker_id, type_col] if type_col and marker_id in meta_by_marker.index else None
        linkage_group = meta_by_marker.loc[marker_id, lg_col] if lg_col and marker_id in meta_by_marker.index else None
        cM = meta_by_marker.loc[marker_id, cm_col] if cm_col and marker_id in meta_by_marker.index else None

        row = {
            "marker_id": marker_id,
            "primer_side": primer_side,
            "qseqid": qseqid,
            "marker_type": marker_type,
            "linkage_group": linkage_group,
            "cM": cM,

            "raw_hits": stat["raw_hits"],
            "full_len_no_gap_hits": stat["full_len_no_gap_hits"],
            "min_mismatch_full_len_no_gap": stat["min_mismatch_full_len_no_gap"],
            "min_mismatch_bin": min_mm_bin(stat["min_mismatch_full_len_no_gap"]),
            "hits_at_min_mismatch": stat["hits_at_min_mismatch"],
            "is_multihit_at_min_mismatch": stat["hits_at_min_mismatch"] > 1,

            "best_full_len_no_gap_sseqid": stat["best_full_len_no_gap_sseqid"],
            "best_full_len_no_gap_sstart": stat["best_full_len_no_gap_sstart"],
            "best_full_len_no_gap_send": stat["best_full_len_no_gap_send"],
            "best_full_len_no_gap_strand": stat["best_full_len_no_gap_strand"],
            "best_full_len_no_gap_bitscore": stat["best_full_len_no_gap_bitscore"],

            "best_any_alignment_length": stat["best_any_alignment_length"],
            "best_any_coverage": stat["best_any_coverage"],
            "best_any_pident": stat["best_any_pident"],
            "best_any_mismatch": stat["best_any_mismatch"],
            "best_any_gapopen": stat["best_any_gapopen"],
            "best_any_sseqid": stat["best_any_sseqid"],
            "best_any_sstart": stat["best_any_sstart"],
            "best_any_send": stat["best_any_send"],
            "best_any_bitscore": stat["best_any_bitscore"],
        }

        for mm in range(MM_EXACT_MAX + 1):
            row[f"full_len_no_gap_mm{mm}_hits"] = stat[f"full_len_no_gap_mm{mm}_hits"]

        row["full_len_no_gap_mm_gt5_hits"] = stat["full_len_no_gap_mm_gt5_hits"]

        for t in THRESHOLDS:
            row[f"full_len_no_gap_mm_le{t}_hits"] = stat[f"full_len_no_gap_mm_le{t}_hits"]
            row[f"has_full_len_no_gap_mm_le{t}"] = stat[f"full_len_no_gap_mm_le{t}_hits"] > 0

        primer_rows.append(row)

    primer_df = pd.DataFrame(primer_rows)
    primer_df.to_csv(PRIMER_SPECTRUM_OUT, sep="\t", index=False)

    # Distribution by primer minimum mismatch
    dist_order = ["0", "1", "2", "3", "4", "5", ">5", "no_full_length_no_gap"]
    dist = (
        primer_df["min_mismatch_bin"]
        .value_counts()
        .reindex(dist_order, fill_value=0)
        .reset_index()
    )
    dist.columns = ["min_mismatch_bin", "n_primers"]
    dist["percent_primers"] = dist["n_primers"] / len(primer_df) * 100
    dist.to_csv(PRIMER_DISTRIBUTION_OUT, sep="\t", index=False)

    # Marker-level table
    marker_rows = []

    for marker_id, g in primer_df.groupby("marker_id"):
        f = g[g["primer_side"] == "F"]
        r = g[g["primer_side"] == "R"]

        if len(f) == 0 or len(r) == 0:
            continue

        f = f.iloc[0]
        r = r.iloc[0]

        f_min = f["min_mismatch_full_len_no_gap"]
        r_min = r["min_mismatch_full_len_no_gap"]

        f_has_full = not pd.isna(f_min)
        r_has_full = not pd.isna(r_min)
        both_have_full = f_has_full and r_has_full

        if both_have_full:
            max_min_mm = max(int(f_min), int(r_min))
            sum_min_mm = int(f_min) + int(r_min)
        else:
            max_min_mm = None
            sum_min_mm = None

        row = {
            "marker_id": marker_id,
            "marker_type": f["marker_type"],
            "linkage_group": f["linkage_group"],
            "cM": f["cM"],

            "F_min_mismatch": f_min,
            "R_min_mismatch": r_min,
            "F_min_mismatch_bin": f["min_mismatch_bin"],
            "R_min_mismatch_bin": r["min_mismatch_bin"],
            "both_primers_have_full_len_no_gap_hit": both_have_full,

            "max_primer_min_mismatch": max_min_mm,
            "sum_primer_min_mismatch": sum_min_mm,

            "F_hits_at_min_mismatch": f["hits_at_min_mismatch"],
            "R_hits_at_min_mismatch": r["hits_at_min_mismatch"],
            "F_multihit_at_min_mismatch": f["is_multihit_at_min_mismatch"],
            "R_multihit_at_min_mismatch": r["is_multihit_at_min_mismatch"],

            "F_best_any_coverage": f["best_any_coverage"],
            "R_best_any_coverage": r["best_any_coverage"],
            "F_best_any_gapopen": f["best_any_gapopen"],
            "R_best_any_gapopen": r["best_any_gapopen"],
        }

        for t in THRESHOLDS:
            row[f"both_primers_min_mm_le{t}"] = (
                both_have_full and int(f_min) <= t and int(r_min) <= t
            )
            row[f"F_has_mm_le{t}"] = bool(f[f"has_full_len_no_gap_mm_le{t}"])
            row[f"R_has_mm_le{t}"] = bool(r[f"has_full_len_no_gap_mm_le{t}"])

        marker_rows.append(row)

    marker_df = pd.DataFrame(marker_rows)
    marker_df.to_csv(MARKER_SPECTRUM_OUT, sep="\t", index=False)

    # Threshold summary
    threshold_rows = []

    previous_n = 0
    previous_ids = set()

    for t in THRESHOLDS:
        col = f"both_primers_min_mm_le{t}"
        eligible = marker_df[marker_df[col]].copy()
        eligible_ids = set(eligible["marker_id"])

        threshold_rows.append({
            "threshold_per_primer": t,
            "n_markers_with_both_primers_min_mm_le_threshold": len(eligible),
            "new_markers_vs_previous_threshold": len(eligible_ids - previous_ids),
            "n_markers_with_any_multihit_at_min_mismatch": int(
                (eligible["F_multihit_at_min_mismatch"] | eligible["R_multihit_at_min_mismatch"]).sum()
            ),
            "note": f"Оба праймера имеют full-length no-gap hit с min mismatch <= {t}",
        })

        previous_n = len(eligible)
        previous_ids = eligible_ids

    threshold_df = pd.DataFrame(threshold_rows)
    threshold_df.to_csv(MARKER_THRESHOLD_OUT, sep="\t", index=False)

    # Threshold by marker type
    type_rows = []
    if "marker_type" in marker_df.columns:
        for t in THRESHOLDS:
            col = f"both_primers_min_mm_le{t}"
            eligible = marker_df[marker_df[col]].copy()
            counts = eligible["marker_type"].fillna("NA").value_counts()
            for marker_type, n in counts.items():
                type_rows.append({
                    "threshold_per_primer": t,
                    "marker_type": marker_type,
                    "n_markers": int(n),
                })

    pd.DataFrame(type_rows).to_csv(MARKER_THRESHOLD_TYPE_OUT, sep="\t", index=False)

    # Rescue candidates for possible mm2/mm3 map
    rescue = marker_df[
        marker_df["both_primers_have_full_len_no_gap_hit"]
        & marker_df["max_primer_min_mismatch"].isin([2, 3])
    ].copy()

    rescue["has_multihit_at_min_mismatch"] = (
        rescue["F_multihit_at_min_mismatch"] | rescue["R_multihit_at_min_mismatch"]
    )

    rescue = rescue.sort_values([
        "max_primer_min_mismatch",
        "sum_primer_min_mismatch",
        "has_multihit_at_min_mismatch",
        "marker_type",
        "linkage_group",
        "cM",
        "marker_id",
    ])

    rescue.to_csv(RESCUE_CANDIDATES_OUT, sep="\t", index=False)

    # Summary markdown
    with open(SUMMARY_OUT, "w", encoding="utf-8") as out:
        out.write("# Cabbage primer mismatch spectrum QC\n\n")

        out.write("## Primer-level minimum mismatch distribution\n\n")
        out.write("```tsv\n" + dist.to_csv(sep="\\t", index=False) + "```")
        out.write("\n\n")

        out.write("## Marker rescue potential by threshold\n\n")
        out.write("```tsv\n" + threshold_df.to_csv(sep="\\t", index=False) + "```")
        out.write("\n\n")

        out.write("## Output files\n\n")
        out.write(f"- Primer spectrum: `{PRIMER_SPECTRUM_OUT}`\n")
        out.write(f"- Primer min mismatch distribution: `{PRIMER_DISTRIBUTION_OUT}`\n")
        out.write(f"- Marker spectrum: `{MARKER_SPECTRUM_OUT}`\n")
        out.write(f"- Marker rescue potential by threshold: `{MARKER_THRESHOLD_OUT}`\n")
        out.write(f"- Marker rescue potential by type: `{MARKER_THRESHOLD_TYPE_OUT}`\n")
        out.write(f"- Candidate markers for mm2/mm3 exploration: `{RESCUE_CANDIDATES_OUT}`\n")

    print("\nDone.")
    print(f"Primer spectrum: {PRIMER_SPECTRUM_OUT}")
    print(f"Primer distribution: {PRIMER_DISTRIBUTION_OUT}")
    print(f"Marker spectrum: {MARKER_SPECTRUM_OUT}")
    print(f"Marker rescue potential: {MARKER_THRESHOLD_OUT}")
    print(f"Marker rescue by type: {MARKER_THRESHOLD_TYPE_OUT}")
    print(f"mm2/mm3 rescue candidates: {RESCUE_CANDIDATES_OUT}")
    print(f"Summary: {SUMMARY_OUT}")

    print("\nPrimer minimum mismatch distribution:")
    print(dist.to_string(index=False))

    print("\nMarker rescue potential:")
    print(threshold_df.to_string(index=False))

    print("\nMarker rescue potential by marker type:")
    if type_rows:
        print(pd.DataFrame(type_rows).to_string(index=False))
    else:
        print("No marker type information available.")

    print("\nTop 30 mm2/mm3 rescue candidates:")
    if len(rescue):
        cols = [
            "marker_id", "marker_type", "linkage_group", "cM",
            "F_min_mismatch", "R_min_mismatch",
            "max_primer_min_mismatch", "sum_primer_min_mismatch",
            "F_hits_at_min_mismatch", "R_hits_at_min_mismatch",
            "has_multihit_at_min_mismatch",
        ]
        print(rescue[cols].head(30).to_string(index=False))
    else:
        print("No mm2/mm3 rescue candidates found.")


if __name__ == "__main__":
    main()

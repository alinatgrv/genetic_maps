import csv
from pathlib import Path
from collections import defaultdict

import pandas as pd


SPECIES_DIR = Path("species/brassica_oleracea_capitata")
ARCHIVE_CABBAGE_DIR = Path("archive_old_structure/results/cabbage")

METADATA_FILE = SPECIES_DIR / "data/metadata/cabbage_markers_metadata.tsv"

BLAST_FILE = ARCHIVE_CABBAGE_DIR / "cabbage_primers_vs_ref.blast.tsv"

FINAL_STRICT_FILE = SPECIES_DIR / "results/final/cabbage_genetic_map.high_confidence_196.tsv"
FINAL_RELAXED_FILE = SPECIES_DIR / "results/final/cabbage_genetic_map.relaxed_mm1_collinear_222.tsv"

QC_DIR = SPECIES_DIR / "results/qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

PRIMER_QC_OUT = QC_DIR / "cabbage_primer_alignment_summary.tsv"
MARKER_QC_OUT = QC_DIR / "cabbage_marker_alignment_summary.tsv"
FUNNEL_OUT = QC_DIR / "cabbage_marker_filtering_funnel.tsv"
GAP_STATS_OUT = QC_DIR / "cabbage_final_map_gap_statistics.by_chromosome.tsv"
LARGE_GAPS_OUT = QC_DIR / "cabbage_final_map_large_gaps.tsv"
SUMMARY_MD_OUT = QC_DIR / "cabbage_qc_summary.md"


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(
            f"Не найден файл: {path}\n"
            f"Проверь путь или восстанови файл из старой рабочей структуры."
        )


def first_existing(paths):
    for path in paths:
        path = Path(path)
        if path.exists():
            return path
    return None


def count_rows_tsv(path: Path) -> int | None:
    if path is None or not path.exists():
        return None
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        n = sum(1 for _ in handle)
    return max(n - 1, 0)


def load_tsv_if_exists(path: Path):
    if path is None or not path.exists():
        return None
    return pd.read_csv(path, sep="\t")


def safe_nunique_marker_id(path: Path):
    df = load_tsv_if_exists(path)
    if df is None:
        return None
    if "marker_id" not in df.columns:
        return None
    return int(df["marker_id"].nunique())


def get_metadata_marker_col(df: pd.DataFrame) -> str:
    if "marker_id" in df.columns:
        return "marker_id"
    for candidate in ["Marker Name", "marker_name", "name"]:
        if candidate in df.columns:
            return candidate
    raise ValueError(
        "Не нашла колонку с ID маркера. Ожидаю marker_id или Marker Name."
    )


def get_metadata_type_col(df: pd.DataFrame) -> str | None:
    for candidate in ["marker_type", "Marker Type", "type"]:
        if candidate in df.columns:
            return candidate
    return None


def get_metadata_lg_col(df: pd.DataFrame) -> str | None:
    for candidate in ["linkage_group", "LG"]:
        if candidate in df.columns:
            return candidate
    return None


def get_metadata_cm_col(df: pd.DataFrame) -> str | None:
    for candidate in ["cM", "Position", "position"]:
        if candidate in df.columns:
            return candidate
    return None


def parse_primer_qseqid(qseqid: str):
    if "__" not in qseqid:
        return None, None
    marker_id, side = qseqid.rsplit("__", 1)
    if side not in {"F", "R"}:
        return marker_id, None
    return marker_id, side


def scan_blast_file(blast_file: Path, expected_primers: set[str]) -> dict:
    """
    Streaming scan of BLAST outfmt 6:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

    Counts:
    - raw_hits: any BLAST hit
    - full_length_hits: alignment length == qlen and no gaps
    - full_length_perfect_hits: full-length, no gaps, mismatch == 0
    - full_length_mm1_hits: full-length, no gaps, mismatch <= 1
    """
    stats = {
        q: {
            "raw_hits": 0,
            "full_length_hits": 0,
            "full_length_perfect_hits": 0,
            "full_length_mm1_hits": 0,
            "min_mismatch_full_length": None,
            "best_bitscore": None,
        }
        for q in expected_primers
    }

    extra_queries = set()

    with open(blast_file, "r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                continue

            qseqid = fields[0]

            if qseqid not in stats:
                extra_queries.add(qseqid)
                stats[qseqid] = {
                    "raw_hits": 0,
                    "full_length_hits": 0,
                    "full_length_perfect_hits": 0,
                    "full_length_mm1_hits": 0,
                    "min_mismatch_full_length": None,
                    "best_bitscore": None,
                }

            try:
                aln_len = int(float(fields[3]))
                mismatch = int(float(fields[4]))
                gapopen = int(float(fields[5]))
                bitscore = float(fields[11])
                qlen = int(float(fields[12]))
            except ValueError:
                continue

            stats[qseqid]["raw_hits"] += 1

            if stats[qseqid]["best_bitscore"] is None or bitscore > stats[qseqid]["best_bitscore"]:
                stats[qseqid]["best_bitscore"] = bitscore

            full_length_no_gap = (aln_len == qlen) and (gapopen == 0)

            if full_length_no_gap:
                stats[qseqid]["full_length_hits"] += 1

                current_min = stats[qseqid]["min_mismatch_full_length"]
                if current_min is None or mismatch < current_min:
                    stats[qseqid]["min_mismatch_full_length"] = mismatch

                if mismatch == 0:
                    stats[qseqid]["full_length_perfect_hits"] += 1

                if mismatch <= 1:
                    stats[qseqid]["full_length_mm1_hits"] += 1

    if extra_queries:
        print(
            f"Warning: found {len(extra_queries)} query IDs in BLAST that were not expected from metadata."
        )

    return stats


def make_primer_qc(metadata: pd.DataFrame, blast_stats: dict, marker_col: str) -> pd.DataFrame:
    rows = []

    type_col = get_metadata_type_col(metadata)
    lg_col = get_metadata_lg_col(metadata)
    cm_col = get_metadata_cm_col(metadata)

    for _, row in metadata.iterrows():
        marker_id = str(row[marker_col])

        marker_type = row[type_col] if type_col else None
        linkage_group = row[lg_col] if lg_col else None
        cM = row[cm_col] if cm_col else None

        for side in ["F", "R"]:
            qseqid = f"{marker_id}__{side}"
            s = blast_stats.get(qseqid, {
                "raw_hits": 0,
                "full_length_hits": 0,
                "full_length_perfect_hits": 0,
                "full_length_mm1_hits": 0,
                "min_mismatch_full_length": None,
                "best_bitscore": None,
            })

            rows.append({
                "marker_id": marker_id,
                "primer_side": side,
                "qseqid": qseqid,
                "linkage_group": linkage_group,
                "cM": cM,
                "marker_type": marker_type,
                "raw_hits": s["raw_hits"],
                "full_length_hits": s["full_length_hits"],
                "full_length_perfect_hits": s["full_length_perfect_hits"],
                "full_length_mm1_hits": s["full_length_mm1_hits"],
                "has_any_blast_hit": s["raw_hits"] > 0,
                "has_full_length_hit": s["full_length_hits"] > 0,
                "has_full_length_perfect_hit": s["full_length_perfect_hits"] > 0,
                "has_full_length_mm1_hit": s["full_length_mm1_hits"] > 0,
                "is_multihit_perfect": s["full_length_perfect_hits"] > 1,
                "is_multihit_mm1": s["full_length_mm1_hits"] > 1,
                "min_mismatch_full_length": s["min_mismatch_full_length"],
                "best_bitscore": s["best_bitscore"],
            })

    return pd.DataFrame(rows)


def make_marker_qc(primer_qc: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for marker_id, g in primer_qc.groupby("marker_id"):
        f = g[g["primer_side"] == "F"]
        r = g[g["primer_side"] == "R"]

        f_row = f.iloc[0] if len(f) else None
        r_row = r.iloc[0] if len(r) else None

        def get(row, col, default=0):
            if row is None:
                return default
            return row[col]

        rows.append({
            "marker_id": marker_id,
            "linkage_group": get(f_row, "linkage_group", None),
            "cM": get(f_row, "cM", None),
            "marker_type": get(f_row, "marker_type", None),

            "F_raw_hits": get(f_row, "raw_hits"),
            "R_raw_hits": get(r_row, "raw_hits"),
            "F_full_length_perfect_hits": get(f_row, "full_length_perfect_hits"),
            "R_full_length_perfect_hits": get(r_row, "full_length_perfect_hits"),
            "F_full_length_mm1_hits": get(f_row, "full_length_mm1_hits"),
            "R_full_length_mm1_hits": get(r_row, "full_length_mm1_hits"),

            "has_F_any_hit": bool(get(f_row, "has_any_blast_hit", False)),
            "has_R_any_hit": bool(get(r_row, "has_any_blast_hit", False)),
            "has_both_any_hit": bool(get(f_row, "has_any_blast_hit", False)) and bool(get(r_row, "has_any_blast_hit", False)),

            "has_F_perfect": bool(get(f_row, "has_full_length_perfect_hit", False)),
            "has_R_perfect": bool(get(r_row, "has_full_length_perfect_hit", False)),
            "has_both_primers_perfect": bool(get(f_row, "has_full_length_perfect_hit", False)) and bool(get(r_row, "has_full_length_perfect_hit", False)),

            "has_F_mm1": bool(get(f_row, "has_full_length_mm1_hit", False)),
            "has_R_mm1": bool(get(r_row, "has_full_length_mm1_hit", False)),
            "has_both_primers_mm1": bool(get(f_row, "has_full_length_mm1_hit", False)) and bool(get(r_row, "has_full_length_mm1_hit", False)),

            "F_multihit_perfect": bool(get(f_row, "is_multihit_perfect", False)),
            "R_multihit_perfect": bool(get(r_row, "is_multihit_perfect", False)),
            "F_multihit_mm1": bool(get(f_row, "is_multihit_mm1", False)),
            "R_multihit_mm1": bool(get(r_row, "is_multihit_mm1", False)),
        })

    return pd.DataFrame(rows)


def add_final_membership(marker_qc: pd.DataFrame) -> pd.DataFrame:
    strict_with_markers = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.strict.clean.collinear_with_markers.tsv",
        SPECIES_DIR / "results/qc/cabbage_genetic_map.strict.clean.collinear_with_markers.tsv",
    ])

    relaxed_with_markers = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.relaxed_primer_mm1.clean.collinear_with_markers.tsv",
        SPECIES_DIR / "results/qc/cabbage_genetic_map.relaxed_primer_mm1.clean.collinear_with_markers.tsv",
    ])

    marker_qc["in_strict_high_confidence_196"] = False
    marker_qc["in_relaxed_mm1_collinear_222"] = False

    if strict_with_markers:
        strict_df = pd.read_csv(strict_with_markers, sep="\t")
        if "marker_id" in strict_df.columns:
            strict_ids = set(strict_df["marker_id"].astype(str))
            marker_qc["in_strict_high_confidence_196"] = marker_qc["marker_id"].astype(str).isin(strict_ids)

    if relaxed_with_markers:
        relaxed_df = pd.read_csv(relaxed_with_markers, sep="\t")
        if "marker_id" in relaxed_df.columns:
            relaxed_ids = set(relaxed_df["marker_id"].astype(str))
            marker_qc["in_relaxed_mm1_collinear_222"] = marker_qc["marker_id"].astype(str).isin(relaxed_ids)

    return marker_qc


def make_funnel(metadata: pd.DataFrame, primer_qc: pd.DataFrame, marker_qc: pd.DataFrame) -> pd.DataFrame:
    strict_amplicons_file = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_amplicons.tsv",
        ARCHIVE_CABBAGE_DIR / "cabbage_amplicons.strict.tsv",
        ARCHIVE_CABBAGE_DIR / "cabbage_amplicons.perfect.tsv",
        ARCHIVE_CABBAGE_DIR / "cabbage_amplicons.strict_perfect.tsv",
    ])

    strict_raw_map_file = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.strict.tsv",
    ])

    strict_clean_unique_file = first_existing([
        SPECIES_DIR / "results/intermediate/cabbage_genetic_map.strict_296.tsv",
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.strict.clean.unique.tsv",
    ])

    relaxed_amplicons_file = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_amplicons.relaxed_primer_mm1.tsv",
    ])

    relaxed_raw_map_file = first_existing([
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.relaxed_primer_mm1.tsv",
    ])

    relaxed_clean_unique_file = first_existing([
        SPECIES_DIR / "results/intermediate/cabbage_genetic_map.relaxed_primer_mm1.clean.unique.tsv",
        ARCHIVE_CABBAGE_DIR / "cabbage_genetic_map.relaxed_primer_mm1.clean.unique.tsv",
    ])

    strict_amplicons = load_tsv_if_exists(strict_amplicons_file)
    relaxed_amplicons = load_tsv_if_exists(relaxed_amplicons_file)

    rows = []

    def add(strategy, stage_order, stage, n, note):
        rows.append({
            "strategy": strategy,
            "stage_order": stage_order,
            "stage": stage,
            "n": n,
            "note": note,
        })

    total_markers = len(metadata)
    total_primers = len(primer_qc)

    add("input", 1, "source_markers", total_markers, "Исходные маркеры из supplementary-таблицы")
    add("input", 2, "source_primers", total_primers, "Forward + reverse primers")

    add("primer_alignment", 3, "primers_with_any_blast_hit", int(primer_qc["has_any_blast_hit"].sum()), "Праймер имеет хотя бы одно BLAST-попадание")
    add("primer_alignment", 4, "primers_with_full_length_perfect_hit", int(primer_qc["has_full_length_perfect_hit"].sum()), "Полное идеальное попадание: full-length, mismatch=0, gapopen=0")
    add("primer_alignment", 5, "primers_with_full_length_mm1_hit", int(primer_qc["has_full_length_mm1_hit"].sum()), "Полное попадание с mismatch <= 1, gapopen=0")

    add("strict", 10, "markers_with_both_primers_perfect", int(marker_qc["has_both_primers_perfect"].sum()), "Оба праймера имеют full-length perfect hit")

    if strict_amplicons is not None and "marker_id" in strict_amplicons.columns:
        add("strict", 11, "candidate_amplicons", len(strict_amplicons), f"Файл: {strict_amplicons_file}")
        add("strict", 12, "markers_with_candidate_amplicon", int(strict_amplicons["marker_id"].nunique()), f"Файл: {strict_amplicons_file}")
        candidates_per_marker = strict_amplicons.groupby("marker_id").size()
        add("strict", 13, "markers_with_exactly_one_candidate_amplicon", int((candidates_per_marker == 1).sum()), "Однозначный strict-ампликон")
    else:
        add("strict", 11, "candidate_amplicons", None, "Файл candidate amplicons не найден")
        add("strict", 12, "markers_with_candidate_amplicon", safe_nunique_marker_id(strict_raw_map_file), f"Оценка по файлу карты: {strict_raw_map_file}")
        add("strict", 13, "markers_with_exactly_one_candidate_amplicon", count_rows_tsv(strict_raw_map_file), f"Оценка по файлу карты: {strict_raw_map_file}")

    add("strict", 14, "strict_clean_unique", count_rows_tsv(strict_clean_unique_file), "После удаления scaffolds, LG/ref conflicts и duplicate positions")
    add("strict", 15, "strict_high_confidence_collinear", count_rows_tsv(FINAL_STRICT_FILE), "Финальная strict high-confidence карта")

    add("relaxed_mm1", 20, "markers_with_both_primers_mm1", int(marker_qc["has_both_primers_mm1"].sum()), "Оба праймера имеют full-length hit с mismatch <= 1")

    if relaxed_amplicons is not None and "marker_id" in relaxed_amplicons.columns:
        add("relaxed_mm1", 21, "candidate_amplicons", len(relaxed_amplicons), f"Файл: {relaxed_amplicons_file}")
        add("relaxed_mm1", 22, "markers_with_candidate_amplicon", int(relaxed_amplicons["marker_id"].nunique()), f"Файл: {relaxed_amplicons_file}")
        candidates_per_marker = relaxed_amplicons.groupby("marker_id").size()
        add("relaxed_mm1", 23, "markers_with_exactly_one_candidate_amplicon", int((candidates_per_marker == 1).sum()), "Однозначный relaxed mm1 ампликон")
    else:
        add("relaxed_mm1", 21, "candidate_amplicons", None, "Файл candidate amplicons не найден")
        add("relaxed_mm1", 22, "markers_with_candidate_amplicon", safe_nunique_marker_id(relaxed_raw_map_file), f"Оценка по файлу карты: {relaxed_raw_map_file}")
        add("relaxed_mm1", 23, "markers_with_exactly_one_candidate_amplicon", count_rows_tsv(relaxed_raw_map_file), f"Оценка по файлу карты: {relaxed_raw_map_file}")

    add("relaxed_mm1", 24, "relaxed_clean_unique", count_rows_tsv(relaxed_clean_unique_file), "После удаления scaffolds, LG/ref conflicts и duplicate positions")
    add("relaxed_mm1", 25, "relaxed_mm1_collinear", count_rows_tsv(FINAL_RELAXED_FILE), "Финальная relaxed mm1 collinear карта")

    funnel = pd.DataFrame(rows)
    return funnel.sort_values(["stage_order"])


def compute_gap_stats_for_map(map_file: Path, map_name: str):
    df = pd.read_csv(map_file, sep="\t")
    required = {"chr", "pos", "cM"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{map_file} не содержит колонки: {missing}")

    df = df.copy()
    df["chr"] = df["chr"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"])
    df["cM"] = pd.to_numeric(df["cM"])

    rows = []
    large_gap_rows = []

    for chrom, g in df.groupby("chr"):
        g = g.sort_values("pos").reset_index(drop=True)
        n = len(g)

        pos_gap = g["pos"].diff()
        cm_gap_abs = g["cM"].diff().abs()

        if n >= 3:
            spearman = g["pos"].rank().corr(g["cM"].rank())
        else:
            spearman = float("nan")

        orientation = "increasing" if pd.isna(spearman) or spearman >= 0 else "decreasing"

        rows.append({
            "map_name": map_name,
            "chr": chrom,
            "n_markers": n,
            "pos_min": int(g["pos"].min()),
            "pos_max": int(g["pos"].max()),
            "physical_span_bp": int(g["pos"].max() - g["pos"].min()) if n > 0 else 0,
            "cM_min": float(g["cM"].min()),
            "cM_max": float(g["cM"].max()),
            "genetic_span_cM": float(g["cM"].max() - g["cM"].min()) if n > 0 else 0,
            "spearman_pos_cM": spearman,
            "orientation": orientation,
            "n_intervals": max(n - 1, 0),
            "mean_physical_gap_bp": float(pos_gap.dropna().mean()) if n > 1 else None,
            "median_physical_gap_bp": float(pos_gap.dropna().median()) if n > 1 else None,
            "max_physical_gap_bp": int(pos_gap.dropna().max()) if n > 1 else None,
            "mean_genetic_gap_cM_abs": float(cm_gap_abs.dropna().mean()) if n > 1 else None,
            "median_genetic_gap_cM_abs": float(cm_gap_abs.dropna().median()) if n > 1 else None,
            "max_genetic_gap_cM_abs": float(cm_gap_abs.dropna().max()) if n > 1 else None,
        })

        for i in range(1, len(g)):
            physical_gap = int(g.loc[i, "pos"] - g.loc[i - 1, "pos"])
            genetic_gap_abs = float(abs(g.loc[i, "cM"] - g.loc[i - 1, "cM"]))

            if physical_gap >= 10_000_000 or genetic_gap_abs >= 10:
                large_gap_rows.append({
                    "map_name": map_name,
                    "chr": chrom,
                    "left_pos": int(g.loc[i - 1, "pos"]),
                    "right_pos": int(g.loc[i, "pos"]),
                    "physical_gap_bp": physical_gap,
                    "left_cM": float(g.loc[i - 1, "cM"]),
                    "right_cM": float(g.loc[i, "cM"]),
                    "genetic_gap_cM_abs": genetic_gap_abs,
                })

    return pd.DataFrame(rows), pd.DataFrame(large_gap_rows)


def write_summary_md(funnel, primer_qc, marker_qc, gap_stats, large_gaps):
    strict_final_n = int(funnel.loc[funnel["stage"] == "strict_high_confidence_collinear", "n"].iloc[0])
    relaxed_final_n = int(funnel.loc[funnel["stage"] == "relaxed_mm1_collinear", "n"].iloc[0])

    with open(SUMMARY_MD_OUT, "w", encoding="utf-8") as out:
        out.write("# QC summary: Brassica oleracea var. capitata\n\n")

        out.write("## Primer alignment QC\n\n")
        out.write(f"- Total primers: {len(primer_qc)}\n")
        out.write(f"- Primers with any BLAST hit: {int(primer_qc['has_any_blast_hit'].sum())}\n")
        out.write(f"- Primers with full-length perfect hit: {int(primer_qc['has_full_length_perfect_hit'].sum())}\n")
        out.write(f"- Primers with full-length <=1 mismatch hit: {int(primer_qc['has_full_length_mm1_hit'].sum())}\n")
        out.write(f"- Primers with multiple perfect full-length hits: {int(primer_qc['is_multihit_perfect'].sum())}\n")
        out.write(f"- Primers with multiple mm1 full-length hits: {int(primer_qc['is_multihit_mm1'].sum())}\n\n")

        out.write("## Marker-level QC\n\n")
        out.write(f"- Total markers: {len(marker_qc)}\n")
        out.write(f"- Markers with both primers perfect: {int(marker_qc['has_both_primers_perfect'].sum())}\n")
        out.write(f"- Markers with both primers <=1 mismatch: {int(marker_qc['has_both_primers_mm1'].sum())}\n")
        out.write(f"- Final strict high-confidence markers: {strict_final_n}\n")
        out.write(f"- Final relaxed mm1 collinear markers: {relaxed_final_n}\n\n")

        out.write("## Final maps\n\n")
        out.write("Final maps were checked for marker number, chromosome coverage, physical gaps and genetic gaps.\n\n")
        out.write(f"- Gap statistics table: `{GAP_STATS_OUT}`\n")
        out.write(f"- Large gaps table: `{LARGE_GAPS_OUT}`\n\n")

        out.write("## Output files\n\n")
        out.write(f"- Primer QC: `{PRIMER_QC_OUT}`\n")
        out.write(f"- Marker QC: `{MARKER_QC_OUT}`\n")
        out.write(f"- Filtering funnel: `{FUNNEL_OUT}`\n")
        out.write(f"- Gap statistics: `{GAP_STATS_OUT}`\n")
        out.write(f"- Large gaps: `{LARGE_GAPS_OUT}`\n")


def main():
    require_file(METADATA_FILE)
    require_file(BLAST_FILE)
    require_file(FINAL_STRICT_FILE)
    require_file(FINAL_RELAXED_FILE)

    metadata = pd.read_csv(METADATA_FILE, sep="\t")
    marker_col = get_metadata_marker_col(metadata)

    expected_primers = set()
    for marker_id in metadata[marker_col].astype(str):
        expected_primers.add(f"{marker_id}__F")
        expected_primers.add(f"{marker_id}__R")

    print("Scanning BLAST file. This may take a few minutes...")
    blast_stats = scan_blast_file(BLAST_FILE, expected_primers)

    primer_qc = make_primer_qc(metadata, blast_stats, marker_col)
    primer_qc.to_csv(PRIMER_QC_OUT, sep="\t", index=False)

    marker_qc = make_marker_qc(primer_qc)
    marker_qc = add_final_membership(marker_qc)
    marker_qc.to_csv(MARKER_QC_OUT, sep="\t", index=False)

    funnel = make_funnel(metadata, primer_qc, marker_qc)
    funnel.to_csv(FUNNEL_OUT, sep="\t", index=False)

    strict_gap_stats, strict_large_gaps = compute_gap_stats_for_map(
        FINAL_STRICT_FILE,
        "strict_high_confidence_196"
    )
    relaxed_gap_stats, relaxed_large_gaps = compute_gap_stats_for_map(
        FINAL_RELAXED_FILE,
        "relaxed_mm1_collinear_222"
    )

    gap_stats = pd.concat([strict_gap_stats, relaxed_gap_stats], ignore_index=True)
    large_gaps = pd.concat([strict_large_gaps, relaxed_large_gaps], ignore_index=True)

    gap_stats["chr_num"] = pd.to_numeric(gap_stats["chr"], errors="coerce")
    gap_stats = gap_stats.sort_values(["map_name", "chr_num", "chr"]).drop(columns=["chr_num"])
    gap_stats.to_csv(GAP_STATS_OUT, sep="\t", index=False)

    if len(large_gaps) > 0:
        large_gaps["chr_num"] = pd.to_numeric(large_gaps["chr"], errors="coerce")
        large_gaps = large_gaps.sort_values(["map_name", "chr_num", "chr", "left_pos"]).drop(columns=["chr_num"])
    large_gaps.to_csv(LARGE_GAPS_OUT, sep="\t", index=False)

    write_summary_md(funnel, primer_qc, marker_qc, gap_stats, large_gaps)

    print("\nDone.")
    print(f"Primer QC: {PRIMER_QC_OUT}")
    print(f"Marker QC: {MARKER_QC_OUT}")
    print(f"Filtering funnel: {FUNNEL_OUT}")
    print(f"Gap statistics: {GAP_STATS_OUT}")
    print(f"Large gaps: {LARGE_GAPS_OUT}")
    print(f"Summary: {SUMMARY_MD_OUT}")

    print("\nFiltering funnel:")
    print(funnel.to_string(index=False))

    print("\nFinal map gap statistics:")
    print(gap_stats.to_string(index=False))


if __name__ == "__main__":
    main()

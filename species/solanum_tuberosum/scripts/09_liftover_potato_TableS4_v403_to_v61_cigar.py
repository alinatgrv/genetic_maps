#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import sys
import pandas as pd


BASE = Path("species/solanum_tuberosum")

TABLES4 = BASE / "data/raw/pgsc_sharma_2013/TableS4.tsv"
PAF = BASE / "results/intermediate/liftover_v403_to_v61/potato_v403_to_v61.minimap2.asm5.paf"

OUT_INTERMEDIATE = BASE / "results/intermediate/liftover_v403_to_v61"
OUT_FINAL = BASE / "results/final"
OUT_QC = BASE / "results/qc"

OUT_INTERMEDIATE.mkdir(parents=True, exist_ok=True)
OUT_FINAL.mkdir(parents=True, exist_ok=True)
OUT_QC.mkdir(parents=True, exist_ok=True)

ALL_OUT = OUT_INTERMEDIATE / "potato_TableS4_liftover_v403_to_v61.cigar.all_markers.tsv"
CHR_AGREE_OUT = OUT_INTERMEDIATE / "potato_TableS4_liftover_v403_to_v61.cigar.chr_agree.tsv"

FINAL_WITH_MARKERS = OUT_FINAL / "potato_genetic_map.liftover_extended_cigar.with_markers.tsv"
FINAL_MAP = OUT_FINAL / "potato_genetic_map.liftover_extended_cigar.tsv"

SUMMARY_OUT = OUT_QC / "potato_liftover_extended_cigar_summary.txt"
QC_BY_CHR_OUT = OUT_QC / "potato_liftover_extended_cigar_qc_by_chr.tsv"
DIRECT_QC_OUT = OUT_QC / "potato_liftover_cigar_vs_direct_spuddb_solcap_qc.tsv"

STRICT_DIRECT = OUT_FINAL / "potato_genetic_map.with_markers.tsv"


def normalize_marker_id(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    x = x.replace("solcap_stsnp_", "solcap_snp_")
    return x.lower()


def extract_attr(attributes, keys):
    for key in keys:
        m = re.search(r"(?:^|[;\s])" + re.escape(key) + r"=([^;]+)", attributes)
        if m:
            return m.group(1).strip()
    return ""


def extract_genetic_position(attributes):
    m = re.search(
        r"Genetic position\s*=\s*Chr0?([0-9]+)\s+([0-9.]+)\s*cM",
        attributes,
        flags=re.IGNORECASE,
    )
    if not m:
        return None, None
    return int(m.group(1)), float(m.group(2))


def detect_delimiter(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#seqid") or line.startswith("seqid"):
                if line.count("\t") >= 8:
                    return "\t"
                return ","
    return "\t"


def read_tableS4(path):
    if not path.exists():
        sys.exit(f"ERROR: TableS4 not found: {path}")

    delimiter = detect_delimiter(path)
    rows = []
    seen_header = False

    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            if not row:
                continue

            first = row[0].strip()

            if first == "#seqid" or first == "seqid":
                seen_header = True
                continue

            if not seen_header:
                continue

            if first.startswith("#"):
                continue

            if len(row) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attributes = row[:9]

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue

            genetic_chr, cm = extract_genetic_position(attributes)
            if genetic_chr is None or cm is None:
                continue

            old_midpoint_1based = (start_i + end_i) // 2
            old_midpoint_0based = old_midpoint_1based - 1

            marker_id = extract_attr(attributes, ["Name", "NAME", "ID", "Alias"])
            if not marker_id:
                marker_id = f"{seqid}:{start_i}-{end_i}:{feature_type}"

            rows.append(
                {
                    "marker_id": marker_id,
                    "marker_id_norm": normalize_marker_id(marker_id),
                    "old_seqid": seqid,
                    "old_start": start_i,
                    "old_end": end_i,
                    "old_midpoint": old_midpoint_1based,
                    "old_midpoint_0based": old_midpoint_0based,
                    "old_source": source,
                    "marker_type": feature_type,
                    "old_strand": strand,
                    "genetic_chr": genetic_chr,
                    "cM": cm,
                    "attributes": attributes,
                }
            )

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(f"ERROR: no markers parsed from {path}")

    return df


def target_chr_to_num(tname):
    m = re.match(r"chr0?([0-9]+)$", str(tname))
    if not m:
        return None
    return int(m.group(1))


def get_cg_tag(fields):
    for x in fields[12:]:
        if x.startswith("cg:Z:"):
            return x[5:]
    return ""


def cigar_ops(cigar):
    return re.findall(r"(\d+)([MIDNSHP=X])", cigar)


def map_position_by_cigar(qpos0, qstart, qend, strand, tstart, tend, cigar):
    """
    Map one 0-based query coordinate to one 1-based target coordinate using PAF cg:Z CIGAR.

    Query is always traversed from qstart to qend.
    For '+' target is traversed from tstart to tend.
    For '-' target is traversed from tend-1 down to tstart.
    """
    q = int(qstart)

    if strand == "+":
        t = int(tstart)
        step = 1
    else:
        t = int(tend) - 1
        step = -1

    for length_s, op in cigar_ops(cigar):
        length = int(length_s)

        if op in ("M", "=", "X"):
            if q <= qpos0 < q + length:
                offset = qpos0 - q
                tpos0 = t + step * offset
                return int(tpos0) + 1, "cigar_aligned_base"

            q += length
            t += step * length

        elif op == "I":
            # Query consumes bases; target does not.
            # A marker here has no exact target base.
            if q <= qpos0 < q + length:
                return pd.NA, "query_insertion_in_old_not_in_new"
            q += length

        elif op in ("D", "N"):
            # Target consumes bases; query does not.
            t += step * length

        elif op == "S":
            if q <= qpos0 < q + length:
                return pd.NA, "query_softclip"
            q += length

        elif op == "H":
            continue

        elif op == "P":
            continue

        else:
            return pd.NA, f"unsupported_cigar_op_{op}"

    return pd.NA, "position_not_reached_in_cigar"


def choose_better(candidate, current):
    if current is None:
        return candidate

    cand_key = (
        candidate["same_chr"],
        candidate["cigar_map_ok"],
        candidate["mapq"],
        candidate["paf_identity"],
        candidate["paf_query_span"],
        candidate["n_match"],
    )
    curr_key = (
        current["same_chr"],
        current["cigar_map_ok"],
        current["mapq"],
        current["paf_identity"],
        current["paf_query_span"],
        current["n_match"],
    )
    return candidate if cand_key > curr_key else current


def liftover_with_paf_cigar(markers, paf_path):
    if not paf_path.exists():
        sys.exit(f"ERROR: PAF not found: {paf_path}")

    markers_by_query = {}
    for idx, row in markers.iterrows():
        markers_by_query.setdefault(row["old_seqid"], []).append(
            (idx, int(row["old_midpoint_0based"]), int(row["genetic_chr"]))
        )

    best = {idx: None for idx in markers.index}
    n_candidates = {idx: 0 for idx in markers.index}

    with open(paf_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            qname = fields[0]
            if qname not in markers_by_query:
                continue

            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tstart = int(fields[7])
            tend = int(fields[8])
            n_match = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])

            cigar = get_cg_tag(fields)
            if not cigar:
                continue

            target_chr_num = target_chr_to_num(tname)
            target_type = "chromosome" if target_chr_num is not None else "scaffold"

            if aln_len <= 0:
                continue

            identity = n_match / aln_len
            qspan = qend - qstart

            for idx, qpos0, genetic_chr in markers_by_query[qname]:
                if not (qstart <= qpos0 < qend):
                    continue

                target_pos, cigar_status = map_position_by_cigar(
                    qpos0=qpos0,
                    qstart=qstart,
                    qend=qend,
                    strand=strand,
                    tstart=tstart,
                    tend=tend,
                    cigar=cigar,
                )

                cigar_map_ok = not pd.isna(target_pos)
                same_chr = int(target_chr_num == genetic_chr) if target_chr_num is not None else 0

                candidate = {
                    "target_seqid": tname,
                    "target_chr": target_chr_num,
                    "target_type": target_type,
                    "target_pos": target_pos,
                    "paf_strand": strand,
                    "paf_qstart": qstart,
                    "paf_qend": qend,
                    "paf_tstart": tstart,
                    "paf_tend": tend,
                    "mapq": mapq,
                    "n_match": n_match,
                    "aln_len": aln_len,
                    "paf_identity": identity,
                    "paf_query_span": qspan,
                    "same_chr": same_chr,
                    "cigar_status": cigar_status,
                    "cigar_map_ok": int(cigar_map_ok),
                }

                n_candidates[idx] += 1
                best[idx] = choose_better(candidate, best[idx])

    lifted_rows = []

    for idx, row in markers.iterrows():
        rec = row.to_dict()
        rec["n_liftover_candidates"] = n_candidates[idx]

        if best[idx] is None:
            rec.update(
                {
                    "target_seqid": "",
                    "target_chr": pd.NA,
                    "target_type": "",
                    "target_pos": pd.NA,
                    "paf_strand": "",
                    "paf_qstart": pd.NA,
                    "paf_qend": pd.NA,
                    "paf_tstart": pd.NA,
                    "paf_tend": pd.NA,
                    "mapq": pd.NA,
                    "n_match": pd.NA,
                    "aln_len": pd.NA,
                    "paf_identity": pd.NA,
                    "paf_query_span": pd.NA,
                    "same_chr": 0,
                    "cigar_status": "no_covering_paf_block",
                    "cigar_map_ok": 0,
                    "liftover_status": "no_covering_paf_block",
                }
            )
        else:
            rec.update(best[idx])

            if rec["cigar_map_ok"] != 1:
                rec["liftover_status"] = "cigar_no_exact_target_base"
            elif rec["target_type"] != "chromosome":
                rec["liftover_status"] = "lifted_to_scaffold"
            elif int(rec["target_chr"]) != int(rec["genetic_chr"]):
                rec["liftover_status"] = "lifted_chr_disagree"
            else:
                rec["liftover_status"] = "lifted_chr_agree"

        lifted_rows.append(rec)

    return pd.DataFrame(lifted_rows)


def find_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def direct_spuddb_qc(lifted_chr_agree):
    if not STRICT_DIRECT.exists():
        return None, "Strict direct file not found"

    direct = pd.read_csv(STRICT_DIRECT, sep="\t")

    id_col = find_col(
        direct,
        [
            "marker_id_norm",
            "normalized_marker_id",
            "marker_id_normalized",
            "marker_id_for_join",
            "snp_id",
            "marker_id",
            "Name",
            "name",
            "ID",
        ],
    )
    chr_col = find_col(direct, ["chr", "target_chr", "physical_chr", "ref_chr", "chromosome"])
    pos_col = find_col(direct, ["pos", "target_pos", "physical_pos", "SNP_POS", "snp_pos"])

    if id_col is None or chr_col is None or pos_col is None:
        msg = f"Could not identify columns in {STRICT_DIRECT}. Columns: {list(direct.columns)}"
        return None, msg

    direct = direct.copy()
    direct["marker_id_norm_for_join"] = direct[id_col].map(normalize_marker_id)
    direct["direct_chr"] = direct[chr_col].astype(str).str.replace("chr", "", regex=False).astype(int)
    direct["direct_pos"] = direct[pos_col].astype(int)

    liftover = lifted_chr_agree.copy()
    liftover["marker_id_norm_for_join"] = liftover["marker_id_norm"].map(normalize_marker_id)

    merged = liftover.merge(
        direct[["marker_id_norm_for_join", "direct_chr", "direct_pos"]],
        on="marker_id_norm_for_join",
        how="inner",
    )

    if merged.empty:
        return merged, "No marker overlap between liftover and strict direct map"

    merged["liftover_chr"] = merged["target_chr"].astype(int)
    merged["liftover_pos"] = merged["target_pos"].astype(int)
    merged["same_chr_direct"] = merged["liftover_chr"] == merged["direct_chr"]
    merged["delta_pos"] = merged["liftover_pos"] - merged["direct_pos"]
    merged["abs_delta_pos"] = merged["delta_pos"].abs()

    return merged, "OK"


def main():
    print(f"Reading TableS4: {TABLES4}")
    markers = read_tableS4(TABLES4)
    print(f"Parsed markers: {len(markers)}")

    print(f"Reading PAF and lifting coordinates by cg:Z CIGAR: {PAF}")
    lifted = liftover_with_paf_cigar(markers, PAF)

    lifted.to_csv(ALL_OUT, sep="\t", index=False)

    chr_agree = lifted[lifted["liftover_status"] == "lifted_chr_agree"].copy()
    chr_agree["chr"] = chr_agree["target_chr"].astype(int)
    chr_agree["pos"] = chr_agree["target_pos"].astype(int)

    chr_agree = chr_agree.sort_values(["chr", "pos", "cM", "marker_id"])
    chr_agree.to_csv(CHR_AGREE_OUT, sep="\t", index=False)

    final_with_markers = chr_agree[
        [
            "chr",
            "pos",
            "cM",
            "marker_id",
            "marker_type",
            "genetic_chr",
            "old_seqid",
            "old_start",
            "old_end",
            "old_midpoint",
            "target_seqid",
            "target_pos",
            "paf_strand",
            "mapq",
            "paf_identity",
            "paf_query_span",
            "n_liftover_candidates",
            "cigar_status",
            "liftover_status",
        ]
    ].copy()

    final_with_markers.to_csv(FINAL_WITH_MARKERS, sep="\t", index=False)

    final_map = final_with_markers[["chr", "pos", "cM"]].drop_duplicates()
    final_map = final_map.sort_values(["chr", "pos", "cM"])
    final_map.to_csv(FINAL_MAP, sep="\t", index=False)

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
    qc_by_chr.to_csv(QC_BY_CHR_OUT, sep="\t", index=False)

    direct_qc, direct_msg = direct_spuddb_qc(chr_agree)
    if direct_qc is not None and not direct_qc.empty:
        direct_qc.to_csv(DIRECT_QC_OUT, sep="\t", index=False)

    with open(SUMMARY_OUT, "w", encoding="utf-8") as out:
        out.write("Potato TableS4 CIGAR liftover from PGSC/DM v4.03 to DM v6.1\n")
        out.write("===========================================================\n\n")
        out.write(f"Input TableS4 markers: {len(markers)}\n")
        out.write(f"PAF: {PAF}\n\n")

        out.write("Liftover status counts:\n")
        out.write(lifted["liftover_status"].value_counts(dropna=False).to_string())
        out.write("\n\n")

        out.write("CIGAR status counts:\n")
        out.write(lifted["cigar_status"].value_counts(dropna=False).to_string())
        out.write("\n\n")

        out.write("Marker type counts in input TableS4:\n")
        out.write(markers["marker_type"].value_counts(dropna=False).to_string())
        out.write("\n\n")

        out.write(f"Lifted chr-agree marker rows: {len(final_with_markers)}\n")
        out.write(f"Final chr-pos-cM rows: {len(final_map)}\n\n")

        out.write("Final rows by chromosome:\n")
        if not qc_by_chr.empty:
            out.write(qc_by_chr.to_string(index=False))
        out.write("\n\n")

        out.write("Direct SpudDB SolCAP QC:\n")
        out.write(f"{direct_msg}\n")
        if direct_qc is not None and not direct_qc.empty:
            out.write(f"Overlapping SolCAP markers: {len(direct_qc)}\n")
            out.write(f"Same chromosome as direct SpudDB: {int(direct_qc['same_chr_direct'].sum())}\n")
            out.write("Absolute position delta quantiles, bp:\n")
            out.write(
                direct_qc["abs_delta_pos"]
                .quantile([0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0])
                .to_string()
            )
            out.write("\n")

    print("Done.")
    print(f"All lifted markers:        {ALL_OUT}")
    print(f"Chr-agree lifted markers:  {CHR_AGREE_OUT}")
    print(f"Final CIGAR extended map:  {FINAL_MAP}")
    print(f"Final CIGAR with IDs:      {FINAL_WITH_MARKERS}")
    print(f"QC by chr:                 {QC_BY_CHR_OUT}")
    print(f"Summary:                   {SUMMARY_OUT}")
    if direct_qc is not None and not direct_qc.empty:
        print(f"Direct SolCAP QC:          {DIRECT_QC_OUT}")


if __name__ == "__main__":
    main()

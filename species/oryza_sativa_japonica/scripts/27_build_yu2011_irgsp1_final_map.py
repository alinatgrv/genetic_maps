#!/usr/bin/env python3

from pathlib import Path
import csv
import re
from collections import Counter, defaultdict

ORIGINAL_BED = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.fasta_seqids.bed")
LIFTED_BED = Path("species/oryza_sativa_japonica/results/liftover/yu2011_tigr6_to_irgsp1/yu2011_bins_irgsp1_lifted.bed")

FINAL_MAP = Path("species/oryza_sativa_japonica/results/final/oryza_sativa_japonica_genetic_map.tsv")
DETAILS = Path("species/oryza_sativa_japonica/results/final/oryza_sativa_japonica_genetic_map.details.tsv")
QC_SUMMARY = Path("species/oryza_sativa_japonica/results/qc/yu2011_irgsp1_liftover_final_map_summary.txt")
QC_BY_CHR = Path("species/oryza_sativa_japonica/results/qc/yu2011_irgsp1_liftover_final_map_by_chr.tsv")

FINAL_MAP.parent.mkdir(parents=True, exist_ok=True)
DETAILS.parent.mkdir(parents=True, exist_ok=True)
QC_SUMMARY.parent.mkdir(parents=True, exist_ok=True)

# IRGSP-1.0 / GCF_001433935.1 primary chromosome accessions.
target_seqid_to_chr = {
    "NC_029256.1": "1",
    "NC_029257.1": "2",
    "NC_029258.1": "3",
    "NC_029259.1": "4",
    "NC_029260.1": "5",
    "NC_029261.1": "6",
    "NC_029262.1": "7",
    "NC_029263.1": "8",
    "NC_029264.1": "9",
    "NC_029265.1": "10",
    "NC_029266.1": "11",
    "NC_029267.1": "12",
}

def parse_original_name(name):
    # Example:
    # Bin1|chr=1|cM=0.000000|old_mid=282464
    parts = name.split("|")
    bin_id = parts[0]
    d = {"bin": bin_id}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            d[k] = v
    return d

# Map original TIGR6 BED coordinate to bin metadata.
original = {}
with ORIGINAL_BED.open() as f:
    for line in f:
        if not line.strip():
            continue
        seqid, start, end, name = line.rstrip("\n").split("\t")[:4]
        meta = parse_original_name(name)
        key = (seqid, int(start), int(end))
        original[key] = {
            "tigr6_seqid": seqid,
            "tigr6_start0": int(start),
            "tigr6_end0": int(end),
            "bin": meta["bin"],
            "source_chr": meta["chr"],
            "cM": meta["cM"],
            "old_mid": meta["old_mid"],
            "original_name": name,
        }

rows = []
unmatched = []
target_unmapped = []

with LIFTED_BED.open() as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise SystemExit(f"Unexpected lifted BED line: {line}")

        target_seqid = fields[0]
        target_start0 = int(fields[1])
        target_end0 = int(fields[2])
        query_key = fields[3]
        strand = fields[5]

        # paftools uses query_id_queryStart_queryEnd in column 4.
        try:
            qseqid, qstart, qend = query_key.rsplit("_", 2)
            qstart = int(qstart)
            qend = int(qend)
        except Exception:
            raise SystemExit(f"Could not parse lifted query key: {query_key}")

        key = (qseqid, qstart, qend)
        if key not in original:
            unmatched.append((query_key, line.rstrip("\n")))
            continue

        meta = original[key]

        if target_seqid not in target_seqid_to_chr:
            target_unmapped.append((target_seqid, line.rstrip("\n")))
            continue

        # BED is 0-based half-open. For 1-bp interval, 1-based coordinate = start + 1 = end.
        target_pos_1based = target_start0 + 1
        target_chr = target_seqid_to_chr[target_seqid]

        rows.append({
            "chr": target_chr,
            "pos": target_pos_1based,
            "cM": float(meta["cM"]),
            "bin": meta["bin"],
            "source_chr": meta["source_chr"],
            "tigr6_seqid": meta["tigr6_seqid"],
            "tigr6_pos": int(meta["old_mid"]),
            "irgsp1_seqid": target_seqid,
            "irgsp1_start0": target_start0,
            "irgsp1_end0": target_end0,
            "strand": strand,
        })

if unmatched:
    raise SystemExit(f"Unmatched lifted rows: {len(unmatched)}. First: {unmatched[0]}")
if target_unmapped:
    raise SystemExit(f"Rows on non-primary/unmapped target seqids: {len(target_unmapped)}. First: {target_unmapped[0]}")

# Sort in natural chromosome order and then physical position.
rows.sort(key=lambda r: (int(r["chr"]), int(r["pos"]), float(r["cM"]), r["bin"]))

# Final minimal map.
with FINAL_MAP.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=["chr", "pos", "cM"])
    writer.writeheader()
    for r in rows:
        writer.writerow({
            "chr": r["chr"],
            "pos": r["pos"],
            "cM": f'{r["cM"]:.6f}',
        })

# Detailed map for traceability.
with DETAILS.open("w", newline="") as w:
    fieldnames = [
        "chr", "pos", "cM", "bin",
        "source_chr", "tigr6_seqid", "tigr6_pos",
        "irgsp1_seqid", "irgsp1_start0", "irgsp1_end0", "strand",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for r in rows:
        out = dict(r)
        out["cM"] = f'{r["cM"]:.6f}'
        writer.writerow(out)

# QC.
by_chr = defaultdict(list)
for r in rows:
    by_chr[r["chr"]].append(r)

with QC_BY_CHR.open("w", newline="") as w:
    fieldnames = [
        "chr",
        "n_bins",
        "pos_min",
        "pos_max",
        "cM_min",
        "cM_max",
        "same_chr_as_source",
        "plus_strand",
        "minus_strand",
        "duplicate_positions",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for c in sorted(by_chr, key=lambda x: int(x)):
        rs = by_chr[c]
        pos_counts = Counter(r["pos"] for r in rs)
        writer.writerow({
            "chr": c,
            "n_bins": len(rs),
            "pos_min": min(r["pos"] for r in rs),
            "pos_max": max(r["pos"] for r in rs),
            "cM_min": f'{min(r["cM"] for r in rs):.6f}',
            "cM_max": f'{max(r["cM"] for r in rs):.6f}',
            "same_chr_as_source": sum(1 for r in rs if r["chr"] == r["source_chr"]),
            "plus_strand": sum(1 for r in rs if r["strand"] == "+"),
            "minus_strand": sum(1 for r in rs if r["strand"] == "-"),
            "duplicate_positions": sum(n - 1 for n in pos_counts.values() if n > 1),
        })

total_duplicates = sum(n - 1 for n in Counter((r["chr"], r["pos"]) for r in rows).values() if n > 1)
same_chr_total = sum(1 for r in rows if r["chr"] == r["source_chr"])
strand_counts = Counter(r["strand"] for r in rows)

with QC_SUMMARY.open("w") as w:
    w.write(f"original_bed_rows\t{len(original)}\n")
    w.write(f"lifted_rows\t{len(rows)}\n")
    w.write(f"final_map_rows\t{len(rows)}\n")
    w.write(f"chromosomes\t{','.join(sorted(by_chr, key=lambda x: int(x)))}\n")
    w.write(f"same_chr_as_source\t{same_chr_total}\n")
    w.write(f"duplicate_chr_pos\t{total_duplicates}\n")
    for strand, n in sorted(strand_counts.items()):
        w.write(f"strand_{strand}\t{n}\n")
    w.write(f"final_map\t{FINAL_MAP}\n")
    w.write(f"details\t{DETAILS}\n")
    w.write(f"qc_by_chr\t{QC_BY_CHR}\n")

print(f"Written final map: {FINAL_MAP}")
print(f"Written details:   {DETAILS}")
print(f"Written QC:        {QC_SUMMARY}")
print(f"Written QC table:  {QC_BY_CHR}")
print(f"Rows: {len(rows)}")
print(f"Same chromosome as source: {same_chr_total}/{len(rows)}")
print(f"Duplicate chr:pos: {total_duplicates}")
print(f"Strands: {dict(strand_counts)}")

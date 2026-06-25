#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import Counter, defaultdict

ORIGINAL_BED = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.fasta_seqids.bed")
LIFTED_BED = Path("species/oryza_sativa_indica/results/liftover/yu2011_tigr6_to_asm465v1/yu2011_bins_asm465v1_lifted.bed")

FINAL_MAP = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.candidate_yu2011_projection.tsv")
DETAILS = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.candidate_yu2011_projection.details.tsv")

QC_SUMMARY = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_projection_summary.txt")
QC_BY_CHR = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_projection_by_chr.tsv")
UNLIFTED = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_unlifted_bins.tsv")
NONPRIMARY = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_nonprimary_lifted_rows.tsv")

for p in [FINAL_MAP, DETAILS, QC_SUMMARY, QC_BY_CHR, UNLIFTED, NONPRIMARY]:
    p.parent.mkdir(parents=True, exist_ok=True)

primary_target_seqids = {str(i) for i in range(1, 13)}

def parse_original_name(name):
    # Bin1|chr=1|cM=0.000000|old_mid=282464
    parts = name.split("|")
    meta = {"bin": parts[0]}
    for p in parts[1:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k] = v
    return meta

original = {}
with ORIGINAL_BED.open() as f:
    for line in f:
        seqid, start, end, name = line.rstrip("\n").split("\t")[:4]
        meta = parse_original_name(name)
        key = (seqid, int(start), int(end))
        original[key] = {
            "tigr6_seqid": seqid,
            "tigr6_start0": int(start),
            "tigr6_end0": int(end),
            "bin": meta["bin"],
            "source_chr": meta["chr"],
            "cM": float(meta["cM"]),
            "old_mid": int(meta["old_mid"]),
            "original_name": name,
        }

lifted_keys = set()
rows = []
nonprimary_rows = []
unmatched = []

with LIFTED_BED.open() as f:
    for line in f:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        target_seqid = fields[0]
        target_start0 = int(fields[1])
        target_end0 = int(fields[2])
        query_key = fields[3]
        strand = fields[5] if len(fields) > 5 else "."

        qseqid, qstart, qend = query_key.rsplit("_", 2)
        qstart = int(qstart)
        qend = int(qend)
        key = (qseqid, qstart, qend)
        lifted_keys.add(key)

        if key not in original:
            unmatched.append(line.rstrip("\n"))
            continue

        meta = original[key]

        if target_seqid not in primary_target_seqids:
            nonprimary_rows.append({
                "target_seqid": target_seqid,
                "target_start0": target_start0,
                "target_end0": target_end0,
                "query_key": query_key,
                "strand": strand,
                **meta,
            })
            continue

        target_chr = target_seqid
        target_pos_1based = target_start0 + 1

        rows.append({
            "chr": target_chr,
            "pos": target_pos_1based,
            "cM": meta["cM"],
            "bin": meta["bin"],
            "source_chr": meta["source_chr"],
            "tigr6_seqid": meta["tigr6_seqid"],
            "tigr6_pos": meta["old_mid"],
            "asm465v1_seqid": target_seqid,
            "asm465v1_start0": target_start0,
            "asm465v1_end0": target_end0,
            "strand": strand,
        })

if unmatched:
    raise SystemExit(f"Unmatched lifted rows: {len(unmatched)}. First: {unmatched[0]}")

unlifted_keys = sorted(set(original) - lifted_keys, key=lambda x: (x[0], x[1], x[2]))

# Write unlifted bins.
with UNLIFTED.open("w", newline="") as w:
    fieldnames = ["bin", "source_chr", "cM", "tigr6_seqid", "tigr6_pos", "tigr6_start0", "tigr6_end0"]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for key in unlifted_keys:
        m = original[key]
        writer.writerow({
            "bin": m["bin"],
            "source_chr": m["source_chr"],
            "cM": f'{m["cM"]:.6f}',
            "tigr6_seqid": m["tigr6_seqid"],
            "tigr6_pos": m["old_mid"],
            "tigr6_start0": m["tigr6_start0"],
            "tigr6_end0": m["tigr6_end0"],
        })

# Write nonprimary rows if any.
with NONPRIMARY.open("w", newline="") as w:
    fieldnames = [
        "target_seqid", "target_start0", "target_end0", "query_key", "strand",
        "bin", "source_chr", "cM", "tigr6_seqid", "old_mid"
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for r in nonprimary_rows:
        writer.writerow({
            "target_seqid": r["target_seqid"],
            "target_start0": r["target_start0"],
            "target_end0": r["target_end0"],
            "query_key": r["query_key"],
            "strand": r["strand"],
            "bin": r["bin"],
            "source_chr": r["source_chr"],
            "cM": f'{r["cM"]:.6f}',
            "tigr6_seqid": r["tigr6_seqid"],
            "old_mid": r["old_mid"],
        })

# Sort candidate map.
rows.sort(key=lambda r: (int(r["chr"]), int(r["pos"]), r["cM"], r["bin"]))

with FINAL_MAP.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=["chr", "pos", "cM"])
    writer.writeheader()
    for r in rows:
        writer.writerow({
            "chr": r["chr"],
            "pos": r["pos"],
            "cM": f'{r["cM"]:.6f}',
        })

with DETAILS.open("w", newline="") as w:
    fieldnames = [
        "chr", "pos", "cM", "bin", "source_chr",
        "tigr6_seqid", "tigr6_pos",
        "asm465v1_seqid", "asm465v1_start0", "asm465v1_end0", "strand",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for r in rows:
        out = dict(r)
        out["cM"] = f'{r["cM"]:.6f}'
        writer.writerow(out)

by_chr = defaultdict(list)
for r in rows:
    by_chr[r["chr"]].append(r)

with QC_BY_CHR.open("w", newline="") as w:
    fieldnames = [
        "chr", "n_bins", "pos_min", "pos_max", "cM_min", "cM_max",
        "same_chr_as_source", "plus_strand", "minus_strand",
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
target_seqid_counts = Counter(r["asm465v1_seqid"] for r in rows)

with QC_SUMMARY.open("w") as w:
    w.write(f"original_bed_rows\t{len(original)}\n")
    w.write(f"lifted_rows_total\t{len(lifted_keys)}\n")
    w.write(f"lifted_rows_primary_1_12\t{len(rows)}\n")
    w.write(f"nonprimary_lifted_rows\t{len(nonprimary_rows)}\n")
    w.write(f"unlifted_rows\t{len(unlifted_keys)}\n")
    w.write(f"candidate_map_rows\t{len(rows)}\n")
    w.write(f"chromosomes\t{','.join(sorted(by_chr, key=lambda x: int(x)))}\n")
    w.write(f"same_chr_as_source\t{same_chr_total}\n")
    w.write(f"duplicate_chr_pos\t{total_duplicates}\n")
    for strand, n in sorted(strand_counts.items()):
        w.write(f"strand_{strand}\t{n}\n")
    for seqid, n in sorted(target_seqid_counts.items(), key=lambda x: int(x[0])):
        w.write(f"target_seqid_{seqid}\t{n}\n")
    w.write(f"candidate_map\t{FINAL_MAP}\n")
    w.write(f"details\t{DETAILS}\n")
    w.write(f"unlifted\t{UNLIFTED}\n")
    w.write(f"nonprimary\t{NONPRIMARY}\n")
    w.write(f"qc_by_chr\t{QC_BY_CHR}\n")

print(f"Written candidate map: {FINAL_MAP}")
print(f"Written details:       {DETAILS}")
print(f"Written QC summary:    {QC_SUMMARY}")
print(f"Written QC by chr:     {QC_BY_CHR}")
print(f"Written unlifted:      {UNLIFTED}")
print(f"Written nonprimary:    {NONPRIMARY}")
print(f"Original rows: {len(original)}")
print(f"Lifted rows total: {len(lifted_keys)}")
print(f"Primary lifted rows: {len(rows)}")
print(f"Unlifted rows: {len(unlifted_keys)}")
print(f"Non-primary rows: {len(nonprimary_rows)}")
print(f"Same chromosome as source: {same_chr_total}/{len(rows)}")
print(f"Duplicate chr:pos: {total_duplicates}")
print(f"Strands: {dict(strand_counts)}")

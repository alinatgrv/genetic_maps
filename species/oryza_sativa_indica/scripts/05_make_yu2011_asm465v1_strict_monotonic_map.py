#!/usr/bin/env python3

from pathlib import Path
import csv
from bisect import bisect_left
from collections import defaultdict, Counter

IN = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.candidate_yu2011_projection.details.tsv")

OUT = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.strict_monotonic_yu2011_projection.tsv")
DETAILS = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.strict_monotonic_yu2011_projection.details.tsv")
EXCLUDED = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_strict_monotonic_excluded_bins.tsv")
QC = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_strict_monotonic_summary.txt")
QC_BY_CHR = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_strict_monotonic_by_chr.tsv")

for p in [OUT, DETAILS, EXCLUDED, QC, QC_BY_CHR]:
    p.parent.mkdir(parents=True, exist_ok=True)

def lis_indices_strict(values):
    """
    Return indices of one longest strictly increasing subsequence.
    values are compared as floats.
    """
    if not values:
        return []

    tails = []
    tails_idx = []
    prev = [-1] * len(values)

    for i, v in enumerate(values):
        j = bisect_left(tails, v)

        if j == len(tails):
            tails.append(v)
            tails_idx.append(i)
        else:
            tails[j] = v
            tails_idx[j] = i

        if j > 0:
            prev[i] = tails_idx[j - 1]

    keep = []
    k = tails_idx[-1]
    while k != -1:
        keep.append(k)
        k = prev[k]

    return list(reversed(keep))

rows = []
fieldnames = None

with IN.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    fieldnames = reader.fieldnames

    for r in reader:
        # Strict filter 1: only same chromosome as source.
        if r["chr"] != r["source_chr"]:
            continue

        rows.append(r)

by_chr = defaultdict(list)
for r in rows:
    by_chr[r["chr"]].append(r)

kept = []
excluded = []

for chrom in sorted(by_chr, key=lambda x: int(x)):
    # Sort by physical position on ASM465v1.
    chrom_rows = sorted(
        by_chr[chrom],
        key=lambda r: (int(r["pos"]), float(r["cM"]), r["bin"])
    )

    values = [float(r["cM"]) for r in chrom_rows]
    keep_idx = set(lis_indices_strict(values))

    for i, r in enumerate(chrom_rows):
        if i in keep_idx:
            kept.append(r)
        else:
            rr = dict(r)
            rr["exclude_reason"] = "breaks_monotonic_cM_order_after_liftover"
            excluded.append(rr)

kept.sort(key=lambda r: (int(r["chr"]), int(r["pos"]), float(r["cM"]), r["bin"]))
excluded.sort(key=lambda r: (int(r["chr"]), int(r["pos"]), float(r["cM"]), r["bin"]))

with OUT.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=["chr", "pos", "cM"])
    writer.writeheader()
    for r in kept:
        writer.writerow({
            "chr": r["chr"],
            "pos": r["pos"],
            "cM": r["cM"],
        })

with DETAILS.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()
    for r in kept:
        writer.writerow(r)

excluded_fieldnames = fieldnames + ["exclude_reason"]
with EXCLUDED.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=excluded_fieldnames)
    writer.writeheader()
    for r in excluded:
        writer.writerow(r)

# QC by chromosome after monotonic filtering.
final_by_chr = defaultdict(list)
for r in kept:
    final_by_chr[r["chr"]].append(r)

with QC_BY_CHR.open("w", newline="") as w:
    fieldnames_qc = [
        "chr",
        "n_bins",
        "pos_min",
        "pos_max",
        "cM_min",
        "cM_max",
        "increasing_cM_steps",
        "decreasing_cM_steps",
        "equal_cM_steps",
        "minus_strand_rows",
        "duplicate_positions",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames_qc)
    writer.writeheader()

    for chrom in sorted(final_by_chr, key=lambda x: int(x)):
        rs = sorted(final_by_chr[chrom], key=lambda r: int(r["pos"]))

        inc = dec = eq = 0
        for a, b in zip(rs, rs[1:]):
            if float(b["cM"]) > float(a["cM"]):
                inc += 1
            elif float(b["cM"]) < float(a["cM"]):
                dec += 1
            else:
                eq += 1

        pos_counts = Counter(r["pos"] for r in rs)

        writer.writerow({
            "chr": chrom,
            "n_bins": len(rs),
            "pos_min": min(int(r["pos"]) for r in rs),
            "pos_max": max(int(r["pos"]) for r in rs),
            "cM_min": f'{min(float(r["cM"]) for r in rs):.6f}',
            "cM_max": f'{max(float(r["cM"]) for r in rs):.6f}',
            "increasing_cM_steps": inc,
            "decreasing_cM_steps": dec,
            "equal_cM_steps": eq,
            "minus_strand_rows": sum(1 for r in rs if r["strand"] == "-"),
            "duplicate_positions": sum(n - 1 for n in pos_counts.values() if n > 1),
        })

total_duplicates = sum(
    n - 1 for n in Counter((r["chr"], r["pos"]) for r in kept).values()
    if n > 1
)

with QC.open("w") as w:
    w.write(f"input_details\t{IN}\n")
    w.write(f"same_chr_input_rows\t{len(rows)}\n")
    w.write(f"kept_monotonic_rows\t{len(kept)}\n")
    w.write(f"excluded_nonmonotonic_rows\t{len(excluded)}\n")
    w.write(f"duplicate_chr_pos\t{total_duplicates}\n")
    w.write(f"output_map\t{OUT}\n")
    w.write(f"output_details\t{DETAILS}\n")
    w.write(f"excluded_bins\t{EXCLUDED}\n")
    w.write(f"qc_by_chr\t{QC_BY_CHR}\n")

print(f"Written: {OUT}")
print(f"Written: {DETAILS}")
print(f"Written: {EXCLUDED}")
print(f"Written: {QC}")
print(f"Written: {QC_BY_CHR}")
print(f"Same-chr input rows: {len(rows)}")
print(f"Kept monotonic rows: {len(kept)}")
print(f"Excluded nonmonotonic rows: {len(excluded)}")
print(f"Duplicate chr:pos: {total_duplicates}")

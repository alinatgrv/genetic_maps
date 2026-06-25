#!/usr/bin/env python3

from pathlib import Path
import csv

IN = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.candidate_yu2011_projection.details.tsv")

OUT = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.strict_yu2011_projection.tsv")
DETAILS = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.strict_yu2011_projection.details.tsv")
QC = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_strict_same_chr_summary.txt")

rows = []
with IN.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r["chr"] != r["source_chr"]:
            continue
        rows.append(r)

rows.sort(key=lambda r: (int(r["chr"]), int(r["pos"]), float(r["cM"]), r["bin"]))

with OUT.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=["chr", "pos", "cM"])
    writer.writeheader()
    for r in rows:
        writer.writerow({
            "chr": r["chr"],
            "pos": r["pos"],
            "cM": r["cM"],
        })

with DETAILS.open("w", newline="") as w:
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=reader.fieldnames)
    writer.writeheader()
    for r in rows:
        writer.writerow(r)

with QC.open("w") as w:
    w.write(f"input_details\t{IN}\n")
    w.write(f"strict_rows\t{len(rows)}\n")
    w.write(f"output_map\t{OUT}\n")
    w.write(f"output_details\t{DETAILS}\n")

print(f"Written: {OUT}")
print(f"Written: {DETAILS}")
print(f"Written: {QC}")
print(f"Rows: {len(rows)}")

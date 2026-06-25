#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import defaultdict

IN = Path("species/oryza_sativa_japonica/results/final/oryza_sativa_japonica_genetic_map.details.tsv")
OUT = Path("species/oryza_sativa_japonica/results/qc/yu2011_irgsp1_final_map_collinearity_qc.tsv")

by_chr = defaultdict(list)

with IN.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        by_chr[r["chr"]].append({
            "chr": r["chr"],
            "pos": int(r["pos"]),
            "cM": float(r["cM"]),
            "bin": r["bin"],
            "strand": r["strand"],
        })

OUT.parent.mkdir(parents=True, exist_ok=True)

with OUT.open("w", newline="") as w:
    fieldnames = [
        "chr",
        "n_bins",
        "pos_min",
        "pos_max",
        "cM_min",
        "cM_max",
        "increasing_cM_steps",
        "decreasing_cM_steps",
        "equal_cM_steps",
        "decreasing_examples",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for chrom in sorted(by_chr, key=lambda x: int(x)):
        rows = sorted(by_chr[chrom], key=lambda r: r["pos"])

        inc = 0
        dec = 0
        eq = 0
        examples = []

        for a, b in zip(rows, rows[1:]):
            if b["cM"] > a["cM"]:
                inc += 1
            elif b["cM"] < a["cM"]:
                dec += 1
                if len(examples) < 5:
                    examples.append(
                        f'{a["bin"]}:{a["pos"]}:{a["cM"]}->'
                        f'{b["bin"]}:{b["pos"]}:{b["cM"]}'
                    )
            else:
                eq += 1

        writer.writerow({
            "chr": chrom,
            "n_bins": len(rows),
            "pos_min": min(r["pos"] for r in rows),
            "pos_max": max(r["pos"] for r in rows),
            "cM_min": f'{min(r["cM"] for r in rows):.6f}',
            "cM_max": f'{max(r["cM"] for r in rows):.6f}',
            "increasing_cM_steps": inc,
            "decreasing_cM_steps": dec,
            "equal_cM_steps": eq,
            "decreasing_examples": ";".join(examples),
        })

print(f"Written: {OUT}")

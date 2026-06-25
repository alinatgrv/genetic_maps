#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import defaultdict

IN = Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.candidate_yu2011_projection.details.tsv")
OUT_ALL = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_projection_collinearity_all_primary.tsv")
OUT_SAME = Path("species/oryza_sativa_indica/results/qc/yu2011_asm465v1_projection_collinearity_same_chr_only.tsv")

def load_rows(same_chr_only=False):
    by_chr = defaultdict(list)
    with IN.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if same_chr_only and r["chr"] != r["source_chr"]:
                continue
            by_chr[r["chr"]].append({
                "chr": r["chr"],
                "pos": int(r["pos"]),
                "cM": float(r["cM"]),
                "bin": r["bin"],
                "source_chr": r["source_chr"],
                "strand": r["strand"],
            })
    return by_chr

def write_qc(by_chr, out):
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w", newline="") as w:
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
            "cross_chr_rows",
            "minus_strand_rows",
            "decreasing_examples",
        ]
        writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for chrom in sorted(by_chr, key=lambda x: int(x)):
            rows = sorted(by_chr[chrom], key=lambda r: r["pos"])

            inc = dec = eq = 0
            examples = []

            for a, b in zip(rows, rows[1:]):
                if b["cM"] > a["cM"]:
                    inc += 1
                elif b["cM"] < a["cM"]:
                    dec += 1
                    if len(examples) < 8:
                        examples.append(
                            f'{a["bin"]}:{a["pos"]}:{a["cM"]}:src{a["source_chr"]}->{b["bin"]}:{b["pos"]}:{b["cM"]}:src{b["source_chr"]}'
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
                "cross_chr_rows": sum(1 for r in rows if r["chr"] != r["source_chr"]),
                "minus_strand_rows": sum(1 for r in rows if r["strand"] == "-"),
                "decreasing_examples": ";".join(examples),
            })

write_qc(load_rows(same_chr_only=False), OUT_ALL)
write_qc(load_rows(same_chr_only=True), OUT_SAME)

print(f"Written: {OUT_ALL}")
print(f"Written: {OUT_SAME}")

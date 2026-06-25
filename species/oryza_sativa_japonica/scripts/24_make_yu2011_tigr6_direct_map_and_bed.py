#!/usr/bin/env python3

from pathlib import Path
import csv

IN = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6.tsv")

DIRECT_MAP = Path("species/oryza_sativa_japonica/results/liftover/yu2011_tigr6_to_irgsp1/yu2011_genetic_map_tigr6_direct.tsv")
BED = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.bed")

DIRECT_MAP.parent.mkdir(parents=True, exist_ok=True)
BED.parent.mkdir(parents=True, exist_ok=True)

with IN.open() as f, DIRECT_MAP.open("w", newline="") as map_out, BED.open("w") as bed_out:
    reader = csv.DictReader(f, delimiter="\t")

    map_writer = csv.DictWriter(
        map_out,
        delimiter="\t",
        fieldnames=["chr", "pos", "cM", "bin", "source"]
    )
    map_writer.writeheader()

    for r in reader:
        chr_num = r["chr"]
        mid_1based = int(r["old_mid_bp_1based"])
        cm = r["cM"]

        map_writer.writerow({
            "chr": chr_num,
            "pos": mid_1based,
            "cM": cm,
            "bin": r["bin"],
            "source": "Yu2011_TIGR6.1_midpoint_not_lifted",
        })

        bed_start = max(0, mid_1based - 1)
        bed_end = mid_1based
        name = f'{r["bin"]}|chr={chr_num}|cM={cm}|old_mid={mid_1based}'
        bed_out.write(f"{chr_num}\t{bed_start}\t{bed_end}\t{name}\n")

print(f"Written: {DIRECT_MAP}")
print(f"Written: {BED}")

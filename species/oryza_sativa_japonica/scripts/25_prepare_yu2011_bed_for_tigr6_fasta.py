#!/usr/bin/env python3

from pathlib import Path
import csv
import re

BED_IN = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.bed")
FAI = Path("species/oryza_sativa_japonica/data/raw/references/tigr6/tigr6_pseudomolecules.fa.fai")

BED_OUT = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.fasta_seqids.bed")
MAP_OUT = Path("species/oryza_sativa_japonica/results/qc/yu2011_tigr6_chr_seqid_mapping.tsv")

seqids = []
with FAI.open() as f:
    for line in f:
        seqid, length, *_ = line.rstrip("\n").split("\t")
        seqids.append(seqid)

def chr_number_from_seqid(seqid):
    s = seqid.lower()
    m = re.search(r'chr0*([0-9]+)', s)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 12:
            return str(n)
    if s in {str(i) for i in range(1, 13)}:
        return s
    return None

chr_to_seqid = {}
for seqid in seqids:
    c = chr_number_from_seqid(seqid)
    if c and c not in chr_to_seqid:
        chr_to_seqid[c] = seqid

missing = [str(i) for i in range(1, 13) if str(i) not in chr_to_seqid]
if missing:
    raise SystemExit(
        "Could not map these chromosomes to FASTA seqids: "
        + ",".join(missing)
        + "\nFASTA seqids: "
        + ",".join(seqids[:20])
    )

MAP_OUT.parent.mkdir(parents=True, exist_ok=True)

with MAP_OUT.open("w", newline="") as w:
    writer = csv.writer(w, delimiter="\t")
    writer.writerow(["chr", "tigr6_seqid"])
    for c in map(str, range(1, 13)):
        writer.writerow([c, chr_to_seqid[c]])

n = 0
with BED_IN.open() as f, BED_OUT.open("w") as w:
    for line in f:
        chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
        new_chrom = chr_to_seqid[chrom]
        w.write(f"{new_chrom}\t{start}\t{end}\t{name}\n")
        n += 1

print(f"Written: {BED_OUT}")
print(f"Written: {MAP_OUT}")
print(f"Rows: {n}")

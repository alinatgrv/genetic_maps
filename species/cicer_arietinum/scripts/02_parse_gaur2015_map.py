#!/usr/bin/env python3
"""Parse Gaur et al. 2015 (Sci Rep 5:13387) Supplementary Information (PDF).

Table S1 gives SNP flanking sequences (…[A/G]… form) for 6144 SNPs.
Table S2 gives the genetic position (linkage group CaLG1-8 and cM) of 6698
markers of the interspecific (ICC4958 x PI489777) chickpea map.

We join by marker name to obtain, per SNP marker, its linkage group, cM and a
probe sequence for alignment. The two tables are distinguished by row shape
(S1 rows carry an [X/Y] SNP bracket; S2 rows start with a CaLG tag), so both
are parsed from the full document text without hard-coded page ranges.

Outputs:
    data/metadata/cicer_snp_markers_metadata.tsv  (marker_id, LG, cM, class, flank)
    data/markers/cicer_snp_probes.fasta           (probe per marker for BLAST)

Run from repo root in the genetic_maps env (requires pypdf).
"""
import re
import csv
from pathlib import Path
from collections import Counter
from pypdf import PdfReader

SPECIES_DIR = Path("species/cicer_arietinum")
PDF_IN = SPECIES_DIR / "data/raw/gaur2015_srep13387/srep13387-s1.pdf"
META_OUT = SPECIES_DIR / "data/metadata/cicer_snp_markers_metadata.tsv"
FASTA_OUT = SPECIES_DIR / "data/markers/cicer_snp_probes.fasta"

# S1: marker followed by a flanking sequence carrying the [A/G] SNP bracket
S1_RE = re.compile(r"(Ca[T]?SNP\d+)\s+([ACGTN]*\[[ACGT]/[ACGT]\][ACGTN]*)")
# S2: "<row#> CaLG<n> <marker> <cM>"
S2_RE = re.compile(r"\b\d+\s+(CaLG[1-8])\s+(\S+)\s+(\d+\.\d+)")
BRACKET_RE = re.compile(r"\[([ACGT])/[ACGT]\]")


def marker_class(mid: str) -> str:
    return "transcript_SNP" if mid.startswith("CaTSNP") else "genomic_SNP"


def main():
    reader = PdfReader(str(PDF_IN))
    text = " ".join((p.extract_text() or "").replace("\n", " ") for p in reader.pages)

    flank = dict(S1_RE.findall(text))                 # marker -> flanking seq
    gmap = {}                                         # marker -> (LG, cM)
    for lg, mk, cm in S2_RE.findall(text):
        gmap.setdefault(mk, (int(lg[-1]), float(cm)))

    joined = {m: (gmap[m][0], gmap[m][1], flank[m]) for m in gmap if m in flank}

    META_OUT.parent.mkdir(parents=True, exist_ok=True)
    FASTA_OUT.parent.mkdir(parents=True, exist_ok=True)

    n_probe = 0
    with open(META_OUT, "w", newline="") as mh, open(FASTA_OUT, "w") as fh:
        w = csv.writer(mh, delimiter="\t")
        w.writerow(["marker_id", "linkage_group", "cM", "marker_class", "flanking_sequence"])
        for mid, (lg, cm, seq) in sorted(joined.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            probe = BRACKET_RE.sub(r"\1", seq).replace("N", "")
            if len(probe) < 40:
                continue
            w.writerow([mid, lg, cm, marker_class(mid), seq])
            fh.write(f">{mid}\n{probe}\n")
            n_probe += 1

    by_lg = Counter(v[0] for v in joined.values())
    by_cls = Counter(marker_class(m) for m in joined)
    print(f"Table S1 flanking sequences: {len(flank)}")
    print(f"Table S2 map markers:        {len(gmap)}")
    print(f"JOIN (cM + flanking seq):    {len(joined)}")
    print(f"probes written:              {n_probe}")
    print(f"metadata: {META_OUT}")
    print(f"probes:   {FASTA_OUT}")
    print("\nby linkage group:", {f"LG{k}": by_lg[k] for k in sorted(by_lg)})
    print("by class:", dict(by_cls))


if __name__ == "__main__":
    main()

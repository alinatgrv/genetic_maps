#!/usr/bin/env python3
"""Build the chickpea candidate genetic map from SNP-probe BLAST hits.

For each marker keep the best BLAST hit (max bitscore) provided:
  - it covers >= 80% of the probe length,
  - it is unique (runner-up bitscore < 0.95 x best),
  - it lands on a placed Ca chromosome (Ca1-Ca8) whose number matches the
    marker's linkage group.
Physical position = midpoint of the hit; cM is taken as-is from the map.

Output (chr/pos/cM + marker metadata):
    results/intermediate/cicer_genetic_map.candidate.raw.tsv
    results/intermediate/cicer_map.candidate.summary.txt
Run from repo root in the genetic_maps env.
"""
import pandas as pd
from pathlib import Path

SPECIES_DIR = Path("species/cicer_arietinum")
BLAST = SPECIES_DIR / "results/intermediate/cicer_snp_probes.blastn.tsv"
META = SPECIES_DIR / "data/metadata/cicer_snp_markers_metadata.tsv"
OUT = SPECIES_DIR / "results/intermediate"

MIN_COV = 0.80
UNIQUE_RATIO = 0.95   # runner-up bitscore must be below this fraction of the best
CHR_MAP = {f"NC_0211{60 + i}.1": i + 1 for i in range(8)}   # Ca1..Ca8

COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]


def main():
    b = pd.read_csv(BLAST, sep="\t", names=COLS)
    b = b[b["length"] >= MIN_COV * b["qlen"]]
    meta = pd.read_csv(META, sep="\t").set_index("marker_id")

    rows, n_unique, n_onchr = [], 0, 0
    for mid, g in b.groupby("qseqid"):
        if mid not in meta.index:
            continue
        g = g.sort_values("bitscore", ascending=False)
        best = g.iloc[0]
        runner = g.iloc[1]["bitscore"] if len(g) > 1 else 0.0
        if runner >= UNIQUE_RATIO * best["bitscore"]:
            continue                          # not a unique placement
        n_unique += 1
        chrom = CHR_MAP.get(best["sseqid"], 0)
        if chrom == 0:
            continue
        n_onchr += 1
        lg = int(meta.loc[mid, "linkage_group"])
        if chrom != lg:
            continue
        rows.append({
            "chr": chrom,
            "pos": int(round((best["sstart"] + best["send"]) / 2)),
            "cM": float(meta.loc[mid, "cM"]),
            "marker_id": mid,
            "linkage_group": lg,
            "marker_class": meta.loc[mid, "marker_class"],
            "pident": round(float(best["pident"]), 2),
        })

    gmap = pd.DataFrame(rows).sort_values(["chr", "pos", "cM"])
    map_out = OUT / "cicer_genetic_map.candidate.raw.tsv"
    gmap.to_csv(map_out, sep="\t", index=False)

    lines = [
        "Chickpea candidate map build summary (SNP probes, Gaur et al. 2015)",
        "=" * 55,
        f"probes with BLAST hits: {b['qseqid'].nunique()}",
        f"unique placements:      {n_unique}",
        f"unique on a Ca chromosome: {n_onchr}",
        f"FINAL placed (unique + chr==LG): {len(gmap)}  "
        f"({100 * len(gmap) / max(1, n_onchr):.0f}% of on-chromosome)",
        "",
        "placed markers by chromosome:",
        gmap["chr"].value_counts().sort_index().to_string(),
        "",
        "placed markers by class:",
        gmap["marker_class"].value_counts().to_string(),
    ]
    (OUT / "cicer_map.candidate.summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n-> {map_out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-6}"

MINIMAP2="/mnt/users/tagirovaa/bin/minimap2"
K8="/mnt/users/grigorieval/miniconda3/bin/k8"
PAFTOOLS="/mnt/users/grigorieval/miniconda3/bin/paftools.js"

TIGR6_REF="species/oryza_sativa_japonica/data/raw/references/tigr6/tigr6_pseudomolecules.fa"
INDICA_REF="/mnt/reference/genomes/oryza_sativa_indica/GCA_000004655.2/Oryza_indica.ASM465v1.dna.toplevel.fa"

BED_IN="species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.fasta_seqids.bed"

OUTDIR="species/oryza_sativa_indica/results/liftover/yu2011_tigr6_to_asm465v1"
mkdir -p "$OUTDIR"

PAF="$OUTDIR/tigr6_to_asm465v1.asm5.paf"
LIFTED_BED="$OUTDIR/yu2011_bins_asm465v1_lifted.bed"
MINIMAP_LOG="$OUTDIR/minimap2_tigr6_to_asm465v1.log"
LIFTOVER_LOG="$OUTDIR/paftools_liftover.log"

echo "THREADS: $THREADS"
echo "TIGR6_REF: $TIGR6_REF"
echo "INDICA_REF: $INDICA_REF"
echo "BED_IN: $BED_IN"

test -s "$MINIMAP2"
test -s "$K8"
test -s "$PAFTOOLS"
test -s "$TIGR6_REF"
test -s "$INDICA_REF"
test -s "$BED_IN"

echo
echo "Running minimap2 TIGR6.1 -> ASM465v1..."
"$MINIMAP2" -t "$THREADS" -cx asm5 "$INDICA_REF" "$TIGR6_REF" > "$PAF" 2> "$MINIMAP_LOG"

echo
echo "Running paftools liftover..."
"$K8" "$PAFTOOLS" liftover -l 1 "$PAF" "$BED_IN" > "$LIFTED_BED" 2> "$LIFTOVER_LOG"

echo
echo "Line counts:"
wc -l "$BED_IN" "$LIFTED_BED" "$PAF"

echo
echo "First lifted rows:"
head "$LIFTED_BED"

echo
echo "Minimap2 log:"
cat "$MINIMAP_LOG"

echo
echo "LiftOver log:"
cat "$LIFTOVER_LOG" || true

#!/usr/bin/env bash
set -euo pipefail

THREADS="${SLURM_NTASKS:-6}"

MINIMAP2="/mnt/users/tagirovaa/bin/minimap2"
K8="/mnt/users/grigorieval/miniconda3/bin/k8"
PAFTOOLS="/mnt/users/grigorieval/miniconda3/bin/paftools.js"

TIGR6_REF="species/oryza_sativa_japonica/data/raw/references/tigr6/tigr6_pseudomolecules.fa"
IRGSP1_REF="/mnt/reference/genomes/oryza_sativa_japonica/GCF_001433935.1/GCF_001433935.1_IRGSP-1.0_genomic.unmasked.fna"

BED_IN="species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6_midpoints.fasta_seqids.bed"

OUTDIR="species/oryza_sativa_japonica/results/liftover/yu2011_tigr6_to_irgsp1"
mkdir -p "$OUTDIR"

PAF="$OUTDIR/tigr6_to_irgsp1.asm5.paf"
LIFTED_BED="$OUTDIR/yu2011_bins_irgsp1_lifted.bed"
LIFTOVER_LOG="$OUTDIR/paftools_liftover.log"
MINIMAP_LOG="$OUTDIR/minimap2_tigr6_to_irgsp1.log"

echo "THREADS: $THREADS"
echo "MINIMAP2: $MINIMAP2"
echo "K8: $K8"
echo "PAFTOOLS: $PAFTOOLS"
echo "TIGR6_REF: $TIGR6_REF"
echo "IRGSP1_REF: $IRGSP1_REF"
echo "BED_IN: $BED_IN"
echo "OUTDIR: $OUTDIR"

test -s "$MINIMAP2"
test -s "$K8"
test -s "$PAFTOOLS"
test -s "$TIGR6_REF"
test -s "$IRGSP1_REF"
test -s "$BED_IN"

echo
echo "Checking paftools liftover help..."
"$K8" "$PAFTOOLS" liftover 2>&1 | head -40 || true

echo
echo "Running minimap2 assembly-to-assembly alignment..."
"$MINIMAP2" -t "$THREADS" -cx asm5 "$IRGSP1_REF" "$TIGR6_REF" > "$PAF" 2> "$MINIMAP_LOG"

echo
echo "PAF written:"
ls -lh "$PAF"

echo
echo "Running paftools liftover..."
"$K8" "$PAFTOOLS" liftover -l 1 "$PAF" "$BED_IN" > "$LIFTED_BED" 2> "$LIFTOVER_LOG"

echo
echo "Lifted BED written:"
ls -lh "$LIFTED_BED"

echo
echo "Line counts:"
wc -l "$BED_IN" "$LIFTED_BED" "$PAF"

echo
echo "First lifted rows:"
head "$LIFTED_BED"

echo
echo "Minimap2 log:"
cat "$MINIMAP_LOG" || true

echo
echo "LiftOver log:"
cat "$LIFTOVER_LOG" || true

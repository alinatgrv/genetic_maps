#!/usr/bin/env bash

set -euo pipefail

THREADS="${THREADS:-8}"

REF="/mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.unmasked.fna"

SPECIES_DIR="species/helianthus_annuus"
QUERY="$SPECIES_DIR/data/markers/sunflower_features_25bp.fasta"

SAI="$SPECIES_DIR/results/intermediate/sunflower_features_25bp.bwa_exact.sai"
SAM="$SPECIES_DIR/results/intermediate/sunflower_features_25bp.bwa_exact.sam"
LOG="$SPECIES_DIR/logs/02_align_sunflower_features_bwa_exact.log"

mkdir -p "$SPECIES_DIR/results/intermediate" "$SPECIES_DIR/logs"

{
  echo "Start: $(date)"
  echo "Host: $(hostname)"
  echo "Threads: $THREADS"
  echo "Reference: $REF"
  echo "Query: $QUERY"

  echo
  echo "Running bwa aln exact..."
  bwa aln \
    -n 0 \
    -o 0 \
    -l 25 \
    -t "$THREADS" \
    "$REF" \
    "$QUERY" \
    > "$SAI"

  echo
  echo "Running bwa samse..."
  bwa samse \
    "$REF" \
    "$SAI" \
    "$QUERY" \
    > "$SAM"

  echo
  echo "Done: $(date)"
  echo "SAI:"
  ls -lh "$SAI"
  echo "SAM:"
  ls -lh "$SAM"
  echo "SAM alignment rows:"
  grep -v '^@' "$SAM" | wc -l
} 2>&1 | tee "$LOG"

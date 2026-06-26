#!/usr/bin/env bash
set -euo pipefail

# Выравнивание последовательностей DArT-клонов на целевой референс
# GCF_020379485.1 (Da-Ae) с помощью BLASTN (megablast).
#
# DArT-клоны — это полноразмерные фрагменты (медиана ~420 bp), поэтому megablast
# подходит лучше, чем bwa/short. Физическую позицию маркера определим по лучшему
# хиту на шаге 05.
#
# Запускать из корня репозитория (в окружении genetic_maps).

THREADS="${THREADS:-8}"

SPECIES_DIR="species/brassica_napus"
REF="${SPECIES_DIR}/data/ref/GCF_020379485.1_genomic.fna"
DB="${SPECIES_DIR}/data/ref/GCF_020379485.1_blastdb"
QUERY="${SPECIES_DIR}/data/markers/bnapus_dart_markers.fasta"

HITS="${SPECIES_DIR}/results/intermediate/bnapus_dart_markers.blastn.tsv"
LOG="${SPECIES_DIR}/logs/04_align_bnapus_dart_markers.log"

mkdir -p "${SPECIES_DIR}/results/intermediate" "${SPECIES_DIR}/logs"

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

{
  echo "Start: $(date)"
  echo "Threads: ${THREADS}"
  echo "Reference: ${REF}"
  echo "Query: ${QUERY}"

  if [[ ! -f "${DB}.nsq" && ! -f "${DB}.nal" ]]; then
    echo
    echo "Building BLAST nucleotide database..."
    makeblastdb -in "${REF}" -dbtype nucl -out "${DB}"
  else
    echo "BLAST database already present: ${DB}"
  fi

  echo
  echo "Running blastn (megablast)..."
  blastn \
    -task megablast \
    -query "${QUERY}" \
    -db "${DB}" \
    -evalue 1e-10 \
    -perc_identity 90 \
    -max_target_seqs 20 \
    -num_threads "${THREADS}" \
    -outfmt "${OUTFMT}" \
    > "${HITS}"

  echo
  echo "Done: $(date)"
  echo "Hits file: ${HITS}"
  echo "Total HSP rows: $(wc -l < "${HITS}")"
  echo "Queries with at least one hit: $(cut -f1 "${HITS}" | sort -u | wc -l)"
} 2>&1 | tee "${LOG}"

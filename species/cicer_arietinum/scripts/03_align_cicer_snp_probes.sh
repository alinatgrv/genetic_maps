#!/usr/bin/env bash
set -euo pipefail

# Выравнивание SNP-зондов (фланки Gaur et al. 2015, ~120 bp) на референс нута
# GCF_000331145.1 (ASM33114v1) с помощью BLASTN (megablast).
#
# Фланки длинные и специфичные, поэтому megablast даёт по сути один уверенный хит
# на маркер. Физическую позицию (середину хита) и уникальность определяет шаг 04.
#
# Запускать из корня репозитория (в окружении genetic_maps).

THREADS="${THREADS:-8}"

SPECIES_DIR="species/cicer_arietinum"
REF="${SPECIES_DIR}/data/ref/GCF_000331145.1_genomic.fna"
DB="${SPECIES_DIR}/data/ref/GCF_000331145.1_blastdb"
QUERY="${SPECIES_DIR}/data/markers/cicer_snp_probes.fasta"

HITS="${SPECIES_DIR}/results/intermediate/cicer_snp_probes.blastn.tsv"
LOG="${SPECIES_DIR}/logs/03_align_cicer_snp_probes.log"

mkdir -p "${SPECIES_DIR}/results/intermediate" "${SPECIES_DIR}/logs"

OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

{
  echo "Start: $(date)"
  echo "Reference: ${REF}"
  echo "Query: ${QUERY}"

  if [[ ! -f "${DB}.nsq" && ! -f "${DB}.nal" ]]; then
    echo "Building BLAST nucleotide database..."
    makeblastdb -in "${REF}" -dbtype nucl -out "${DB}"
  else
    echo "BLAST database already present: ${DB}"
  fi

  echo "Running blastn (megablast)..."
  blastn \
    -task megablast \
    -query "${QUERY}" \
    -db "${DB}" \
    -evalue 1e-15 \
    -perc_identity 90 \
    -max_target_seqs 5 \
    -num_threads "${THREADS}" \
    -outfmt "${OUTFMT}" \
    > "${HITS}"

  echo "Done: $(date)"
  echo "Hits file: ${HITS}"
  echo "Total HSP rows: $(wc -l < "${HITS}")"
  echo "Probes with at least one hit: $(cut -f1 "${HITS}" | sort -u | wc -l)"
} 2>&1 | tee "${LOG}"

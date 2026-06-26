#!/usr/bin/env bash
set -euo pipefail

# Перенос bin-серединок Yu 2011 с IRGSP-1.0 (Nipponbare) на indica-родителя ZS97RS3
# через короткие фланкирующие окна (быстрый локальный liftover).
#
# Полногеномное base-level выравнивание двух разошедшихся рисовых геномов в режиме
# asm5 оказалось патологически медленным (десятки минут, без прогресса). Эквивалентный
# по результату, но быстрый путь: вокруг каждой точки берём окно ±500 bp на IRGSP-1.0,
# выравниваем окна на ZS97RS3 (minimap2) и переносим центр окна через paftools.js liftover.
#
# Запускать из корня репозитория в env genetic_maps (minimap2, samtools, k8, paftools.js).

THREADS="${THREADS:-6}"

SPECIES_DIR="species/oryza_sativa_indica"
TARGET="${SPECIES_DIR}/data/ref/GCA_001623345.3_genomic.fna"   # ZS97RS3
QUERY_REF="${SPECIES_DIR}/data/ref/GCF_001433935.1_genomic.fna" # IRGSP-1.0 (источник окон)

OUT_DIR="${SPECIES_DIR}/results/liftover/yu2011_irgsp1_to_zs97rs3"
WINDOWS="${OUT_DIR}/windows_irgsp1_flank.fasta"
CENTER_BED="${OUT_DIR}/windows_center.bed"
PAF="${OUT_DIR}/windows_irgsp1_to_zs97rs3.paf"
LIFTED="${OUT_DIR}/yu2011_bins_zs97rs3_lifted.bed"
LOG="${SPECIES_DIR}/logs/09_liftover_irgsp1_to_zs97rs3.log"

mkdir -p "${OUT_DIR}" "${SPECIES_DIR}/logs"

{
  echo "Start: $(date)"

  echo
  echo "Indexing IRGSP-1.0 (faidx, if needed)..."
  [[ -f "${QUERY_REF}.fai" ]] || samtools faidx "${QUERY_REF}"

  echo "Extracting flanking windows around bin midpoints..."
  python3 "${SPECIES_DIR}/scripts/09_make_flank_windows.py"

  echo
  echo "Aligning windows to ZS97RS3 (minimap2 asm5, primary only)..."
  minimap2 -t "${THREADS}" -cx asm5 --secondary=no "${TARGET}" "${WINDOWS}" > "${PAF}"
  echo "PAF lines: $(wc -l < "${PAF}")"

  echo
  echo "Lifting window centers with paftools.js liftover..."
  paftools.js liftover -l 1 "${PAF}" "${CENTER_BED}" > "${LIFTED}"

  echo
  echo "Done: $(date)"
  echo "Lifted BED: ${LIFTED}"
  echo "Input bins: $(wc -l < "${CENTER_BED}") ; lifted rows: $(wc -l < "${LIFTED}")"
} 2>&1 | tee "${LOG}"

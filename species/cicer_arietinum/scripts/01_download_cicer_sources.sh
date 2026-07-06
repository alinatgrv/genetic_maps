#!/usr/bin/env bash
set -euo pipefail

# Загрузка исходных данных для построения генетической карты нута (Cicer arietinum).
#
# Два источника:
#   1. Целевой референс GCF_000331145.1 (ASM33114v1, CDC Frontier, kabuli),
#      8 псевдохромосом Ca1-Ca8 (NC_021160.1..NC_021167.1)      -> data/ref/ (gitignored)
#   2. Высокоплотная SNP-карта, Gaur et al. 2015, Sci Rep 5:13387,
#      Supplementary Information (srep13387-s1.pdf):
#        Table S1 — фланкирующие последовательности 6144 SNP,
#        Table S2 — cM/LG для 6698 маркеров межвидовой карты
#                                                              -> data/raw/gaur2015_srep13387/
#
# Запускать из корня репозитория (в окружении genetic_maps).
# Тяжёлую загрузку референса можно отключить через DOWNLOAD_REF=0.

SPECIES_DIR="species/cicer_arietinum"
REF_DIR="${SPECIES_DIR}/data/ref"
RAW_DIR="${SPECIES_DIR}/data/raw/gaur2015_srep13387"

mkdir -p "${REF_DIR}" "${RAW_DIR}"

DOWNLOAD_SUPP="${DOWNLOAD_SUPP:-1}"
DOWNLOAD_REF="${DOWNLOAD_REF:-1}"

# --- 2. Supplementary Information (Gaur et al. 2015) ---
if [[ "${DOWNLOAD_SUPP}" == "1" ]]; then
  echo "Downloading Gaur et al. 2015 Supplementary Information (srep13387-s1.pdf)..."
  curl -fL --retry 3 \
    "https://static-content.springer.com/esm/art%3A10.1038%2Fsrep13387/MediaObjects/41598_2015_BFsrep13387_MOESM1_ESM.pdf" \
    -o "${RAW_DIR}/srep13387-s1.pdf"
fi

# --- 1. Целевой референс GCF_000331145.1 (ASM33114v1) ---
if [[ "${DOWNLOAD_REF}" == "1" ]]; then
  echo "Downloading reference genome GCF_000331145.1 (ASM33114v1) via NCBI datasets..."
  TMP_ZIP="${REF_DIR}/GCF_000331145.1.zip"
  datasets download genome accession GCF_000331145.1 \
    --include genome \
    --filename "${TMP_ZIP}"
  ( cd "${REF_DIR}" && unzip -o "GCF_000331145.1.zip" )
  FNA=$(find "${REF_DIR}/ncbi_dataset" -name "*_genomic.fna" | head -1)
  cp "${FNA}" "${REF_DIR}/GCF_000331145.1_genomic.fna"
  echo "Reference FASTA: ${REF_DIR}/GCF_000331145.1_genomic.fna"
fi

echo
echo "Downloaded files:"
ls -lh "${RAW_DIR}" 2>/dev/null || true
ls -lh "${REF_DIR}"/*.fna 2>/dev/null || true

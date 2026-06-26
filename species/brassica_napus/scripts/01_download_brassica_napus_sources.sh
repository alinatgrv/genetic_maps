#!/usr/bin/env bash
set -euo pipefail

# Загрузка исходных данных для построения генетической карты рапса (Brassica napus).
#
# Три источника:
#   1. Целевой референс GCF_020379485.1 (Da-Ae), хромосомный уровень  -> data/ref/ (gitignored)
#   2. Консенсус-карта DArT, Raman et al. 2014, BMC Genomics 14:277,
#      Additional file 1 (1471-2164-14-277-S1.xls)                    -> data/raw/raman2014_dart/
#   3. Последовательности DArT-маркеров рода Brassica (DArT_Brassica.fasta,
#      Diversity Arrays Technology, бесплатно по T&C; не перераспространять) -> data/markers/ (gitignored)
#
# Запускать из корня репозитория. Тяжёлые загрузки (референс, FASTA) можно
# включать выборочно через переменные окружения DOWNLOAD_REF / DOWNLOAD_FASTA / DOWNLOAD_MAP.

SPECIES_DIR="species/brassica_napus"
REF_DIR="${SPECIES_DIR}/data/ref"
RAMAN_DIR="${SPECIES_DIR}/data/raw/raman2014_dart"
MARKERS_DIR="${SPECIES_DIR}/data/markers"

mkdir -p "${REF_DIR}" "${RAMAN_DIR}" "${MARKERS_DIR}"

DOWNLOAD_MAP="${DOWNLOAD_MAP:-1}"
DOWNLOAD_FASTA="${DOWNLOAD_FASTA:-1}"
DOWNLOAD_REF="${DOWNLOAD_REF:-1}"

# --- 2. Консенсус-карта Raman et al. 2014 (Additional file 1) ---
if [[ "${DOWNLOAD_MAP}" == "1" ]]; then
  echo "Downloading Raman et al. 2014 consensus map (Additional file 1)..."
  # Прямой CDN Springer (PMC отдаёт анти-бот HTML-заглушку вместо файла).
  curl -fL --retry 3 \
    "https://static-content.springer.com/esm/art%3A10.1186%2F1471-2164-14-277/MediaObjects/12864_2012_4970_MOESM1_ESM.xls" \
    -o "${RAMAN_DIR}/1471-2164-14-277-S1.xls"
fi

# --- 3. Последовательности DArT-маркеров (Brassica) ---
if [[ "${DOWNLOAD_FASTA}" == "1" ]]; then
  echo "Downloading DArT Brassica marker sequences (Diversity Arrays Technology)..."
  curl -fL --retry 3 \
    "https://www.diversityarrays.com/files/fasta/DArT_Brassica.fasta" \
    -o "${MARKERS_DIR}/DArT_Brassica.fasta"
fi

# --- 1. Целевой референс GCF_020379485.1 (Da-Ae) ---
if [[ "${DOWNLOAD_REF}" == "1" ]]; then
  echo "Downloading reference genome GCF_020379485.1 (Da-Ae) via NCBI datasets..."
  TMP_ZIP="${REF_DIR}/GCF_020379485.1.zip"
  datasets download genome accession GCF_020379485.1 \
    --include genome \
    --filename "${TMP_ZIP}"
  ( cd "${REF_DIR}" && unzip -o "GCF_020379485.1.zip" )
  # Кладём геном на удобный путь.
  FNA=$(find "${REF_DIR}/ncbi_dataset" -name "*_genomic.fna" | head -1)
  cp "${FNA}" "${REF_DIR}/GCF_020379485.1_genomic.fna"
  echo "Reference FASTA: ${REF_DIR}/GCF_020379485.1_genomic.fna"
fi

echo
echo "Downloaded files:"
ls -lh "${RAMAN_DIR}" "${MARKERS_DIR}" 2>/dev/null || true

echo
echo "File types:"
file "${RAMAN_DIR}"/* "${MARKERS_DIR}"/*.fasta 2>/dev/null || true

#!/usr/bin/env bash
set -euo pipefail

# Загрузка референсов для пересборки indica-карты на родительском геноме кросса Yu 2011.
#
# Цель (target): ZS97RS3 — Zhenshan 97, реальный родитель RIL-популяции Yu et al. 2011,
#   GCA_001623345.3 (Complete Genome, без гэпов).
# Источник (query) для liftover: IRGSP-1.0 — современная сборка Nipponbare,
#   GCF_001433935.1 (на ней уже лежат bin-позиции из japonica-карты).
#
# Качаем напрямую с NCBI FTP (curl, резюмируемо). Оба генома gitignored (data/ref/).
# Запускать из корня репозитория.

SPECIES_DIR="species/oryza_sativa_indica"
REF_DIR="${SPECIES_DIR}/data/ref"
mkdir -p "${REF_DIR}"

ZS97_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/623/345/GCA_001623345.3_ZS97RS3/GCA_001623345.3_ZS97RS3_genomic.fna.gz"
IRGSP_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/433/935/GCF_001433935.1_IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic.fna.gz"

fetch () {
  local url="$1" out="$2" label="$3"
  echo "Downloading ${label} ..."
  curl -fL -C - --retry 5 --retry-delay 3 "${url}" -o "${out}"
}

# Последовательная загрузка с резюмированием (-C -): надёжнее параллельной на
# нестабильном соединении к NCBI (одно активное соединение за раз).
fetch "${ZS97_URL}"  "${REF_DIR}/GCA_001623345.3_genomic.fna.gz" "ZS97RS3 (Zhenshan 97, indica target)"
fetch "${IRGSP_URL}" "${REF_DIR}/GCF_001433935.1_genomic.fna.gz" "IRGSP-1.0 (Nipponbare, query)"

echo "Decompressing..."
gunzip -kf "${REF_DIR}/GCA_001623345.3_genomic.fna.gz"
gunzip -kf "${REF_DIR}/GCF_001433935.1_genomic.fna.gz"

echo
echo "Reference FASTAs:"
ls -lh "${REF_DIR}"/*.fna

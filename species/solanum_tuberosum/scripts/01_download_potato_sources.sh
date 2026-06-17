#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="species/solanum_tuberosum/data/raw"
SPUDDB_DIR="${RAW_DIR}/spuddb_dm_v6_1"
PGSC_DIR="${RAW_DIR}/pgsc_sharma_2013"

mkdir -p "${SPUDDB_DIR}" "${PGSC_DIR}"

BASE_URL="https://spuddb.uga.edu/data/dm_v61"

echo "Removing previous failed HTML downloads..."
rm -f "${SPUDDB_DIR}/solcap_69k_SNPs_DM_v6_1_pos.txt.gz"
rm -f "${SPUDDB_DIR}/potvar_SNPs_DM_v6_1_pos.txt.gz"
rm -f "${SPUDDB_DIR}/solcap_69k_potvar_SNP_pos_DM_v6_1.xlsx"

echo "Downloading potato SNP physical positions on DM v6.1 from SpudDB..."

curl -L \
  "${BASE_URL}/solcap_69k_SNPs_DM_v6_1_pos.txt.gz" \
  -o "${SPUDDB_DIR}/solcap_69k_SNPs_DM_v6_1_pos.txt.gz"

curl -L \
  "${BASE_URL}/potvar_SNPs_DM_v6_1_pos.txt.gz" \
  -o "${SPUDDB_DIR}/potvar_SNPs_DM_v6_1_pos.txt.gz"

curl -L \
  "${BASE_URL}/solcap_69k_potvar_SNP_pos_DM_v6_1.xlsx" \
  -o "${SPUDDB_DIR}/solcap_69k_potvar_SNP_pos_DM_v6_1.xlsx"

echo
echo "Downloaded files:"
ls -lh "${SPUDDB_DIR}"

echo
echo "File types:"
file "${SPUDDB_DIR}"/*

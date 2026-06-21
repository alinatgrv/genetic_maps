#!/usr/bin/env bash
set -euo pipefail

BASE="species/solanum_tuberosum"
REF_DIR="${BASE}/data/ref/liftover_v403_to_v61"
mkdir -p "${REF_DIR}"

# Старый референс, на котором лежат координаты маркеров из Sharma et al. 2013 TableS4.
OLD_URL="https://spuddb.uga.edu/data/pgsc_legacy/PGSC_DM_v4.03_pseudomolecules.fasta.zip"
OLD_ZIP="${REF_DIR}/PGSC_DM_v4.03_pseudomolecules.fasta.zip"

# Новый целевой референс уже лежит на сервере.
NEW_FA="/mnt/reference/genomes/solanum_tuberosum/DM_1-3_516_R44_potato.v6.1.fa"

echo "Checking target DM v6.1 reference..."
if [ ! -s "${NEW_FA}" ]; then
  echo "ERROR: target reference not found or empty:"
  echo "${NEW_FA}"
  exit 1
fi

echo "Target reference found:"
ls -lh "${NEW_FA}"

echo
echo "Downloading old PGSC/DM v4.03 pseudomolecules..."
if [ ! -s "${OLD_ZIP}" ]; then
  curl -L "${OLD_URL}" -o "${OLD_ZIP}"
else
  echo "Already exists: ${OLD_ZIP}"
fi

echo
echo "Checking downloaded old reference archive..."
file "${OLD_ZIP}"
ls -lh "${OLD_ZIP}"

echo
echo "Unpacking old v4.03 zip..."
rm -rf "${REF_DIR}/v403_unzipped"
mkdir -p "${REF_DIR}/v403_unzipped"
unzip -o "${OLD_ZIP}" -d "${REF_DIR}/v403_unzipped"

echo
echo "Finding old v4.03 FASTA..."
OLD_FA=$(find "${REF_DIR}/v403_unzipped" -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | head -n 1)

if [ -z "${OLD_FA}" ]; then
  echo "ERROR: old v4.03 FASTA not found after unzip"
  find "${REF_DIR}/v403_unzipped" -type f | head -n 50
  exit 1
fi

cat > "${REF_DIR}/liftover_reference_paths.env" <<PATHS
OLD_FA=${OLD_FA}
NEW_FA=${NEW_FA}
PATHS

echo
echo "Reference paths:"
cat "${REF_DIR}/liftover_reference_paths.env"

echo
echo "FASTA headers:"
echo "old v4.03:"
grep '^>' "${OLD_FA}" | head
echo
echo "new v6.1:"
grep '^>' "${NEW_FA}" | head

echo
echo "Done."

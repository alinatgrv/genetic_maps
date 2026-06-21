#!/usr/bin/env bash
set -euo pipefail

BASE="species/solanum_tuberosum"
REF_DIR="${BASE}/data/ref/liftover_v403_to_v61"
OUT_DIR="${BASE}/results/intermediate/liftover_v403_to_v61"
LOG_DIR="${BASE}/logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

export PATH="$HOME/bin:$PATH"

source "${REF_DIR}/liftover_reference_paths.env"

THREADS="${THREADS:-6}"

PAF="${OUT_DIR}/potato_v403_to_v61.minimap2.asm5.paf"
LOG="${LOG_DIR}/07_align_potato_v403_to_v61_minimap2.log"

echo "minimap2: $(command -v minimap2)"
minimap2 --version

echo "Old query genome: ${OLD_FA}"
echo "New target genome: ${NEW_FA}"
echo "Threads: ${THREADS}"
echo "Output PAF: ${PAF}"
echo "Log: ${LOG}"

if [ ! -s "${OLD_FA}" ]; then
  echo "ERROR: old v4.03 FASTA not found: ${OLD_FA}"
  exit 1
fi

if [ ! -s "${NEW_FA}" ]; then
  echo "ERROR: new v6.1 FASTA not found: ${NEW_FA}"
  exit 1
fi

/usr/bin/time -p minimap2 \
  -x asm5 \
  -t "${THREADS}" \
  -c \
  --cs=long \
  --secondary=no \
  "${NEW_FA}" \
  "${OLD_FA}" \
  > "${PAF}" \
  2> "${LOG}"

echo "Done."
ls -lh "${PAF}"
tail -n 30 "${LOG}"
head "${PAF}"
wc -l "${PAF}"

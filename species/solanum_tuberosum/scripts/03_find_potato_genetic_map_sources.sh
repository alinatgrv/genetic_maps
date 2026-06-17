#!/usr/bin/env bash
set -euo pipefail

BASE="species/solanum_tuberosum"
RAW_DIR="${BASE}/data/raw/pgsc_sharma_2013"
mkdir -p "${RAW_DIR}"

ARTICLE_URL="https://academic.oup.com/g3journal/article/3/11/2031/6025374"

echo "Downloading article page..."
curl -L "${ARTICLE_URL}" -o "${RAW_DIR}/sharma_2013_article_page.html"

echo
echo "Searching for supplementary/data/xlsx/zip links..."
grep -Eio 'https?://[^"]+|href="[^"]+"' "${RAW_DIR}/sharma_2013_article_page.html" \
  | sed 's/^href="//; s/"$//' \
  | grep -Ei 'supp|suppl|xlsx|xls|zip|csv|tsv|txt|doc|data|DC1|g3\.113\.007153' \
  | sort -u \
  > "${RAW_DIR}/candidate_supplementary_links.txt" || true

echo
echo "Candidate links:"
cat "${RAW_DIR}/candidate_supplementary_links.txt"

echo
echo "Saved:"
ls -lh "${RAW_DIR}"

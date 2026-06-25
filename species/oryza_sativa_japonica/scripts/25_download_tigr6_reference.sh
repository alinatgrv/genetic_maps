#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="species/oryza_sativa_japonica/data/raw/references/tigr6"
mkdir -p "$OUT_DIR"

OUT="$OUT_DIR/tigr6_pseudomolecules.fa"

urls=(
"http://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.con"
"ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.con"
"https://rice.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_6.1/all.dir/all.con"
)

for url in "${urls[@]}"; do
    echo
    echo "Trying: $url"
    tmp="${OUT}.tmp"

    if wget -O "$tmp" "$url"; then
        echo "Downloaded candidate: $tmp"
        file "$tmp"
        head -5 "$tmp" || true

        if grep -q "^>" "$tmp"; then
            mv "$tmp" "$OUT"
            echo "Saved FASTA: $OUT"
            exit 0
        else
            echo "Downloaded file does not look like FASTA, removing."
            rm -f "$tmp"
        fi
    else
        echo "Failed: $url"
        rm -f "$tmp"
    fi
done

echo "ERROR: could not download TIGR6.1 pseudomolecules FASTA from tested URLs."
exit 1

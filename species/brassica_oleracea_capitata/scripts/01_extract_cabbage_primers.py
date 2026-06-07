import pandas as pd
from pathlib import Path
import re

input_file = Path("data/raw/cabbage/cabbage_markers.xls")

metadata_out = Path("data/metadata/cabbage/cabbage_markers_metadata.tsv")
fasta_out = Path("data/markers/cabbage/cabbage_primers.fasta")

metadata_out.parent.mkdir(parents=True, exist_ok=True)
fasta_out.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_excel(input_file, sheet_name="Sheet1")

# Rename columns to safer technical names
df = df.rename(columns={
    "Marker Name": "marker_id",
    "Marker Type": "marker_type",
    "Forward Primer": "forward_primer",
    "Reverse Primer": "reverse_primer",
    "LG": "linkage_group",
    "Position": "position",
    "cM": "cM"
})

required = [
    "linkage_group",
    "position",
    "cM",
    "marker_id",
    "marker_type",
    "forward_primer",
    "reverse_primer"
]

missing = [col for col in required if col not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

df = df[required].copy()

# Clean text fields
for col in ["linkage_group", "marker_id", "marker_type", "forward_primer", "reverse_primer"]:
    df[col] = df[col].astype(str).str.strip()

df["forward_primer"] = df["forward_primer"].str.upper()
df["reverse_primer"] = df["reverse_primer"].str.upper()

# Check primer characters
valid_re = re.compile(r"^[ACGTN]+$")

bad = df[
    ~df["forward_primer"].apply(lambda x: bool(valid_re.match(x))) |
    ~df["reverse_primer"].apply(lambda x: bool(valid_re.match(x)))
]

if len(bad) > 0:
    print("WARNING: primers with non-ACGTN characters found:")
    print(bad[["marker_id", "forward_primer", "reverse_primer"]].to_string(index=False))

# Save metadata
df.to_csv(metadata_out, sep="\t", index=False)

# Save FASTA
with open(fasta_out, "w") as out:
    for _, row in df.iterrows():
        marker_id = row["marker_id"]
        fwd = row["forward_primer"]
        rev = row["reverse_primer"]

        out.write(f">{marker_id}__F\n{fwd}\n")
        out.write(f">{marker_id}__R\n{rev}\n")

print(f"Input markers: {len(df)}")
print(f"Primer FASTA records: {len(df) * 2}")
print(f"Metadata written to: {metadata_out}")
print(f"Primer FASTA written to: {fasta_out}")

print("\nMarkers by linkage group:")
print(df["linkage_group"].value_counts().sort_index().to_string())

print("\nMarkers by type:")
print(df["marker_type"].value_counts().to_string())

print("\nPrimer length summary:")
primer_lengths = pd.concat([
    df["forward_primer"].str.len(),
    df["reverse_primer"].str.len()
])
print(primer_lengths.describe().to_string())

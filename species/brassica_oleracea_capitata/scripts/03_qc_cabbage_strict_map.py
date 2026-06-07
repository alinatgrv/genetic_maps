import pandas as pd
from pathlib import Path

map_file = Path("results/cabbage/cabbage_genetic_map.strict.tsv")

df = pd.read_csv(map_file, sep="\t")

print("Strict map rows:", len(df))
print("Unique markers:", df["marker_id"].nunique())
print("Unique ref_seq:", df["chr"].nunique())
print("Unique linkage groups:", df["linkage_group"].nunique())

print("\nMarkers per ref_seq:")
print(df["chr"].value_counts().to_string())

print("\nMarkers per linkage group:")
print(df["linkage_group"].value_counts().sort_index().to_string())

print("\nCross-tab: ref_seq x linkage_group")
ct = pd.crosstab(df["chr"], df["linkage_group"])
print(ct.to_string())

# Save cross-tab
out = Path("results/cabbage/cabbage_strict_refseq_vs_LG.tsv")
ct.to_csv(out, sep="\t")
print(f"\nCross-tab saved to: {out}")

# Check duplicate physical coordinates
dup = df[df.duplicated(subset=["chr", "pos"], keep=False)].copy()
dup = dup.sort_values(["chr", "pos", "marker_id"])

dup_out = Path("results/cabbage/cabbage_strict_duplicate_positions.tsv")
dup.to_csv(dup_out, sep="\t", index=False)

print("\nDuplicate physical positions:", len(dup))
if len(dup) > 0:
    print(dup[["chr", "pos", "cM", "marker_id", "linkage_group", "marker_type"]].head(30).to_string(index=False))
    print(f"Duplicate positions saved to: {dup_out}")

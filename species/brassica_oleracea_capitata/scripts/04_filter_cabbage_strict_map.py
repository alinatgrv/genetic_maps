import pandas as pd
from pathlib import Path

input_file = Path("results/cabbage/cabbage_genetic_map.strict.tsv")

clean_with_markers_out = Path("results/cabbage/cabbage_genetic_map.strict.clean_with_markers.tsv")
final_unique_out = Path("results/cabbage/cabbage_genetic_map.strict.clean.unique.tsv")
excluded_scaffolds_out = Path("results/cabbage/cabbage_strict_excluded_scaffolds.tsv")
excluded_conflicts_out = Path("results/cabbage/cabbage_strict_excluded_LG_refseq_conflicts.tsv")
duplicates_out = Path("results/cabbage/cabbage_strict_clean_duplicate_positions.tsv")
summary_out = Path("results/cabbage/cabbage_strict_clean_filtering.summary.txt")

# Mapping inferred from ref_seq x linkage_group QC
ref_to_chr = {
    "gb|CM031006.1|": "1",
    "gb|CM031007.1|": "2",
    "gb|CM031008.1|": "3",
    "gb|CM031009.1|": "4",
    "gb|CM031010.1|": "5",
    "gb|CM031011.1|": "6",
    "gb|CM031012.1|": "7",
    "gb|CM031013.1|": "8",
    "gb|CM031014.1|": "9",
}

ref_to_lg = {
    "gb|CM031006.1|": "C01",
    "gb|CM031007.1|": "C02",
    "gb|CM031008.1|": "C03",
    "gb|CM031009.1|": "C04",
    "gb|CM031010.1|": "C05",
    "gb|CM031011.1|": "C06",
    "gb|CM031012.1|": "C07",
    "gb|CM031013.1|": "C08",
    "gb|CM031014.1|": "C09",
}

df = pd.read_csv(input_file, sep="\t")

# Keep original accession before renaming
df = df.rename(columns={"chr": "ref_seq"})
df["expected_lg"] = df["ref_seq"].map(ref_to_lg)
df["chr"] = df["ref_seq"].map(ref_to_chr)

# 1. Exclude non-chromosome scaffolds/contigs
scaffolds = df[df["chr"].isna()].copy()

# 2. Among chromosome-level hits, exclude LG/ref_seq conflicts
on_chromosomes = df[df["chr"].notna()].copy()
conflicts = on_chromosomes[on_chromosomes["linkage_group"] != on_chromosomes["expected_lg"]].copy()

# 3. Clean = chromosome-level and LG agrees with expected chromosome
clean = on_chromosomes[on_chromosomes["linkage_group"] == on_chromosomes["expected_lg"]].copy()

# Sort by chromosome and physical position
clean["chr_order"] = clean["chr"].astype(int)
clean = clean.sort_values(["chr_order", "pos", "cM", "marker_id"])

# 4. Detect duplicate physical positions after cleaning
duplicates = clean[clean.duplicated(subset=["chr", "pos"], keep=False)].copy()
duplicates = duplicates.sort_values(["chr_order", "pos", "cM", "marker_id"])

# 5. Strictest final table: remove all duplicated physical positions
unique_positions = clean[~clean.duplicated(subset=["chr", "pos"], keep=False)].copy()

# Final output requested by task: chr, pos, cM
final_unique = unique_positions[["chr", "pos", "cM"]].copy()
final_unique = final_unique.sort_values(["chr", "pos"])

# Save detailed clean table with marker information
clean_with_markers = clean[[
    "chr",
    "pos",
    "cM",
    "marker_id",
    "linkage_group",
    "marker_type",
    "ref_seq",
    "amp_start",
    "amp_end",
    "amp_len"
]].copy()

clean_with_markers.to_csv(clean_with_markers_out, sep="\t", index=False)
final_unique.to_csv(final_unique_out, sep="\t", index=False)
scaffolds.to_csv(excluded_scaffolds_out, sep="\t", index=False)
conflicts.to_csv(excluded_conflicts_out, sep="\t", index=False)
duplicates.to_csv(duplicates_out, sep="\t", index=False)

with open(summary_out, "w") as out:
    out.write("Cabbage strict map filtering summary\n")
    out.write("===================================\n\n")
    out.write(f"Input strict markers: {len(df)}\n")
    out.write(f"Markers on chromosome-level CM accessions: {len(on_chromosomes)}\n")
    out.write(f"Markers on scaffold/contig accessions excluded: {len(scaffolds)}\n")
    out.write(f"LG/ref_seq conflict markers excluded: {len(conflicts)}\n")
    out.write(f"Clean markers after chromosome and LG filtering: {len(clean)}\n")
    out.write(f"Duplicate-position rows in clean map: {len(duplicates)}\n")
    out.write(f"Final unique chr-pos-cM rows: {len(final_unique)}\n\n")

    out.write("Clean markers by chromosome:\n")
    out.write(clean["chr"].value_counts().sort_index(key=lambda x: x.astype(int)).to_string())
    out.write("\n\n")

    out.write("Final unique markers by chromosome:\n")
    out.write(final_unique["chr"].value_counts().sort_index(key=lambda x: x.astype(int)).to_string())
    out.write("\n\n")

    out.write("RefSeq accession to chromosome mapping used:\n")
    for ref, chrom in ref_to_chr.items():
        out.write(f"{ref}\tchr{chrom}\t{ref_to_lg[ref]}\n")

print(f"Clean map with marker info: {clean_with_markers_out}")
print(f"Final unique chr-pos-cM map: {final_unique_out}")
print(f"Excluded scaffolds: {excluded_scaffolds_out}")
print(f"Excluded LG/ref_seq conflicts: {excluded_conflicts_out}")
print(f"Duplicate positions: {duplicates_out}")
print(f"Summary: {summary_out}")

print("\nSummary:")
print(f"Input strict markers: {len(df)}")
print(f"Markers on scaffold/contig accessions excluded: {len(scaffolds)}")
print(f"LG/ref_seq conflict markers excluded: {len(conflicts)}")
print(f"Clean markers after filtering: {len(clean)}")
print(f"Duplicate-position rows in clean map: {len(duplicates)}")
print(f"Final unique chr-pos-cM rows: {len(final_unique)}")

print("\nClean markers by chromosome:")
print(clean["chr"].value_counts().sort_index(key=lambda x: x.astype(int)).to_string())

print("\nFirst rows of final unique map:")
print(final_unique.head(10).to_string(index=False))

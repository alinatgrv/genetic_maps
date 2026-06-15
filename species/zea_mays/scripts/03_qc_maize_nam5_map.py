from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SPECIES_DIR = PROJECT_ROOT / "species" / "zea_mays"

INTERMEDIATE = SPECIES_DIR / "results" / "intermediate" / "maizegdb_genetic1_nam5_with_positions.tsv"
FINAL = SPECIES_DIR / "results" / "final" / "zea_mays_genetic_map.tsv"
QC_DIR = SPECIES_DIR / "results" / "qc"

OUT_SUMMARY = QC_DIR / "maize_nam5_map_qc_summary.txt"
OUT_DUP_POS = QC_DIR / "maize_nam5_duplicate_chr_pos.tsv"
OUT_CHR_LG = QC_DIR / "maize_nam5_chr_vs_linkage_group.tsv"
OUT_INVALID = QC_DIR / "maize_nam5_rows_without_valid_position.tsv"
OUT_CORR = QC_DIR / "maize_nam5_physical_vs_genetic_correlation.tsv"

QC_DIR.mkdir(parents=True, exist_ok=True)

def pearson_corr(a, b):
    tmp = pd.DataFrame({
        "a": pd.to_numeric(a, errors="coerce"),
        "b": pd.to_numeric(b, errors="coerce"),
    }).dropna()

    if len(tmp) < 3:
        return pd.NA
    if tmp["a"].nunique() < 2 or tmp["b"].nunique() < 2:
        return pd.NA

    return float(tmp["a"].corr(tmp["b"], method="pearson"))

def spearman_corr_without_scipy(a, b):
    tmp = pd.DataFrame({
        "a": pd.to_numeric(a, errors="coerce"),
        "b": pd.to_numeric(b, errors="coerce"),
    }).dropna()

    if len(tmp) < 3:
        return pd.NA
    if tmp["a"].nunique() < 2 or tmp["b"].nunique() < 2:
        return pd.NA

    ar = tmp["a"].rank(method="average")
    br = tmp["b"].rank(method="average")
    return float(ar.corr(br, method="pearson"))

df = pd.read_csv(INTERMEDIATE, sep="\t")
final = pd.read_csv(FINAL, sep="\t")

# has_nam5_position sometimes can be read as bool, sometimes as string
has_pos = df["has_nam5_position"].astype(str).str.lower().isin(["true", "1", "yes"])
valid = df[has_pos].copy()
invalid = df[~has_pos].copy()

for col in ["chr", "start", "end", "cM", "linkage_group"]:
    valid[col] = pd.to_numeric(valid[col], errors="coerce")

valid = valid.dropna(subset=["chr", "start", "end", "cM", "linkage_group"]).copy()

valid["pos"] = ((valid["start"] + valid["end"]) / 2).round().astype("Int64")
valid["chr"] = valid["chr"].astype("Int64")
valid["linkage_group"] = valid["linkage_group"].astype("Int64")

dup_pos = (
    valid[valid.duplicated(["chr", "pos"], keep=False)]
    .sort_values(["chr", "pos", "cM", "locus"])
)

chr_lg = pd.crosstab(valid["chr"], valid["linkage_group"])

rows = []
for chr_id, sub in valid.groupby("chr"):
    sub = sub.dropna(subset=["pos", "cM"]).sort_values("pos")
    if len(sub) < 3:
        continue

    dcM = sub["cM"].diff()
    decreases = int((dcM < 0).sum())

    rows.append({
        "chr": int(chr_id),
        "n": len(sub),
        "spearman_pos_cM": spearman_corr_without_scipy(sub["pos"], sub["cM"]),
        "pearson_pos_cM": pearson_corr(sub["pos"], sub["cM"]),
        "cM_decreases_after_sort_by_pos": decreases,
        "min_pos": int(sub["pos"].min()),
        "max_pos": int(sub["pos"].max()),
        "min_cM": float(sub["cM"].min()),
        "max_cM": float(sub["cM"].max()),
    })

corr_df = pd.DataFrame(rows).sort_values("chr")

dup_pos.to_csv(OUT_DUP_POS, sep="\t", index=False)
chr_lg.to_csv(OUT_CHR_LG, sep="\t")
invalid.to_csv(OUT_INVALID, sep="\t", index=False)
corr_df.to_csv(OUT_CORR, sep="\t", index=False)

with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
    out.write("Maize NAM-5.0 map QC summary\n")
    out.write("============================\n\n")

    out.write(f"Intermediate file: {INTERMEDIATE}\n")
    out.write(f"Final file: {FINAL}\n\n")

    out.write(f"Rows in intermediate table: {len(df)}\n")
    out.write(f"Rows with valid NAM-5.0 position: {len(valid)}\n")
    out.write(f"Rows without valid NAM-5.0 position: {len(invalid)}\n")
    out.write(f"Rows in final chr-pos-cM table: {len(final)}\n\n")

    out.write("Final map rows by chromosome:\n")
    out.write(final.groupby("chr").size().to_string())
    out.write("\n\n")

    out.write("Duplicate chr-pos rows before final deduplication:\n")
    out.write(f"{len(dup_pos)} rows\n\n")

    out.write("Chromosome vs linkage group crosstab:\n")
    out.write(chr_lg.to_string())
    out.write("\n\n")

    out.write("Physical vs genetic coordinate correlation:\n")
    out.write(corr_df.to_string(index=False))
    out.write("\n\n")

    out.write("QC output files:\n")
    out.write(f"- {OUT_DUP_POS}\n")
    out.write(f"- {OUT_CHR_LG}\n")
    out.write(f"- {OUT_INVALID}\n")
    out.write(f"- {OUT_CORR}\n")

print("Done.")
print(f"QC summary: {OUT_SUMMARY}")
print(f"Duplicate positions: {OUT_DUP_POS}")
print(f"Invalid rows: {OUT_INVALID}")
print(f"Correlation table: {OUT_CORR}")

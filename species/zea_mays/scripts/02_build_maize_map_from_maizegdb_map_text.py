from pathlib import Path
import re
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SPECIES_DIR = PROJECT_ROOT / "species" / "zea_mays"

RAW_DIR = SPECIES_DIR / "data" / "raw" / "maizegdb_map_1203637" / "map_text_browser"

INTERMEDIATE_DIR = SPECIES_DIR / "results" / "intermediate"
FINAL_DIR = SPECIES_DIR / "results" / "final"
QC_DIR = SPECIES_DIR / "results" / "qc"

for d in [INTERMEDIATE_DIR, FINAL_DIR, QC_DIR]:
    d.mkdir(parents=True, exist_ok=True)

OUT_ALL = INTERMEDIATE_DIR / "maizegdb_genetic1_all_parsed_rows.tsv"
OUT_NAM = INTERMEDIATE_DIR / "maizegdb_genetic1_nam5_with_positions.tsv"
OUT_FINAL = FINAL_DIR / "zea_mays_genetic_map.tsv"
OUT_COUNTS = QC_DIR / "maize_marker_counts_by_chr.tsv"
OUT_SUMMARY = QC_DIR / "maizegdb_genetic1_nam5_summary.txt"

TARGET = "Zm-B73-REFERENCE-NAM-5.0"

REQUIRED = [
    "Locus",
    "Coordinate",
    f"{TARGET}_gene_model",
    f"{TARGET}_chr",
    f"{TARGET}_start",
    f"{TARGET}_end",
]

def read_maizegdb_text(path: Path) -> pd.DataFrame:
    """
    MaizeGDB map_text files have:
    line 1: title
    line 2: tab-delimited header
    following lines: tab-delimited data
    """
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()

    if len(lines) < 3:
        raise ValueError(f"Too few lines in {path}")

    header = [x.strip() for x in lines[1].split("\t")]
    rows = []

    for line in lines[2:]:
        if not line.strip():
            continue

        parts = line.split("\t")

        # pad short rows caused by empty trailing fields
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))

        # truncate long rows; target columns are before possible Sequence text
        if len(parts) > len(header):
            parts = parts[:len(header)]

        rows.append(parts)

    df = pd.DataFrame(rows, columns=header)
    return df

def clean_int(x):
    if pd.isna(x):
        return pd.NA
    s = str(x).strip().replace(",", "")
    if s == "" or s.lower() in {"nan", "none"}:
        return pd.NA
    s = re.sub(r"[^0-9]", "", s)
    if s == "":
        return pd.NA
    try:
        return int(s)
    except ValueError:
        return pd.NA

def clean_chr(x):
    if pd.isna(x):
        return pd.NA
    s = str(x).strip()
    if s == "":
        return pd.NA
    s = re.sub(r"^chr", "", s, flags=re.IGNORECASE)
    if not re.fullmatch(r"\d+", s):
        return pd.NA
    return int(s)

files = sorted(RAW_DIR.glob("maizegdb_genetic1_*_map_text_BROWSER.txt"))
if not files:
    raise FileNotFoundError(f"No browser map_text files found in {RAW_DIR}")

all_frames = []
problems = []

for path in files:
    m = re.search(r"maizegdb_genetic1_(\d+)_map_text_BROWSER\.txt$", path.name)
    map_id = m.group(1) if m else path.stem

    lg = int(map_id) - 1203636 if map_id.isdigit() else pd.NA

    df = read_maizegdb_text(path)
    df.columns = [str(c).strip() for c in df.columns]

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        problems.append((path.name, missing))
        continue

    df["source_file"] = path.name
    df["map_id"] = map_id
    df["linkage_group"] = lg

    all_frames.append(df)

if problems:
    msg = "\n".join([f"{fn}: missing {missing}" for fn, missing in problems])
    raise ValueError("Some files do not contain required columns:\n" + msg)

all_df = pd.concat(all_frames, ignore_index=True)
all_df.to_csv(OUT_ALL, sep="\t", index=False)

nam = pd.DataFrame({
    "source_file": all_df["source_file"],
    "map_id": all_df["map_id"],
    "linkage_group": all_df["linkage_group"],
    "locus": all_df["Locus"],
    "cM": pd.to_numeric(all_df["Coordinate"], errors="coerce"),
    "bin": all_df.get("Bin", pd.NA),
    "gene_model": all_df[f"{TARGET}_gene_model"],
    "chr": all_df[f"{TARGET}_chr"].map(clean_chr),
    "start": all_df[f"{TARGET}_start"].map(clean_int),
    "end": all_df[f"{TARGET}_end"].map(clean_int),
})

# Basic filtering
nam["has_nam5_position"] = (
    nam["chr"].notna()
    & nam["start"].notna()
    & nam["end"].notna()
    & nam["cM"].notna()
)

# Avoid impossible corrupted coordinates if any parsing issue appears
MAX_REASONABLE_POS = 400_000_000
nam.loc[nam["start"] > MAX_REASONABLE_POS, "has_nam5_position"] = False
nam.loc[nam["end"] > MAX_REASONABLE_POS, "has_nam5_position"] = False

nam_valid = nam[nam["has_nam5_position"]].copy()

# If start > end, normalize interval
swap_mask = nam_valid["start"] > nam_valid["end"]
nam_valid.loc[swap_mask, ["start", "end"]] = nam_valid.loc[swap_mask, ["end", "start"]].to_numpy()

nam_valid["pos"] = ((nam_valid["start"] + nam_valid["end"]) / 2).round().astype("int64")
nam_valid["chr"] = nam_valid["chr"].astype("int64")

nam.to_csv(OUT_NAM, sep="\t", index=False)

final = (
    nam_valid[["chr", "pos", "cM"]]
    .drop_duplicates()
    .sort_values(["chr", "pos", "cM"])
    .reset_index(drop=True)
)

final.to_csv(OUT_FINAL, sep="\t", index=False)

counts = (
    nam_valid
    .groupby("chr", as_index=False)
    .agg(
        markers=("locus", "count"),
        unique_loci=("locus", "nunique"),
        unique_positions=("pos", "nunique"),
        min_pos=("pos", "min"),
        max_pos=("pos", "max"),
        min_cM=("cM", "min"),
        max_cM=("cM", "max"),
    )
    .sort_values("chr")
)

counts.to_csv(OUT_COUNTS, sep="\t", index=False)

with open(OUT_SUMMARY, "w", encoding="utf-8") as out:
    out.write("MaizeGDB Genetic 1 -> NAM-5.0 genetic map summary\n")
    out.write("=================================================\n\n")
    out.write(f"Input directory: {RAW_DIR}\n")
    out.write(f"Input files: {len(files)}\n")
    out.write(f"Target assembly columns: {TARGET}\n\n")

    out.write(f"All parsed MaizeGDB rows: {len(all_df)}\n")
    out.write(f"Rows with numeric cM: {nam['cM'].notna().sum()}\n")
    out.write(f"Rows with NAM-5.0 chr/start/end/cM: {len(nam_valid)}\n")
    out.write(f"Final unique chr-pos-cM rows: {len(final)}\n\n")

    out.write("Rows by source map id / linkage group:\n")
    tmp = nam.groupby(["map_id", "linkage_group"], as_index=False).size()
    out.write(tmp.to_string(index=False))
    out.write("\n\n")

    out.write("Valid NAM-5.0 rows by chromosome:\n")
    out.write(counts.to_string(index=False))
    out.write("\n\n")

    out.write("Output files:\n")
    out.write(f"- {OUT_ALL}\n")
    out.write(f"- {OUT_NAM}\n")
    out.write(f"- {OUT_FINAL}\n")
    out.write(f"- {OUT_COUNTS}\n")

print("Done.")
print(f"All parsed rows: {len(all_df)}")
print(f"Rows with NAM-5.0 positions: {len(nam_valid)}")
print(f"Final unique chr-pos-cM rows: {len(final)}")
print()
print(f"Final map: {OUT_FINAL}")
print(f"Summary:   {OUT_SUMMARY}")

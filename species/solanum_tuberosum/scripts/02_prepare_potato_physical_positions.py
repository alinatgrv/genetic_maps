#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE = Path("species/solanum_tuberosum")

RAW_DIR = BASE / "data" / "raw" / "spuddb_dm_v6_1"
OUT_DIR = BASE / "results" / "intermediate"
QC_DIR = BASE / "results" / "qc"

OUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

inputs = [
    ("solcap_69k", RAW_DIR / "solcap_69k_SNPs_DM_v6_1_pos.txt.gz"),
    ("potvar", RAW_DIR / "potvar_SNPs_DM_v6_1_pos.txt.gz"),
]

frames = []

for source, path in inputs:
    df = pd.read_csv(path, sep=r"\s+", compression="gzip")
    df["source"] = source

    df = df.rename(
        columns={
            "SNP_ID": "marker_id",
            "REF": "chr",
            "SNP_POS": "pos",
            "NUM_MISMATCH": "num_mismatch",
            "MULTIMAP": "multimap",
        }
    )

    df["marker_id"] = df["marker_id"].astype(str)
    df["chr"] = df["chr"].astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df["num_mismatch"] = pd.to_numeric(df["num_mismatch"], errors="coerce").astype("Int64")
    df["multimap"] = df["multimap"].astype(str)

    frames.append(df)

all_pos = pd.concat(frames, ignore_index=True)

# Базовая таблица всех позиций, включая мультимэпперы.
all_pos = all_pos[
    ["marker_id", "chr", "pos", "num_mismatch", "multimap", "source"]
].sort_values(["source", "marker_id", "chr", "pos"])

all_out = OUT_DIR / "potato_physical_positions.all.tsv"
all_pos.to_csv(all_out, sep="\t", index=False)

# Надежные физические координаты: маркер не мультимэппер.
unique_pos = all_pos[all_pos["multimap"] == "N"].copy()

# Проверка: после удаления мультимэпперов каждый marker_id должен иметь одну позицию.
marker_counts = unique_pos.groupby("marker_id").size().reset_index(name="n_positions")
duplicated_unique = marker_counts[marker_counts["n_positions"] > 1]

duplicated_out = QC_DIR / "potato_physical_positions.unique_marker_duplicates.tsv"
duplicated_unique.to_csv(duplicated_out, sep="\t", index=False)

unique_out = OUT_DIR / "potato_physical_positions.unique.tsv"
unique_pos.to_csv(unique_out, sep="\t", index=False)

# Самая строгая версия: только уникальные и без mismatch.
strict_exact = unique_pos[unique_pos["num_mismatch"] == 0].copy()
strict_out = OUT_DIR / "potato_physical_positions.strict_exact.tsv"
strict_exact.to_csv(strict_out, sep="\t", index=False)

# QC summary.
summary_lines = []
summary_lines.append("Potato physical SNP positions on DM v6.1")
summary_lines.append("=" * 45)
summary_lines.append("")
summary_lines.append(f"Input rows total: {len(all_pos)}")
summary_lines.append(f"Unique-position rows, MULTIMAP == N: {len(unique_pos)}")
summary_lines.append(f"Strict exact rows, MULTIMAP == N and NUM_MISMATCH == 0: {len(strict_exact)}")
summary_lines.append(f"Duplicated marker_id among MULTIMAP == N rows: {len(duplicated_unique)}")
summary_lines.append("")

summary_lines.append("Rows by source:")
summary_lines.append(all_pos.groupby("source").size().to_string())
summary_lines.append("")

summary_lines.append("Rows by source and multimap:")
summary_lines.append(all_pos.groupby(["source", "multimap"]).size().to_string())
summary_lines.append("")

summary_lines.append("Rows by source and num_mismatch:")
summary_lines.append(all_pos.groupby(["source", "num_mismatch"]).size().to_string())
summary_lines.append("")

summary_lines.append("Unique-position rows by chromosome:")
summary_lines.append(unique_pos.groupby("chr").size().sort_index().to_string())
summary_lines.append("")

summary_file = QC_DIR / "potato_physical_positions_summary.txt"
summary_file.write_text("\n".join(summary_lines) + "\n")

print("Done.")
print(f"All positions:          {all_out}")
print(f"Unique positions:       {unique_out}")
print(f"Strict exact positions: {strict_out}")
print(f"QC summary:             {summary_file}")
print(f"Duplicate check:        {duplicated_out}")

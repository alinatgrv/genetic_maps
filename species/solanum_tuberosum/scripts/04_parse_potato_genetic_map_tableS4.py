#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import pandas as pd

BASE = Path("species/solanum_tuberosum")

RAW_DIR = BASE / "data" / "raw" / "pgsc_sharma_2013"
PHYSICAL_FILE = BASE / "results" / "intermediate" / "potato_physical_positions.unique.tsv"

OUT_DIR = BASE / "results" / "intermediate"
QC_DIR = BASE / "results" / "qc"

OUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Поддерживаем оба варианта: TSV или CSV.
candidate_files = [
    RAW_DIR / "TableS4.tsv",
    RAW_DIR / "TableS4.csv",
    RAW_DIR / "supp_g3.113.007153_TableS4.tsv",
    RAW_DIR / "supp_g3.113.007153_TableS4.csv",
]

table_file = None
for p in candidate_files:
    if p.exists():
        table_file = p
        break

if table_file is None:
    raise FileNotFoundError(
        "Не найден TableS4.tsv/TableS4.csv в species/solanum_tuberosum/data/raw/pgsc_sharma_2013/"
    )

print(f"Using genetic map file: {table_file}")

def detect_delimiter(path: Path) -> str:
    # Если явно TSV — берем tab.
    if path.suffix.lower() == ".tsv":
        return "\t"

    # Иначе пытаемся угадать.
    with path.open("r", encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(4096)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except csv.Error:
        return ","

def parse_attrs(attr: str) -> dict:
    result = {}
    if not isinstance(attr, str):
        return result

    # GFF attributes обычно разделены ;
    parts = [x.strip() for x in attr.split(";") if x.strip()]
    for part in parts:
        if "=" in part:
            k, v = part.split("=", 1)
            result[k.strip()] = v.strip()

    return result

def normalize_chr(x: str) -> str:
    if pd.isna(x):
        return ""

    x = str(x).strip()

    # Chr01 / chr01 / chromosome 01 -> chr01
    m = re.search(r"(\d+)", x)
    if m:
        return f"chr{int(m.group(1)):02d}"

    return x

def normalize_marker_for_join(marker_id: str) -> str:
    """
    В TableS4 SolCAP SNP часто записаны как solcap_stsnp_*,
    а в v6.1 physical file как solcap_snp_*.
    """
    marker_id = str(marker_id).strip()
    marker_id = marker_id.replace("solcap_stsnp_", "solcap_snp_")
    return marker_id

delimiter = detect_delimiter(table_file)

rows = []
with table_file.open("r", encoding="utf-8-sig", errors="replace", newline="") as f:
    reader = csv.reader(f, delimiter=delimiter)

    for fields in reader:
        if not fields:
            continue

        # Убираем пустые хвостовые колонки после экспорта из Excel.
        fields = [x for x in fields if str(x).strip() != ""]

        if not fields:
            continue

        first = str(fields[0]).strip()

        # Комментарии и заголовки пропускаем, кроме реальных GFF-строк.
        if first.startswith("#"):
            continue

        if len(fields) < 9:
            continue

        seqid, source, marker_type, start, end, score, strand, phase, attrs = fields[:9]

        attrs_dict = parse_attrs(attrs)

        marker_id = attrs_dict.get("ID") or attrs_dict.get("Name") or attrs_dict.get("NAME")
        marker_name = attrs_dict.get("Name") or attrs_dict.get("NAME") or marker_id

        if marker_id is None:
            continue

        # Ищем "Genetic position = Chr01 8.698 cM" в attributes.
        m = re.search(
            r"Genetic position\s*=\s*(Chr\s*0?\d+|Chr0?\d+|chr0?\d+)\s+([0-9]+(?:\.[0-9]+)?)\s*cM",
            attrs,
            flags=re.IGNORECASE,
        )

        if not m:
            continue

        genetic_chr_raw = m.group(1).replace(" ", "")
        cM = float(m.group(2))

        rows.append(
            {
                "marker_id_original": str(marker_id).strip(),
                "marker_id": normalize_marker_for_join(marker_id),
                "marker_name": str(marker_name).strip(),
                "genetic_chr": normalize_chr(genetic_chr_raw),
                "cM": cM,
                "tableS4_seqid_v403": seqid,
                "tableS4_start_v403": start,
                "tableS4_end_v403": end,
                "marker_type": marker_type,
                "attributes": attrs,
            }
        )

genetic = pd.DataFrame(rows)

if genetic.empty:
    raise RuntimeError("Не удалось извлечь ни одной строки с genetic position из TableS4.")

genetic = genetic.sort_values(["genetic_chr", "cM", "marker_id_original"])

genetic_out = OUT_DIR / "potato_genetic_map.sharma2013_TableS4.tsv"
genetic.to_csv(genetic_out, sep="\t", index=False)

# Отдельно SolCAP SNP, потому что именно они должны пересекаться с physical v6.1.
genetic_solcap = genetic[
    genetic["marker_id_original"].str.contains("solcap", case=False, na=False)
].copy()

genetic_solcap_out = OUT_DIR / "potato_genetic_map.sharma2013_TableS4_solcap.tsv"
genetic_solcap.to_csv(genetic_solcap_out, sep="\t", index=False)

# Загружаем физические позиции на DM v6.1.
physical = pd.read_csv(PHYSICAL_FILE, sep="\t")
physical["marker_id"] = physical["marker_id"].astype(str)
physical["chr"] = physical["chr"].map(normalize_chr)

# Проверка пересечения.
genetic_ids_original = set(genetic["marker_id_original"].astype(str))
genetic_ids_normalized = set(genetic["marker_id"].astype(str))
physical_ids = set(physical["marker_id"].astype(str))

direct_overlap = genetic_ids_original & physical_ids
normalized_overlap = genetic_ids_normalized & physical_ids

merged = genetic.merge(
    physical,
    on="marker_id",
    how="inner",
    suffixes=("_genetic", "_physical"),
)

merged["chr_agrees"] = merged["genetic_chr"] == merged["chr"]

merged_out = OUT_DIR / "potato_genetic_map.sharma2013_TableS4_merged_with_DM_v6_1_physical.tsv"
merged.to_csv(merged_out, sep="\t", index=False)

merged_chr_agree = merged[merged["chr_agrees"]].copy()
merged_chr_agree_out = OUT_DIR / "potato_genetic_map.sharma2013_TableS4_merged_chr_agree.tsv"
merged_chr_agree.to_csv(merged_chr_agree_out, sep="\t", index=False)

# Финальная предварительная карта: chr pos cM, только где genetic chr == physical chr.
final_prelim = merged_chr_agree[["chr", "pos", "cM"]].copy()
final_prelim = final_prelim.drop_duplicates().sort_values(["chr", "pos", "cM"])

final_prelim_out = OUT_DIR / "potato_genetic_map.preliminary_chr_pos_cM.tsv"
final_prelim.to_csv(final_prelim_out, sep="\t", index=False)

summary = []
summary.append("Potato Sharma et al. 2013 TableS4 genetic map parsing")
summary.append("=" * 60)
summary.append("")
summary.append(f"Input TableS4 file: {table_file}")
summary.append(f"Parsed genetic map rows: {len(genetic)}")
summary.append(f"Parsed SolCAP rows: {len(genetic_solcap)}")
summary.append(f"Physical DM v6.1 unique-position rows: {len(physical)}")
summary.append("")
summary.append(f"Direct marker_id overlap with physical v6.1: {len(direct_overlap)}")
summary.append(f"Normalized marker_id overlap with physical v6.1: {len(normalized_overlap)}")
summary.append(f"Merged rows after normalized marker_id join: {len(merged)}")
summary.append(f"Merged rows with genetic_chr == physical chr: {len(merged_chr_agree)}")
summary.append(f"Preliminary final chr-pos-cM rows: {len(final_prelim)}")
summary.append("")

summary.append("Marker types in parsed genetic map:")
summary.append(genetic.groupby("marker_type").size().sort_values(ascending=False).to_string())
summary.append("")

summary.append("Parsed genetic rows by chromosome:")
summary.append(genetic.groupby("genetic_chr").size().sort_index().to_string())
summary.append("")

summary.append("Merged rows by physical chromosome:")
if not merged.empty:
    summary.append(merged.groupby("chr").size().sort_index().to_string())
else:
    summary.append("No merged rows")
summary.append("")

summary.append("Merged rows with chromosome agreement by chromosome:")
if not merged_chr_agree.empty:
    summary.append(merged_chr_agree.groupby("chr").size().sort_index().to_string())
else:
    summary.append("No chr-agree rows")
summary.append("")

summary_file = QC_DIR / "potato_sharma2013_TableS4_overlap_summary.txt"
summary_file.write_text("\n".join(summary) + "\n")

print("Done.")
print(f"Parsed genetic map:          {genetic_out}")
print(f"Parsed SolCAP genetic map:   {genetic_solcap_out}")
print(f"Merged with physical v6.1:   {merged_out}")
print(f"Merged chr-agree only:       {merged_chr_agree_out}")
print(f"Preliminary chr-pos-cM map:  {final_prelim_out}")
print(f"QC summary:                  {summary_file}")

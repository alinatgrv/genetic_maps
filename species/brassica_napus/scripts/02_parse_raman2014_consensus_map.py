#!/usr/bin/env python3

# Парсинг консенсус-карты DArT рапса (Raman et al. 2014, BMC Genomics 14:277,
# Additional file 1, лист "Consensus Map").
#
# Из листа извлекаем для каждого маркера: имя, хромосому (A01..A10/C1..C9),
# генетическую позицию cM и класс (DArT / Non-DArT).
#
# Названия хромосом нормализуем к стилю целевого референса Da-Ae (A1..A10, C1..C9).
# Для DArT-маркеров формируем join-id (нижний регистр, без ведущего 'X'),
# по которому затем подбираются последовательности из DArT_Brassica.fasta.
#
# Запускать из корня репозитория.

from pathlib import Path
from collections import Counter
import re
import xlrd

SPECIES_DIR = Path("species/brassica_napus")

xls_file = SPECIES_DIR / "data/raw/raman2014_dart/1471-2164-14-277-S1.xls"
metadata_out = SPECIES_DIR / "data/metadata/bnapus_consensus_map_metadata.tsv"
dart_ids_out = SPECIES_DIR / "data/metadata/bnapus_consensus_dart_join_ids.txt"
summary_out = SPECIES_DIR / "results/qc/bnapus_consensus_map_summary.txt"

metadata_out.parent.mkdir(parents=True, exist_ok=True)
summary_out.parent.mkdir(parents=True, exist_ok=True)

CONSENSUS_SHEET = "Consensus Map"
# В листе "Consensus Map": строка 0 — заголовок-описание, строка 1 — шапка,
# данные начинаются со строки 2 (0-индексация).
DATA_START_ROW = 2


def normalize_chr(raw):
    """A01 -> A1, C09 -> C9. Возвращает None для не-хромосомных значений."""
    s = str(raw).strip()
    m = re.fullmatch(r"([ACac])0*([0-9]+)", s)
    if not m:
        return None
    return f"{m.group(1).upper()}{int(m.group(2))}"


def join_id(marker):
    """Нормализованный ключ для джойна с FASTA: нижний регистр, без ведущего 'x'."""
    return re.sub(r"^x", "", str(marker).strip().lower())


book = xlrd.open_workbook(xls_file)
if CONSENSUS_SHEET not in book.sheet_names():
    raise ValueError(f"Лист '{CONSENSUS_SHEET}' не найден; есть: {book.sheet_names()}")
sheet = book.sheet_by_name(CONSENSUS_SHEET)

rows = []
skipped = 0
for r in range(DATA_START_ROW, sheet.nrows):
    marker = str(sheet.cell_value(r, 0)).strip()
    chrom_raw = sheet.cell_value(r, 1)
    pos_raw = sheet.cell_value(r, 2)
    marker_class = str(sheet.cell_value(r, 3)).strip()

    if not marker:
        skipped += 1
        continue

    chrom = normalize_chr(chrom_raw)
    if chrom is None:
        skipped += 1
        continue

    try:
        cM = float(pos_raw)
    except (TypeError, ValueError):
        skipped += 1
        continue

    jid = join_id(marker)
    is_dart = (marker_class.lower() == "dart") or jid.startswith("brpb-")

    rows.append({
        "marker_id": marker,
        "marker_join_id": jid,
        "chr": chrom,
        "cM": f"{cM:.6f}".rstrip("0").rstrip(".") if "." in f"{cM:.6f}" else f"{cM:.6f}",
        "cM_value": cM,
        "marker_class": "DArT" if is_dart else "Non-DArT",
        "is_dart": "1" if is_dart else "0",
    })

# Запись метаданных карты (все маркеры).
header = ["marker_id", "marker_join_id", "chr", "cM", "marker_class", "is_dart"]
with open(metadata_out, "w", encoding="utf-8") as out:
    out.write("\t".join(header) + "\n")
    for row in rows:
        out.write("\t".join(str(row[h]) for h in header) + "\n")

# Уникальные join-id DArT-маркеров для извлечения последовательностей.
dart_join_ids = sorted({row["marker_join_id"] for row in rows if row["is_dart"] == "1"})
with open(dart_ids_out, "w", encoding="utf-8") as out:
    for jid in dart_join_ids:
        out.write(jid + "\n")

# Сводка.
chr_order = [f"A{i}" for i in range(1, 11)] + [f"C{i}" for i in range(1, 10)]
chr_counts = Counter(row["chr"] for row in rows)
class_counts = Counter(row["marker_class"] for row in rows)
dart_chr_counts = Counter(row["chr"] for row in rows if row["is_dart"] == "1")

with open(summary_out, "w", encoding="utf-8") as out:
    out.write("Brassica napus DArT consensus map (Raman et al. 2014, Additional file 1)\n")
    out.write("=" * 70 + "\n\n")
    out.write(f"Input file: {xls_file}\n")
    out.write(f"Sheet: {CONSENSUS_SHEET}\n")
    out.write(f"Parsed marker rows: {len(rows)}\n")
    out.write(f"Skipped rows: {skipped}\n")
    out.write(f"Unique DArT join-ids: {len(dart_join_ids)}\n\n")

    out.write("Markers by class:\n")
    for cls, n in class_counts.most_common():
        out.write(f"  {cls}\t{n}\n")

    out.write("\nMarkers by chromosome (all / DArT):\n")
    for ch in chr_order:
        out.write(f"  {ch}\t{chr_counts.get(ch, 0)}\t{dart_chr_counts.get(ch, 0)}\n")

print("Done.")
print(f"Map metadata:   {metadata_out}")
print(f"DArT join-ids:  {dart_ids_out}")
print(f"Summary:        {summary_out}")
print(f"Total markers: {len(rows)} | DArT: {len(dart_join_ids)}")

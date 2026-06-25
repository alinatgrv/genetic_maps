#!/usr/bin/env python3

from pathlib import Path
import csv
import re
import xlrd

XLS = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/pone.0017595.s004.xls")
OUT = Path("species/oryza_sativa_japonica/data/raw/public_sources/published_maps/yu_2011_plosone/yu2011_bins_tigr6.tsv")
QC = Path("species/oryza_sativa_japonica/results/qc/yu2011_bins_extract_summary.txt")

OUT.parent.mkdir(parents=True, exist_ok=True)
QC.parent.mkdir(parents=True, exist_ok=True)

book = xlrd.open_workbook(str(XLS))
sheet = book.sheet_by_index(0)

print(f"Workbook: {XLS}")
print(f"Sheet: {sheet.name}")
print(f"Rows: {sheet.nrows}, Cols: {sheet.ncols}")

# Ищем строку с заголовками.
header_row = None
header = None

for i in range(min(sheet.nrows, 50)):
    vals = [str(sheet.cell_value(i, j)).strip() for j in range(sheet.ncols)]
    joined = "\t".join(vals).lower()
    if "start" in joined and "stop" in joined:
        header_row = i
        header = vals
        break

if header_row is None:
    raise SystemExit("Could not find header row with Start/Stop")

def norm(x):
    return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")

headers = [norm(x) for x in header]
print("Header row:", header_row + 1)
print("Headers:", headers)

def find_col(patterns):
    for p in patterns:
        for idx, h in enumerate(headers):
            if re.search(p, h):
                return idx
    return None

bin_col = find_col([r"^bin$"])
chr_col = find_col([r"^chr", r"chrom"])
start_col = find_col([r"start"])
stop_col = find_col([r"stop", r"end"])
cm_col = find_col([r"genetic.*position", r"genetic", r"cm", r"position"])

# Иногда в старых Excel заголовки могут быть неполными, поэтому задаем fallback:
# ожидаемая структура Table S2: Bin, Chr, Start, Stop, Length, Genetic position
if bin_col is None:
    bin_col = 0
if chr_col is None:
    chr_col = 1
if start_col is None:
    start_col = 2
if stop_col is None:
    stop_col = 3
if cm_col is None:
    cm_col = 5

print("Columns used:")
print("bin_col", bin_col, "chr_col", chr_col, "start_col", start_col, "stop_col", stop_col, "cm_col", cm_col)

rows = []

with OUT.open("w", newline="") as w:
    fieldnames = [
        "bin",
        "chr",
        "start_mb",
        "stop_mb",
        "length_mb",
        "cM",
        "old_start_bp_1based",
        "old_stop_bp_1based",
        "old_mid_bp_1based",
    ]
    writer = csv.DictWriter(w, delimiter="\t", fieldnames=fieldnames)
    writer.writeheader()

    for i in range(header_row + 1, sheet.nrows):
        vals = [sheet.cell_value(i, j) for j in range(sheet.ncols)]

        bin_id = str(vals[bin_col]).strip()
        if not bin_id or bin_id.lower() == "nan":
            continue
        if not bin_id.lower().startswith("bin"):
            continue

        try:
            chr_raw = vals[chr_col]
            start_mb = float(vals[start_col])
            stop_mb = float(vals[stop_col])
            cm = float(vals[cm_col])
        except Exception:
            continue

        chr_num = str(int(float(chr_raw)))
        old_start = int(round(start_mb * 1_000_000)) + 1
        old_stop = int(round(stop_mb * 1_000_000))
        old_mid = int(round(((start_mb + stop_mb) / 2) * 1_000_000))

        row = {
            "bin": bin_id,
            "chr": chr_num,
            "start_mb": f"{start_mb:.6f}",
            "stop_mb": f"{stop_mb:.6f}",
            "length_mb": f"{(stop_mb - start_mb):.6f}",
            "cM": f"{cm:.6f}",
            "old_start_bp_1based": old_start,
            "old_stop_bp_1based": old_stop,
            "old_mid_bp_1based": old_mid,
        }
        writer.writerow(row)
        rows.append(row)

if not rows:
    raise SystemExit("No bin rows extracted")

with QC.open("w") as w:
    w.write(f"input_xls\t{XLS}\n")
    w.write(f"sheet\t{sheet.name}\n")
    w.write(f"header_row_1based\t{header_row + 1}\n")
    w.write(f"n_bins\t{len(rows)}\n")
    w.write(f"chromosomes\t{','.join(sorted(set(r['chr'] for r in rows), key=lambda x: int(x)))}\n")
    w.write(f"cM_min\t{min(float(r['cM']) for r in rows):.6f}\n")
    w.write(f"cM_max\t{max(float(r['cM']) for r in rows):.6f}\n")

print(f"Written: {OUT}")
print(f"Bins: {len(rows)}")
print(f"QC: {QC}")

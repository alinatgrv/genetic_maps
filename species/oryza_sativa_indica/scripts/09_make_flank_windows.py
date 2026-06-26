#!/usr/bin/env python3

# Извлечение коротких фланкирующих окон вокруг bin-серединок на IRGSP-1.0.
#
# Вместо медленного полногеномного выравнивания двух разошедшихся рисовых геномов
# переносим 1619 точек локально: вокруг каждой точки берём окно ±W bp на IRGSP-1.0,
# затем эти окна выравниваются на ZS97RS3 (minimap2), и центр окна переносится
# через paftools.js liftover. Результат — те же координаты bin на ZS97RS3, но быстро.
#
# На выходе:
#   windows_irgsp1_flank.fasta  — окна (заголовок = имя bin: "Bin5|src=1|cM=...")
#   windows_center.bed          — позиция центра (исходной точки) внутри окна
#
# Запускать из корня репозитория (нужен samtools на PATH; IRGSP-1.0 проиндексирован faidx).

from pathlib import Path
import subprocess

SPECIES_DIR = Path("species/oryza_sativa_indica")
IRGSP_FNA = SPECIES_DIR / "data/ref/GCF_001433935.1_genomic.fna"
BIN_BED = SPECIES_DIR / "data/raw/yu2011_bins_irgsp1_midpoints.bed"

OUT_DIR = SPECIES_DIR / "results/liftover/yu2011_irgsp1_to_zs97rs3"
OUT_DIR.mkdir(parents=True, exist_ok=True)
WINDOWS_FA = OUT_DIR / "windows_irgsp1_flank.fasta"
CENTER_BED = OUT_DIR / "windows_center.bed"
REGIONS_TXT = OUT_DIR / "windows_regions.txt"

W = 500  # половина ширины окна; окно = 2*W+1 bp (с обрезкой у концов хромосом)

# Длины хромосом из .fai (для обрезки правого края окна).
fai = {}
with open(str(IRGSP_FNA) + ".fai") as f:
    for line in f:
        c = line.split("\t")
        fai[c[0]] = int(c[1])

# Готовим регионы (1-based) в порядке bins.
entries = []  # (region_str, name, center_offset0)
with BIN_BED.open() as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 4:
            continue
        seqid, end0, name = c[0], int(c[2]), c[3]
        pos1 = end0  # 1-based позиция исходной точки
        if seqid not in fai:
            continue
        start1 = max(1, pos1 - W)
        end1 = min(fai[seqid], pos1 + W)
        region = f"{seqid}:{start1}-{end1}"
        center_offset0 = pos1 - start1  # 0-based смещение точки внутри окна
        entries.append((region, name, center_offset0))

REGIONS_TXT.write_text("".join(r + "\n" for r, _, _ in entries))

# Батч-извлечение окон одним вызовом samtools faidx -r.
proc = subprocess.run(
    ["samtools", "faidx", str(IRGSP_FNA), "-r", str(REGIONS_TXT)],
    check=True, capture_output=True, text=True,
)

# Парсим вывод samtools (заголовки = region-строки), сопоставляем по порядку.
records = []
name = None
seq = []
for line in proc.stdout.splitlines():
    if line.startswith(">"):
        if name is not None:
            records.append((name, "".join(seq)))
        name = line[1:].strip()
        seq = []
    else:
        seq.append(line.strip())
if name is not None:
    records.append((name, "".join(seq)))

assert len(records) == len(entries), f"records {len(records)} != regions {len(entries)}"

n_short = 0
with WINDOWS_FA.open("w") as fa, CENTER_BED.open("w") as bed:
    for (region, name, center_offset0), (hdr, sequence) in zip(entries, records):
        assert hdr.split()[0] == region, f"order mismatch: {hdr} != {region}"
        fa.write(f">{name}\n{sequence}\n")
        # центр (исходная точка) как 1bp интервал в координатах окна
        if center_offset0 >= len(sequence):
            center_offset0 = len(sequence) - 1
            n_short += 1
        bed.write(f"{name}\t{center_offset0}\t{center_offset0 + 1}\n")

print("Done.")
print(f"Windows FASTA: {WINDOWS_FA}  ({len(records)} windows, half-width {W} bp)")
print(f"Center BED:    {CENTER_BED}")
print(f"Windows clamped at chr end: {n_short}")

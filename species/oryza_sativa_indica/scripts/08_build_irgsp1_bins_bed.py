#!/usr/bin/env python3

# Подготовка BED-файла bin-серединок Yu 2011 в координатах IRGSP-1.0.
#
# Берём готовую japonica-карту (она же — те же bins Yu 2011, уже чисто перенесённые
# Nipponbare TIGR6.1 -> Nipponbare IRGSP-1.0, 100%). Для каждого bin есть позиция на
# IRGSP-1.0 (irgsp1_seqid / irgsp1_start0 / irgsp1_end0). Эти точки далее переносятся
# с IRGSP-1.0 на indica-родителя ZS97RS3.
#
# Имя в BED кодирует bin, исходную хромосому и cM, чтобы потом собрать карту.
# Запускать из корня репозитория.

from pathlib import Path
import csv

JAPONICA_DETAILS = Path("species/oryza_sativa_japonica/results/final/oryza_sativa_japonica_genetic_map.details.tsv")
OUT_BED = Path("species/oryza_sativa_indica/data/raw/yu2011_bins_irgsp1_midpoints.bed")
OUT_BED.parent.mkdir(parents=True, exist_ok=True)

n = 0
with JAPONICA_DETAILS.open() as fin, OUT_BED.open("w", newline="") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    for r in reader:
        seqid = r["irgsp1_seqid"]
        start0 = int(r["irgsp1_start0"])
        end0 = int(r["irgsp1_end0"])
        name = f'{r["bin"]}|src={r["source_chr"]}|cM={r["cM"]}'
        fout.write(f"{seqid}\t{start0}\t{end0}\t{name}\n")
        n += 1

print("Done.")
print(f"BED: {OUT_BED}")
print(f"Bins written: {n}")

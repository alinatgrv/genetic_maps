#!/usr/bin/env python3

# Сборка indica-карты на ZS97RS3 из перенесённых bin-серединок Yu 2011.
#
# Воронка (как у исходной indica-карты, но честный QC):
#   lifted -> primary chr 1..12 -> same_chr (target chr == source chr)
#   -> [QC: РЕАЛЬНАЯ коллинеарность Spearman ДО фильтра] -> LIS strict monotonic.
#
# Важно: "0 decreasing cM steps" в strict-карте навязано LIS-фильтром, поэтому
# отдельно сохраняем коллинеарность same-chr НАБОРА (prefilter) — это честная мера.
#
# Запускать из корня репозитория.

from pathlib import Path
from bisect import bisect_left
from collections import defaultdict, Counter
import csv
import re

SPECIES_DIR = Path("species/oryza_sativa_indica")
LIFTED = SPECIES_DIR / "results/liftover/yu2011_irgsp1_to_zs97rs3/yu2011_bins_zs97rs3_lifted.bed"

OUT_CAND = SPECIES_DIR / "results/final/oryza_sativa_indica_genetic_map.zs97rs3_projection.candidate.tsv"
OUT_CAND_DET = SPECIES_DIR / "results/final/oryza_sativa_indica_genetic_map.zs97rs3_projection.candidate.details.tsv"
OUT_STRICT = SPECIES_DIR / "results/final/oryza_sativa_indica_genetic_map.zs97rs3_projection.strict_monotonic.tsv"
OUT_STRICT_DET = SPECIES_DIR / "results/final/oryza_sativa_indica_genetic_map.zs97rs3_projection.strict_monotonic.details.tsv"
OUT_STD = SPECIES_DIR / "results/final/oryza_sativa_indica_genetic_map.zs97rs3.tsv"

QC_SUMMARY = SPECIES_DIR / "results/qc/yu2011_zs97rs3_projection_summary.txt"
QC_PREFILTER = SPECIES_DIR / "results/qc/yu2011_zs97rs3_collinearity_prefilter.tsv"
QC_BY_CHR = SPECIES_DIR / "results/qc/yu2011_zs97rs3_strict_monotonic_by_chr.tsv"
QC_EXCLUDED = SPECIES_DIR / "results/qc/yu2011_zs97rs3_strict_monotonic_excluded_bins.tsv"

for p in [OUT_CAND, QC_SUMMARY]:
    p.parent.mkdir(parents=True, exist_ok=True)
QC_SUMMARY.parent.mkdir(parents=True, exist_ok=True)

# ZS97RS3 accession -> номер хромосомы.
ACC_TO_CHR = {f"CP0560{52+i}.1": str(i + 1) for i in range(12)}
CHR_ORDER = [str(i) for i in range(1, 13)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_ORDER)}

# paftools.js liftover дописывает к имени суффикс "_<start>_<end>" (интервал центра),
# поэтому он опционален в конце.
NAME_RE = re.compile(r"^(?P<bin>[^|]+)\|src=(?P<src>\d+)\|cM=(?P<cM>[-0-9.]+)(?:_\d+_\d+)?$")


def lis_indices_strict(values):
    """Индексы одной самой длинной строго возрастающей подпоследовательности."""
    if not values:
        return []
    tails, tails_idx, prev = [], [], [-1] * len(values)
    for i, v in enumerate(values):
        j = bisect_left(tails, v)
        if j == len(tails):
            tails.append(v); tails_idx.append(i)
        else:
            tails[j] = v; tails_idx[j] = i
        if j > 0:
            prev[i] = tails_idx[j - 1]
    keep, k = [], tails_idx[-1]
    while k != -1:
        keep.append(k); k = prev[k]
    return list(reversed(keep))


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    def rank(a):
        order = sorted(range(n), key=lambda i: a[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and a[order[j + 1]] == a[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    if dx == 0 or dy == 0:
        return 0.0
    return num / (dx * dy)


# --- Чтение lifted BED ---
raw = []
nonprimary = 0
parse_fail = 0
with LIFTED.open() as f:
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 4:
            continue
        seqid, start0, end0, name = c[0], c[1], c[2], c[3]
        m = NAME_RE.match(name)
        if not m:
            parse_fail += 1
            continue
        chrom = ACC_TO_CHR.get(seqid)
        if chrom is None:
            nonprimary += 1
            continue
        pos = int(end0)  # 1-based позиция середины bin
        raw.append({
            "bin": m.group("bin"), "source_chr": m.group("src"), "cM": m.group("cM"),
            "chr": chrom, "pos": pos, "target_seqid": seqid,
        })

# Уникализация по bin: если bin перенёсся в несколько мест — предпочитаем хит на
# хромосоме, совпадающей с исходной.
by_bin = defaultdict(list)
for r in raw:
    by_bin[r["bin"]].append(r)

candidate = []
multi_hit_bins = 0
for b, rs in by_bin.items():
    if len(rs) > 1:
        multi_hit_bins += 1
        same = [r for r in rs if r["chr"] == r["source_chr"]]
        chosen = sorted(same or rs, key=lambda r: (CHR_INDEX[r["chr"]], r["pos"]))[0]
    else:
        chosen = rs[0]
    candidate.append(chosen)

candidate.sort(key=lambda r: (CHR_INDEX[r["chr"]], r["pos"]))
same_chr = [r for r in candidate if r["chr"] == r["source_chr"]]

# --- Честный QC: коллинеарность same-chr набора ДО монотонного фильтра ---
prefilter_rows = []
by_chr_same = defaultdict(list)
for r in same_chr:
    by_chr_same[r["chr"]].append(r)

for chrom in CHR_ORDER:
    rs = sorted(by_chr_same.get(chrom, []), key=lambda r: (r["pos"], float(r["cM"])))
    n = len(rs)
    if n == 0:
        continue
    cms = [float(r["cM"]) for r in rs]
    poss = [r["pos"] for r in rs]
    sp = spearman(poss, cms)
    diffs = [cms[i + 1] - cms[i] for i in range(n - 1)]
    inc = sum(1 for d in diffs if d > 0)
    dec = sum(1 for d in diffs if d < 0)
    eq = sum(1 for d in diffs if d == 0)
    prefilter_rows.append({
        "chr": chrom, "n_bins": n,
        "spearman_pos_cM": f"{sp:.4f}" if sp is not None else "NA",
        "increasing_cM_steps": inc, "decreasing_cM_steps": dec, "equal_cM_steps": eq,
    })

# --- LIS strict monotonic ---
kept, excluded = [], []
for chrom in CHR_ORDER:
    rs = sorted(by_chr_same.get(chrom, []), key=lambda r: (r["pos"], float(r["cM"]), r["bin"]))
    keep_idx = set(lis_indices_strict([float(r["cM"]) for r in rs]))
    for i, r in enumerate(rs):
        (kept if i in keep_idx else excluded).append(r)

kept.sort(key=lambda r: (CHR_INDEX[r["chr"]], r["pos"]))


def write_map(path, rows):
    with path.open("w", newline="") as w:
        w.write("chr\tpos\tcM\n")
        for r in rows:
            w.write(f'{r["chr"]}\t{r["pos"]}\t{r["cM"]}\n')


def write_details(path, rows):
    cols = ["chr", "pos", "cM", "bin", "source_chr", "target_seqid"]
    with path.open("w", newline="") as w:
        w.write("\t".join(cols) + "\n")
        for r in rows:
            w.write("\t".join(str(r[c]) for c in cols) + "\n")


write_map(OUT_CAND, candidate)
write_details(OUT_CAND_DET, candidate)
write_map(OUT_STRICT, kept)
write_details(OUT_STRICT_DET, kept)
write_map(OUT_STD, kept)

with QC_PREFILTER.open("w", newline="") as w:
    cols = ["chr", "n_bins", "spearman_pos_cM", "increasing_cM_steps", "decreasing_cM_steps", "equal_cM_steps"]
    w.write("\t".join(cols) + "\n")
    for r in prefilter_rows:
        w.write("\t".join(str(r[c]) for c in cols) + "\n")

final_by_chr = defaultdict(list)
for r in kept:
    final_by_chr[r["chr"]].append(r)
with QC_BY_CHR.open("w", newline="") as w:
    w.write("chr\tn_bins\tpos_min\tpos_max\tcM_min\tcM_max\n")
    for chrom in CHR_ORDER:
        rs = final_by_chr.get(chrom, [])
        if not rs:
            continue
        cms = [float(r["cM"]) for r in rs]
        w.write(f'{chrom}\t{len(rs)}\t{min(r["pos"] for r in rs)}\t{max(r["pos"] for r in rs)}\t{min(cms):.6f}\t{max(cms):.6f}\n')

with QC_EXCLUDED.open("w", newline="") as w:
    cols = ["chr", "pos", "cM", "bin", "source_chr", "target_seqid"]
    w.write("\t".join(cols) + "\texclude_reason\n")
    for r in sorted(excluded, key=lambda r: (CHR_INDEX[r["chr"]], r["pos"])):
        w.write("\t".join(str(r[c]) for c in cols) + "\tbreaks_monotonic_cM_order_after_liftover\n")

# Честные агрегаты коллинеарности (по same-chr, ДО фильтра).
abs_sp = [abs(float(r["spearman_pos_cM"])) for r in prefilter_rows if r["spearman_pos_cM"] != "NA"]
mean_abs = sum(abs_sp) / len(abs_sp) if abs_sp else float("nan")
total_dec_prefilter = sum(r["decreasing_cM_steps"] for r in prefilter_rows)

with QC_SUMMARY.open("w") as w:
    w.write("Oryza sativa indica — Yu2011 projection onto ZS97RS3 (parent genome)\n")
    w.write("=" * 66 + "\n\n")
    w.write(f"input_bins_irgsp1\t1619\n")
    w.write(f"lifted_rows_raw\t{len(raw)}\n")
    w.write(f"bins_lifted_unique\t{len(by_bin)}\n")
    w.write(f"bins_multi_hit\t{multi_hit_bins}\n")
    w.write(f"nonprimary_lifted_rows\t{nonprimary}\n")
    w.write(f"candidate_primary_1_12\t{len(candidate)}\n")
    w.write(f"same_chr_as_source\t{len(same_chr)}\n")
    w.write(f"off_chr\t{len(candidate) - len(same_chr)}\n")
    w.write(f"strict_monotonic_kept\t{len(kept)}\n")
    w.write(f"strict_monotonic_excluded\t{len(excluded)}\n\n")
    w.write(f"lift_rate_pct\t{100.0 * len(by_bin) / 1619:.1f}\n")
    w.write(f"same_chr_pct_of_candidate\t{100.0 * len(same_chr) / max(len(candidate),1):.1f}\n\n")
    w.write("HONEST collinearity (same-chr set, BEFORE monotonic filter):\n")
    w.write(f"  mean_abs_spearman_over_chr\t{mean_abs:.4f}\n")
    w.write(f"  total_decreasing_cM_steps_prefilter\t{total_dec_prefilter}\n")
    w.write("  NOTE: strict_monotonic map has 0 decreasing steps BY CONSTRUCTION (LIS filter).\n")

print("Done.")
print(f"candidate: {len(candidate)} | same_chr: {len(same_chr)} | strict_monotonic: {len(kept)}")
print(f"lift rate: {100.0*len(by_bin)/1619:.1f}%  | mean |Spearman| prefilter: {mean_abs:.4f}")
print(f"Std map: {OUT_STD}")

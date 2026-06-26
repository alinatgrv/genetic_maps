#!/usr/bin/env python3

# Построение генетической карты рапса (chr / pos / cM) из BLAST-выравнивания
# DArT-клонов на референс Da-Ae (GCF_020379485.1).
#
# Логика (как у подсолнечника/капусты — уникальный хит на ожидаемой хромосоме):
#   * берём только хиты на 19 хромосом (NC_063434.1 .. NC_063452.1);
#   * для каждого маркера лучший HSP по bitscore определяет позицию
#     (pos = середина интервала на референсе);
#   * маркер уникален, если нет конкурирующего хита в ДРУГОМ локусе с
#     bitscore >= 0.95 * лучшего (так отсекаются гомеологичные копии A<->C);
#   * в строгую карту берём маркеры, у которых лучший хит уникален и лежит
#     на хромосоме, совпадающей с группой сцепления из консенсус-карты.
#
# Также пишем "relaxed" карту (лучший хит на ожидаемой хромосоме без требования
# уникальности) и промежуточные таблицы/QC.
#
# Запускать из корня репозитория.

from pathlib import Path
from collections import defaultdict, Counter

SPECIES_DIR = Path("species/brassica_napus")

hits_file = SPECIES_DIR / "results/intermediate/bnapus_dart_markers.blastn.tsv"
metadata_file = SPECIES_DIR / "data/metadata/bnapus_consensus_map_metadata.tsv"

out_best = SPECIES_DIR / "results/intermediate/bnapus_marker_best_hits.tsv"
out_final = SPECIES_DIR / "results/final/brassica_napus_genetic_map.tsv"
out_final_markers = SPECIES_DIR / "results/final/brassica_napus_genetic_map.with_markers.tsv"
out_relaxed = SPECIES_DIR / "results/final/brassica_napus_genetic_map.relaxed_best_hit.tsv"
out_excluded = SPECIES_DIR / "results/qc/bnapus_excluded_markers.tsv"
out_summary = SPECIES_DIR / "results/qc/bnapus_map_build_summary.txt"

for p in [out_best, out_final, out_relaxed, out_excluded, out_summary]:
    p.parent.mkdir(parents=True, exist_ok=True)

# Accession -> имя хромосомы (Da-Ae).
ACC_TO_CHR = {
    "NC_063434.1": "A1", "NC_063435.1": "A2", "NC_063436.1": "A3",
    "NC_063437.1": "A4", "NC_063438.1": "A5", "NC_063439.1": "A6",
    "NC_063440.1": "A7", "NC_063441.1": "A8", "NC_063442.1": "A9",
    "NC_063443.1": "A10",
    "NC_063444.1": "C1", "NC_063445.1": "C2", "NC_063446.1": "C3",
    "NC_063447.1": "C4", "NC_063448.1": "C5", "NC_063449.1": "C6",
    "NC_063450.1": "C7", "NC_063451.1": "C8", "NC_063452.1": "C9",
}

CHR_ORDER = [f"A{i}" for i in range(1, 11)] + [f"C{i}" for i in range(1, 10)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_ORDER)}

UNIQUE_BITSCORE_FRAC = 0.95   # конкурент с bitscore >= 0.95*best => не уникален
SAME_LOCUS_WINDOW = 2000      # хиты ближе этого на той же последовательности — один локус


# --- Метаданные карты: join_id -> {chr: cM} ---
meta = defaultdict(dict)
with open(metadata_file, encoding="utf-8") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    for line in handle:
        f = line.rstrip("\n").split("\t")
        if f[idx["is_dart"]] != "1":
            continue
        meta[f[idx["marker_join_id"]]][f[idx["chr"]]] = f[idx["cM"]]


# --- HSP по маркерам (только хромосомные) ---
COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]

hits_by_marker = defaultdict(list)
total_hsp = 0
chromosomal_hsp = 0
with open(hits_file, encoding="utf-8") as handle:
    for line in handle:
        f = line.rstrip("\n").split("\t")
        if len(f) < len(COLS):
            continue
        total_hsp += 1
        rec = dict(zip(COLS, f))
        chrom = ACC_TO_CHR.get(rec["sseqid"])
        if chrom is None:
            continue  # неразмещённый скэффолд NW_* — пропускаем
        chromosomal_hsp += 1
        rec["chr"] = chrom
        rec["pident"] = float(rec["pident"])
        rec["length"] = int(rec["length"])
        rec["sstart"] = int(rec["sstart"])
        rec["send"] = int(rec["send"])
        rec["bitscore"] = float(rec["bitscore"])
        rec["pos"] = (rec["sstart"] + rec["send"]) // 2
        hits_by_marker[rec["qseqid"]].append(rec)


def is_different_locus(a, b):
    if a["sseqid"] != b["sseqid"]:
        return True
    return abs(a["sstart"] - b["sstart"]) > SAME_LOCUS_WINDOW


best_rows = []
final_markers = []          # строгая карта
relaxed_markers = []        # лучший хит на ожидаемой хромосоме
excluded = []               # причины исключения

for marker, hsps in hits_by_marker.items():
    hsps.sort(key=lambda r: (r["bitscore"], r["pident"], r["length"]), reverse=True)
    best = hsps[0]

    # Лучший конкурент в другом локусе.
    competitor = None
    for h in hsps[1:]:
        if is_different_locus(h, best):
            competitor = h
            break
    unique = competitor is None or competitor["bitscore"] < UNIQUE_BITSCORE_FRAC * best["bitscore"]

    expected = meta.get(marker, {})
    expected_chrs = set(expected.keys())
    chr_agree = best["chr"] in expected_chrs
    cM = expected.get(best["chr"]) or (next(iter(expected.values())) if expected else "")

    best_rows.append({
        "marker_id": marker,
        "best_chr": best["chr"],
        "pos": str(best["pos"]),
        "cM": cM,
        "expected_chr": ",".join(sorted(expected_chrs, key=lambda c: CHR_INDEX.get(c, 99))),
        "pident": f'{best["pident"]:.2f}',
        "aln_len": str(best["length"]),
        "bitscore": f'{best["bitscore"]:.1f}',
        "competitor_chr": competitor["chr"] if competitor else "",
        "competitor_bitscore": f'{competitor["bitscore"]:.1f}' if competitor else "",
        "unique": "1" if unique else "0",
        "chr_agree": "1" if chr_agree else "0",
    })

    if chr_agree:
        relaxed_markers.append((best["chr"], best["pos"], cM, marker, best))
    if unique and chr_agree:
        final_markers.append((best["chr"], best["pos"], cM, marker, best))
    else:
        reason = []
        if not chr_agree:
            reason.append("best_hit_off_expected_chr")
        if not unique:
            reason.append("non_unique_homoeolog")
        excluded.append({
            "marker_id": marker,
            "best_chr": best["chr"],
            "expected_chr": ",".join(sorted(expected_chrs, key=lambda c: CHR_INDEX.get(c, 99))),
            "pos": str(best["pos"]),
            "pident": f'{best["pident"]:.2f}',
            "bitscore": f'{best["bitscore"]:.1f}',
            "competitor_chr": competitor["chr"] if competitor else "",
            "competitor_bitscore": f'{competitor["bitscore"]:.1f}' if competitor else "",
            "reason": ";".join(reason),
        })

# Маркеры без единого хромосомного хита.
no_hit = sorted(set(meta) - set(hits_by_marker))
for marker in no_hit:
    excluded.append({
        "marker_id": marker,
        "best_chr": "", "expected_chr": ",".join(sorted(meta[marker], key=lambda c: CHR_INDEX.get(c, 99))),
        "pos": "", "pident": "", "bitscore": "", "competitor_chr": "", "competitor_bitscore": "",
        "reason": "no_chromosomal_hit",
    })


def sort_key(t):
    return (CHR_INDEX[t[0]], t[1], float(t[2]) if t[2] != "" else 0.0, t[3])


def dedup_positions(markers):
    """Убираем дубли (chr,pos,cM) и позиции (chr,pos) с конфликтом cM."""
    markers = sorted(markers, key=sort_key)
    seen = set()
    dedup = []
    for t in markers:
        key = (t[0], t[1], t[2])
        if key not in seen:
            seen.add(key)
            dedup.append(t)
    cm_by_pos = defaultdict(set)
    for t in dedup:
        cm_by_pos[(t[0], t[1])].add(t[2])
    conflict = {k for k, v in cm_by_pos.items() if len(v) > 1}
    clean = [t for t in dedup if (t[0], t[1]) not in conflict]
    return clean, len(conflict)


final_clean, final_conflicts = dedup_positions(final_markers)
relaxed_clean, relaxed_conflicts = dedup_positions(relaxed_markers)


def write_tsv(path, header, records):
    with open(path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for rec in records:
            out.write("\t".join(str(rec[h]) for h in header) + "\n")


best_rows.sort(key=lambda r: (CHR_INDEX.get(r["best_chr"], 99), int(r["pos"])))
write_tsv(out_best,
          ["marker_id", "best_chr", "pos", "cM", "expected_chr", "pident", "aln_len",
           "bitscore", "competitor_chr", "competitor_bitscore", "unique", "chr_agree"],
          best_rows)

# Финальная строгая карта chr/pos/cM.
with open(out_final, "w", encoding="utf-8") as out:
    out.write("chr\tpos\tcM\n")
    for chrom, pos, cM, marker, best in final_clean:
        out.write(f"{chrom}\t{pos}\t{cM}\n")

# Та же карта с маркерами и метриками выравнивания.
with open(out_final_markers, "w", encoding="utf-8") as out:
    out.write("chr\tpos\tcM\tmarker_id\tpident\taln_len\tbitscore\n")
    for chrom, pos, cM, marker, best in final_clean:
        out.write(f'{chrom}\t{pos}\t{cM}\t{marker}\t{best["pident"]:.2f}\t{best["length"]}\t{best["bitscore"]:.1f}\n')

# Relaxed-карта.
with open(out_relaxed, "w", encoding="utf-8") as out:
    out.write("chr\tpos\tcM\n")
    for chrom, pos, cM, marker, best in relaxed_clean:
        out.write(f"{chrom}\t{pos}\t{cM}\n")

excluded.sort(key=lambda r: (r["reason"], r["marker_id"]))
write_tsv(out_excluded,
          ["marker_id", "best_chr", "expected_chr", "pos", "pident", "bitscore",
           "competitor_chr", "competitor_bitscore", "reason"],
          excluded)

# Сводка.
final_by_chr = Counter(t[0] for t in final_clean)
relaxed_by_chr = Counter(t[0] for t in relaxed_clean)
reason_counts = Counter(r["reason"] for r in excluded)

with open(out_summary, "w", encoding="utf-8") as out:
    out.write("Brassica napus genetic map build summary\n")
    out.write("=" * 45 + "\n\n")
    out.write(f"DArT markers in consensus map (target set): {len(meta)}\n")
    out.write(f"Markers with >=1 chromosomal hit: {len(hits_by_marker)}\n")
    out.write(f"Total HSP rows: {total_hsp} (chromosomal: {chromosomal_hsp})\n\n")
    out.write(f"Strict map markers (unique + chr agree, deduplicated): {len(final_clean)}\n")
    out.write(f"  positions dropped by cM conflict: {final_conflicts}\n")
    out.write(f"Relaxed map markers (best hit on expected chr): {len(relaxed_clean)}\n")
    out.write(f"  positions dropped by cM conflict: {relaxed_conflicts}\n\n")
    out.write("Excluded markers by reason:\n")
    for reason, n in reason_counts.most_common():
        out.write(f"  {reason}\t{n}\n")
    out.write("\nStrict map markers by chromosome (strict / relaxed):\n")
    for ch in CHR_ORDER:
        out.write(f"  {ch}\t{final_by_chr.get(ch, 0)}\t{relaxed_by_chr.get(ch, 0)}\n")

print("Done.")
print(f"Best hits:        {out_best}")
print(f"Final strict map: {out_final}  ({len(final_clean)} markers)")
print(f"Relaxed map:      {out_relaxed}  ({len(relaxed_clean)} markers)")
print(f"Excluded:         {out_excluded}")
print(f"Summary:          {out_summary}")

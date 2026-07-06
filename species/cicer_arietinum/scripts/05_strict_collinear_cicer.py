#!/usr/bin/env python3
"""Turn the chickpea candidate map into a strictly collinear final map.

Per chromosome:
  1) orient: pick forward vs reverse by whichever keeps more markers (LNDS);
     if reverse, flip the genetic axis cM' = cM_max(chr) - cM so cM increases
     with physical position (consistent orientation across chromosomes);
  2) collapse duplicate physical positions to one representative (median cM);
  3) keep the Longest Non-Decreasing Subsequence of cM ordered by pos
     (= minimum markers removed for strict collinearity).

Input : results/intermediate/cicer_genetic_map.candidate.raw.tsv
Output: results/final/cicer_arietinum_genetic_map.candidate.tsv      (all placed)
        results/final/cicer_arietinum_genetic_map.tsv                (strict, canonical)
        results/final/cicer_arietinum_genetic_map.with_markers.tsv
Run from repo root in the genetic_maps env.
"""
import bisect
import csv
from pathlib import Path
from decimal import Decimal
from collections import defaultdict

SPECIES_DIR = Path("species/cicer_arietinum")
RAW = SPECIES_DIR / "results/intermediate/cicer_genetic_map.candidate.raw.tsv"
FINAL = SPECIES_DIR / "results/final"


def lnds_keep(vals):
    tails, tails_idx, prev = [], [], [-1] * len(vals)
    for i, x in enumerate(vals):
        j = bisect.bisect_right(tails, x)
        if j == len(tails):
            tails.append(x); tails_idx.append(i)
        else:
            tails[j] = x; tails_idx[j] = i
        prev[i] = tails_idx[j - 1] if j > 0 else -1
    keep, k = [], (tails_idx[-1] if tails_idx else -1)
    while k != -1:
        keep.append(k); k = prev[k]
    return set(keep)


def lnds_len(vals):
    tails = []
    for x in vals:
        j = bisect.bisect_right(tails, x)
        if j == len(tails):
            tails.append(x)
        else:
            tails[j] = x
    return len(tails)


def median(vals):
    s = sorted(vals)
    return s[(len(s) - 1) // 2]


def load():
    by_chr = defaultdict(list)
    with open(RAW, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_chr[int(r["chr"])].append((int(r["pos"]), Decimal(r["cM"]), r["marker_id"]))
    return by_chr


def strict(by_chr):
    out, out_full, stats = [], [], []
    for c in sorted(by_chr):
        rows = by_chr[c]
        n_in = len(rows)
        cm_max = max(m for _, m, _ in rows)
        fwd = sorted(rows, key=lambda t: (t[0], t[1]))
        rev = sorted(rows, key=lambda t: (t[0], -t[1]))
        flip = lnds_len([-m for _, m, _ in rev]) > lnds_len([m for _, m, _ in fwd])
        if flip:
            rows = [(p, cm_max - m, mid) for p, m, mid in rows]
        by_pos = defaultdict(list)
        for p, m, mid in rows:
            by_pos[p].append((m, mid))
        uniq, collapsed = [], 0
        for p, items in by_pos.items():
            if len(items) > 1:
                collapsed += len(items) - 1
            med = median([m for m, _ in items])
            mid = next((mid for m, mid in items if m == med), items[0][1])
            uniq.append((p, med, mid))
        uniq.sort(key=lambda t: (t[0], t[1]))
        keep = lnds_keep([m for _, m, _ in uniq])
        kept = [uniq[i] for i in sorted(keep)]
        for p, m, mid in kept:
            out.append((c, p, m)); out_full.append((c, p, m, mid))
        stats.append((c, n_in, int(flip), collapsed, len(uniq) - len(kept), len(kept)))
    return out, out_full, stats


def write_map(rows, path, with_markers=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chr", "pos", "cM", "marker_id"] if with_markers else ["chr", "pos", "cM"])
        for row in rows:
            w.writerow([row[0], row[1], str(row[2])] + ([row[3]] if with_markers else []))


def main():
    by_chr = load()
    # candidate (all placed), sorted
    cand = sorted(((c, p, m) for c in by_chr for p, m, _ in by_chr[c]),
                  key=lambda t: (t[0], t[1], t[2]))
    write_map(cand, FINAL / "cicer_arietinum_genetic_map.candidate.tsv")

    out, out_full, stats = strict(by_chr)
    write_map(out, FINAL / "cicer_arietinum_genetic_map.tsv")
    write_map(out_full, FINAL / "cicer_arietinum_genetic_map.with_markers.tsv", with_markers=True)

    n_in = sum(s[1] for s in stats)
    n_out = sum(s[5] for s in stats)
    flips = [str(s[0]) for s in stats if s[2]]
    print(f"candidate (all placed): {len(cand)} markers")
    print(f"STRICT collinear:       {n_in} -> {n_out} markers "
          f"(flipped chr: {flips or '-'}; collapsed: {sum(s[3] for s in stats)}; "
          f"removed non-collinear: {sum(s[4] for s in stats)})")
    print(f"   {'chr':>3}{'in':>6}{'flip':>5}{'collapse':>9}{'removed':>8}{'out':>6}")
    for c, n_in_c, flip, collapsed, removed, n_out_c in stats:
        print(f"   {c:>3}{n_in_c:>6}{flip:>5}{collapsed:>9}{removed:>8}{n_out_c:>6}")


if __name__ == "__main__":
    main()

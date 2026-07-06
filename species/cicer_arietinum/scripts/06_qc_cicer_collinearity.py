#!/usr/bin/env python3
"""Collinearity QC for the chickpea maps.

For the candidate map (all placed markers) and the strict final map, compute
per-chromosome Spearman(pos, cM), orientation, adjacent cM-step counts and the
monotonic fraction. The candidate QC exposes how well each ASM33114v1
pseudomolecule agrees with the genetic map; the strict QC confirms the delivered
map is monotonic (monotonic_fraction = 1.0).

Outputs:
    results/qc/cicer_candidate_map_collinearity_qc.tsv
    results/qc/cicer_final_map_collinearity_qc.tsv
    results/qc/cicer_final_map_summary.txt
Run from repo root in the genetic_maps env.
"""
import csv
from pathlib import Path
from collections import defaultdict

SPECIES_DIR = Path("species/cicer_arietinum")
FINAL = SPECIES_DIR / "results/final"
QC = SPECIES_DIR / "results/qc"


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = sum((rx[i] - mx) ** 2 for i in range(n)) ** 0.5
    dy = sum((ry[i] - my) ** 2 for i in range(n)) ** 0.5
    return num / (dx * dy) if dx and dy else float("nan")


def load(path):
    by = defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by[int(r["chr"])].append((int(r["pos"]), float(r["cM"])))
    return by


def qc_rows(by):
    rows = []
    for c in sorted(by):
        recs = sorted(by[c], key=lambda t: t[0])
        pos = [p for p, _ in recs]
        cm = [m for _, m in recs]
        rho = spearman(pos, cm)
        orient = "increasing" if (rho == rho and rho >= 0) else "decreasing"
        seq = cm if orient == "increasing" else cm[::-1]
        inc = sum(1 for i in range(1, len(seq)) if seq[i] > seq[i - 1])
        dec = sum(1 for i in range(1, len(seq)) if seq[i] < seq[i - 1])
        eq = sum(1 for i in range(1, len(seq)) if seq[i] == seq[i - 1])
        steps = max(1, len(seq) - 1)
        rows.append({
            "chr": c, "n_markers": len(recs),
            "pos_min": min(pos), "pos_max": max(pos),
            "physical_span_Mb": round((max(pos) - min(pos)) / 1e6, 3),
            "cM_min": round(min(cm), 3), "cM_max": round(max(cm), 3),
            "genetic_span_cM": round(max(cm) - min(cm), 3),
            "spearman_pos_cM": round(rho, 4), "abs_spearman": round(abs(rho), 4),
            "orientation": orient,
            "increasing_cM_steps": inc, "decreasing_cM_steps": dec, "equal_cM_steps": eq,
            "monotonic_fraction": round((inc + eq) / steps, 3),
        })
    return rows


def write_tsv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    cand = load(FINAL / "cicer_arietinum_genetic_map.candidate.tsv")
    strict = load(FINAL / "cicer_arietinum_genetic_map.tsv")
    cand_rows = qc_rows(cand)
    strict_rows = qc_rows(strict)
    write_tsv(cand_rows, QC / "cicer_candidate_map_collinearity_qc.tsv")
    write_tsv(strict_rows, QC / "cicer_final_map_collinearity_qc.tsv")

    n_cand = sum(r["n_markers"] for r in cand_rows)
    n_strict = sum(r["n_markers"] for r in strict_rows)
    med_abs = sorted(r["abs_spearman"] for r in strict_rows)[len(strict_rows) // 2]
    all_mono = all(r["monotonic_fraction"] == 1.0 for r in strict_rows)
    lines = [
        "Chickpea (Cicer arietinum) genetic map — collinearity QC summary",
        "=" * 62,
        "Reference: GCF_000331145.1 (ASM33114v1, CDC Frontier), chromosomes Ca1-Ca8",
        "Source:    Gaur et al. 2015, Sci Rep 5:13387 (SNP flanks + interspecific cM)",
        "",
        f"Candidate map (all uniquely placed SNPs): {n_cand} markers",
        f"Strict collinear final map:               {n_strict} markers",
        f"Strict map: every chromosome monotonic (monotonic_fraction==1.0): {all_mono}",
        f"Strict map: median |Spearman| over chromosomes: {med_abs}",
        "",
        "Candidate map per-chromosome |Spearman| (assembly-vs-map agreement):",
    ]
    for r in cand_rows:
        lines.append(f"  Ca{r['chr']}: n={r['n_markers']:4d}  spearman={r['spearman_pos_cM']:+.3f}  "
                     f"monotonic_fraction={r['monotonic_fraction']}")
    lines += ["", "Strict final map markers per chromosome:"]
    for r in strict_rows:
        lines.append(f"  Ca{r['chr']}: {r['n_markers']}")
    (QC / "cicer_final_map_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

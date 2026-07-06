#!/usr/bin/env python3

# Визуализация финальной генетической карты нута (Cicer arietinum):
#   * покрытие маркерами по cM (genetic map coverage);
#   * покрытие маркерами по физической позиции (physical map coverage);
#   * плотность маркеров в бинах по 10 cM.
#
# Плюс сводные таблицы по хромосомам и по плотности. SVG генерируется без внешних
# зависимостей (как в остальных культурах проекта). Запускать из корня репозитория.

from pathlib import Path
import csv
import math
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
SPECIES_DIR = SCRIPT_DIR.parent

MAPS = [
    {
        "name": "strict",
        "label": "Gaur 2015 SNP map, strict collinear (final)",
        "path": SPECIES_DIR / "results/final/cicer_arietinum_genetic_map.tsv",
    },
]

FIG_DIR = SPECIES_DIR / "results/figures"
QC_DIR = SPECIES_DIR / "results/qc"
FIG_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

CHR_ORDER = [str(i) for i in range(1, 9)]
CHR_INDEX = {c: i for i, c in enumerate(CHR_ORDER)}


def load_map(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({"chr": str(row["chr"]), "pos": float(row["pos"]), "cM": float(row["cM"])})
    return rows


def chr_sort_key(chrom):
    return CHR_INDEX.get(chrom, 99)


def svg_escape(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_svg(path, width, height, body):
    with open(path, "w", encoding="utf-8") as out:
        out.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n')
        out.write('<rect width="100%" height="100%" fill="white"/>\n')
        out.write(body)
        out.write("\n</svg>\n")


def make_chromosome_summary(all_data):
    rows = []
    for map_info in MAPS:
        by_chr = defaultdict(list)
        for row in all_data[map_info["name"]]:
            by_chr[row["chr"]].append(row)
        for chrom in sorted(by_chr, key=chr_sort_key):
            chrom_rows = by_chr[chrom]
            pos_values = [r["pos"] for r in chrom_rows]
            cm_values = [r["cM"] for r in chrom_rows]
            physical_span_bp = max(pos_values) - min(pos_values)
            genetic_span_cm = max(cm_values) - min(cm_values)
            rows.append({
                "map_name": map_info["name"], "chr": chrom, "n_markers": len(chrom_rows),
                "pos_min": int(min(pos_values)), "pos_max": int(max(pos_values)),
                "physical_span_bp": int(physical_span_bp),
                "cm_min": min(cm_values), "cm_max": max(cm_values),
                "genetic_span_cM": genetic_span_cm,
                "markers_per_Mbp": round(len(chrom_rows) / (physical_span_bp / 1_000_000), 3) if physical_span_bp > 0 else "",
                "markers_per_cM": round(len(chrom_rows) / genetic_span_cm, 3) if genetic_span_cm > 0 else "",
            })
    out_tsv = QC_DIR / "cicer_final_map_chromosome_summary.tsv"
    with open(out_tsv, "w", encoding="utf-8", newline="") as handle:
        fieldnames = ["map_name", "chr", "n_markers", "pos_min", "pos_max", "physical_span_bp",
                      "cm_min", "cm_max", "genetic_span_cM", "markers_per_Mbp", "markers_per_cM"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return out_tsv


def make_rug_map_svg(all_data, value_col, unit_label, out_path, title):
    chroms = sorted({row["chr"] for data in all_data.values() for row in data}, key=chr_sort_key)
    width = 1500
    panel_height = 560
    top_margin, left_margin, right_margin, bottom_margin, panel_gap = 80, 140, 90, 80, 60
    height = top_margin + len(MAPS) * panel_height + (len(MAPS) - 1) * panel_gap + bottom_margin
    plot_width = width - left_margin - right_margin
    chr_step = plot_width / len(chroms)

    all_values = [row[value_col] for data in all_data.values() for row in data]
    global_min = 0
    global_max = (math.ceil(max(all_values) / 10_000_000) * 10_000_000
                  if value_col == "pos" else math.ceil(max(all_values) / 10) * 10)
    if global_max == global_min:
        global_max = global_min + 1

    def y_scale(value, panel_top):
        plot_top = panel_top + 60
        plot_bottom = panel_top + panel_height - 80
        return plot_top + (value - global_min) / (global_max - global_min) * (plot_bottom - plot_top)

    body = [f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">{svg_escape(title)}</text>\n']
    for panel_i, map_info in enumerate(MAPS):
        data = all_data[map_info["name"]]
        panel_top = top_margin + panel_i * (panel_height + panel_gap)
        plot_top = panel_top + 60
        plot_bottom = panel_top + panel_height - 80
        body.append(f'<text x="{left_margin}" y="{panel_top + 25}" font-family="Arial" font-size="20" font-weight="bold">{svg_escape(map_info["label"])}</text>\n')
        body.append(f'<line x1="{left_margin-25}" y1="{plot_top}" x2="{left_margin-25}" y2="{plot_bottom}" stroke="black" stroke-width="1"/>\n')
        for t in range(6):
            value = global_min + (global_max - global_min) * t / 5
            y = y_scale(value, panel_top)
            label = f"{value/1_000_000:.0f}" if value_col == "pos" else f"{value:.0f}"
            body.append(f'<line x1="{left_margin-30}" y1="{y}" x2="{left_margin-25}" y2="{y}" stroke="black" stroke-width="1"/>\n')
            body.append(f'<text x="{left_margin-38}" y="{y+4}" text-anchor="end" font-family="Arial" font-size="12">{label}</text>\n')
        body.append(f'<text x="{left_margin-88}" y="{(plot_top+plot_bottom)/2}" transform="rotate(-90 {left_margin-88},{(plot_top+plot_bottom)/2})" text-anchor="middle" font-family="Arial" font-size="14">{svg_escape(unit_label)}</text>\n')
        by_chr = defaultdict(list)
        for row in data:
            by_chr[row["chr"]].append(row)
        for i, chrom in enumerate(chroms):
            x = left_margin + chr_step * (i + 0.5)
            vals = [r[value_col] for r in by_chr[chrom]]
            y1 = y_scale(min(vals), panel_top) if vals else plot_top
            y2 = y_scale(max(vals), panel_top) if vals else plot_bottom
            body.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" stroke="black" stroke-width="2"/>\n')
            for r in sorted(by_chr[chrom], key=lambda z: z[value_col]):
                y = y_scale(r[value_col], panel_top)
                body.append(f'<line x1="{x-10}" y1="{y}" x2="{x+10}" y2="{y}" stroke="black" stroke-width="0.55"/>\n')
            body.append(f'<text x="{x}" y="{plot_bottom + 28}" text-anchor="middle" font-family="Arial" font-size="14">Ca{svg_escape(chrom)}</text>\n')
            body.append(f'<text x="{x}" y="{plot_bottom + 48}" text-anchor="middle" font-family="Arial" font-size="12">{len(by_chr[chrom])}</text>\n')
        body.append(f'<text x="{width - right_margin}" y="{plot_bottom + 48}" text-anchor="end" font-family="Arial" font-size="12">numbers = markers per chromosome</text>\n')
    write_svg(out_path, width, height, "".join(body))


def make_density_10cm(all_data):
    bin_size = 10.0
    density_rows = []
    for map_info in MAPS:
        by_chr = defaultdict(list)
        for row in all_data[map_info["name"]]:
            by_chr[row["chr"]].append(row["cM"])
        for chrom, values in sorted(by_chr.items(), key=lambda x: chr_sort_key(x[0])):
            if not values:
                continue
            cmin = math.floor(min(values) / bin_size) * bin_size
            cmax = math.ceil(max(values) / bin_size) * bin_size
            current = cmin
            while current < cmax:
                count = sum(1 for v in values if current <= v < current + bin_size or (v == cmax and current + bin_size == cmax))
                density_rows.append({"map_name": map_info["name"], "chr": chrom,
                                     "bin_start_cM": current, "bin_end_cM": current + bin_size, "n_markers": count})
                current += bin_size
    out_tsv = QC_DIR / "cicer_marker_density_10cM_bins.tsv"
    with open(out_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["map_name", "chr", "bin_start_cM", "bin_end_cM", "n_markers"], delimiter="\t")
        writer.writeheader()
        writer.writerows(density_rows)
    make_density_svg(density_rows, FIG_DIR / "cicer_marker_density_10cM.svg")
    return out_tsv


def make_density_svg(density_rows, out_path):
    chroms = sorted({r["chr"] for r in density_rows}, key=chr_sort_key)
    max_bin = max(r["bin_end_cM"] for r in density_rows)
    max_count = max(r["n_markers"] for r in density_rows) if density_rows else 1
    bin_size = 10
    bins = list(range(0, int(math.ceil(max_bin / bin_size) * bin_size), bin_size))
    cell_w, cell_h, left_margin, top_margin, panel_gap, title_h = 70, 26, 150, 85, 80, 42
    width = left_margin + len(chroms) * cell_w + 190
    panel_height = title_h + len(bins) * cell_h + 55
    height = top_margin + len(MAPS) * panel_height + (len(MAPS)-1) * panel_gap + 85
    body = [f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">Chickpea (Cicer arietinum) final genetic map: marker density in 10 cM bins</text>\n']
    by_key = {(r["map_name"], r["chr"], int(r["bin_start_cM"])): r["n_markers"] for r in density_rows}
    for panel_i, map_info in enumerate(MAPS):
        panel_top = top_margin + panel_i * (panel_height + panel_gap)
        body.append(f'<text x="{left_margin}" y="{panel_top}" font-family="Arial" font-size="20" font-weight="bold">{svg_escape(map_info["label"])}</text>\n')
        for i, chrom in enumerate(chroms):
            x = left_margin + i * cell_w
            body.append(f'<text x="{x + cell_w/2}" y="{panel_top + 30}" text-anchor="middle" font-family="Arial" font-size="13">Ca{svg_escape(chrom)}</text>\n')
        for j, bin_start in enumerate(bins):
            y = panel_top + title_h + j * cell_h
            body.append(f'<text x="{left_margin - 15}" y="{y + 17}" text-anchor="end" font-family="Arial" font-size="12">{bin_start}-{bin_start+bin_size}</text>\n')
            for i, chrom in enumerate(chroms):
                x = left_margin + i * cell_w
                count = by_key.get((map_info["name"], chrom, bin_start), 0)
                shade = max(35, min(255, 255 - int(200 * (count / max_count)) if max_count else 255))
                body.append(f'<rect x="{x}" y="{y}" width="{cell_w-4}" height="{cell_h-3}" fill="rgb({shade},{shade},{shade})" stroke="white"/>\n')
                if count > 0:
                    text_color = "white" if shade < 120 else "black"
                    body.append(f'<text x="{x + (cell_w-4)/2}" y="{y + 17}" text-anchor="middle" font-family="Arial" font-size="12" fill="{text_color}">{count}</text>\n')
        body.append(f'<text x="{left_margin - 86}" y="{panel_top + title_h + len(bins)*cell_h/2}" transform="rotate(-90 {left_margin-86},{panel_top + title_h + len(bins)*cell_h/2})" text-anchor="middle" font-family="Arial" font-size="14">cM interval</text>\n')
    legend_x, legend_y = width - 135, 85
    body.append(f'<text x="{legend_x}" y="{legend_y}" font-family="Arial" font-size="13">n markers</text>\n')
    for k in range(6):
        count = round(max_count * k / 5)
        shade = max(35, min(255, 255 - int(200 * (count / max_count)) if max_count else 255))
        body.append(f'<rect x="{legend_x}" y="{legend_y + 15 + k*24}" width="25" height="18" fill="rgb({shade},{shade},{shade})" stroke="black" stroke-width="0.5"/>\n')
        body.append(f'<text x="{legend_x + 35}" y="{legend_y + 29 + k*24}" font-family="Arial" font-size="12">{count}</text>\n')
    write_svg(out_path, width, height, "".join(body))


def main():
    all_data = {}
    for map_info in MAPS:
        if not map_info["path"].exists():
            raise FileNotFoundError(map_info["path"])
        all_data[map_info["name"]] = load_map(map_info["path"])
    summary_tsv = make_chromosome_summary(all_data)
    make_rug_map_svg(all_data, "cM", "Genetic position, cM",
                     FIG_DIR / "cicer_genetic_map_coverage.svg",
                     "Chickpea (Cicer arietinum) final genetic map: marker coverage by cM")
    make_rug_map_svg(all_data, "pos", "Physical position, Mb",
                     FIG_DIR / "cicer_physical_map_coverage.svg",
                     "Chickpea (Cicer arietinum) final genetic map: marker coverage by physical position")
    density_tsv = make_density_10cm(all_data)
    print("Done. Figures + tables written to results/figures and results/qc.")
    print(summary_tsv)
    print(density_tsv)


if __name__ == "__main__":
    main()

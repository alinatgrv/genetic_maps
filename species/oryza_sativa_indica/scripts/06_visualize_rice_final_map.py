#!/usr/bin/env python3

from pathlib import Path
import csv
import math
from collections import defaultdict


CONFIGS = [
    {
        "species_dir": Path("species/oryza_sativa_japonica"),
        "species_label": "Oryza sativa japonica",
        "map_label": "Final map on IRGSP-1.0, 1,619 bins",
        "prefix": "oryza_sativa_japonica",
        "map_path": Path("species/oryza_sativa_japonica/results/final/oryza_sativa_japonica_genetic_map.tsv"),
    },
    {
        "species_dir": Path("species/oryza_sativa_indica"),
        "species_label": "Oryza sativa indica",
        "map_label": "Strict monotonic projection on ASM465v1, 1,303 bins",
        "prefix": "oryza_sativa_indica",
        "map_path": Path("species/oryza_sativa_indica/results/final/oryza_sativa_indica_genetic_map.tsv"),
    },
]


def chr_sort_key(chrom):
    try:
        return int(chrom)
    except ValueError:
        return chrom


def svg_escape(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def write_svg(path, width, height, body):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as out:
        out.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n')
        out.write('<rect width="100%" height="100%" fill="white"/>\n')
        out.write(body)
        out.write("\n</svg>\n")


def load_map(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({
                "chr": str(row["chr"]),
                "pos": float(row["pos"]),
                "cM": float(row["cM"]),
            })
    return rows


def make_rug_map_svg(data, value_col, unit_label, out_path, title, panel_label):
    chroms = sorted({row["chr"] for row in data}, key=chr_sort_key)

    width = 1500
    height = 650
    top_margin = 90
    left_margin = 130
    right_margin = 80
    bottom_margin = 90

    plot_width = width - left_margin - right_margin
    plot_top = top_margin + 55
    plot_bottom = height - bottom_margin
    chr_step = plot_width / len(chroms)

    all_values = [row[value_col] for row in data]
    global_min = 0
    global_max = max(all_values)

    if value_col == "pos":
        global_max = math.ceil(global_max / 10_000_000) * 10_000_000
    else:
        global_max = math.ceil(global_max / 10) * 10

    if global_max <= global_min:
        global_max = global_min + 1

    def y_scale(value):
        return plot_top + (value - global_min) / (global_max - global_min) * (plot_bottom - plot_top)

    by_chr = defaultdict(list)
    for row in data:
        by_chr[row["chr"]].append(row)

    body = []
    body.append(f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">{svg_escape(title)}</text>\n')
    body.append(f'<text x="{left_margin}" y="{top_margin}" font-family="Arial" font-size="20" font-weight="bold">{svg_escape(panel_label)}</text>\n')

    body.append(f'<line x1="{left_margin-25}" y1="{plot_top}" x2="{left_margin-25}" y2="{plot_bottom}" stroke="black" stroke-width="1"/>\n')

    n_ticks = 5
    for t in range(n_ticks + 1):
        value = global_min + (global_max - global_min) * t / n_ticks
        y = y_scale(value)
        if value_col == "pos":
            label = f"{value/1_000_000:.0f}"
        else:
            label = f"{value:.0f}"

        body.append(f'<line x1="{left_margin-30}" y1="{y}" x2="{left_margin-25}" y2="{y}" stroke="black" stroke-width="1"/>\n')
        body.append(f'<text x="{left_margin-38}" y="{y+4}" text-anchor="end" font-family="Arial" font-size="12">{label}</text>\n')

    body.append(f'<text x="{left_margin-85}" y="{(plot_top+plot_bottom)/2}" transform="rotate(-90 {left_margin-85},{(plot_top+plot_bottom)/2})" text-anchor="middle" font-family="Arial" font-size="14">{svg_escape(unit_label)}</text>\n')

    for i, chrom in enumerate(chroms):
        x = left_margin + chr_step * (i + 0.5)
        vals = [r[value_col] for r in by_chr[chrom]]

        if vals:
            chr_min = min(vals)
            chr_max = max(vals)
            y1 = y_scale(chr_min)
            y2 = y_scale(chr_max)
        else:
            y1, y2 = plot_top, plot_bottom

        body.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" stroke="black" stroke-width="2"/>\n')

        for r in sorted(by_chr[chrom], key=lambda z: z[value_col]):
            y = y_scale(r[value_col])
            body.append(f'<line x1="{x-9}" y1="{y}" x2="{x+9}" y2="{y}" stroke="black" stroke-width="1"/>\n')

        body.append(f'<text x="{x}" y="{plot_bottom + 30}" text-anchor="middle" font-family="Arial" font-size="14">chr{svg_escape(chrom)}</text>\n')
        body.append(f'<text x="{x}" y="{plot_bottom + 50}" text-anchor="middle" font-family="Arial" font-size="12">{len(by_chr[chrom])}</text>\n')

    body.append(f'<text x="{width - right_margin}" y="{plot_bottom + 50}" text-anchor="end" font-family="Arial" font-size="12">numbers = markers per chromosome</text>\n')

    write_svg(out_path, width, height, "".join(body))


def make_density_10cm(data, prefix, qc_dir, fig_dir, title):
    bin_size = 10.0
    density_rows = []

    by_chr = defaultdict(list)
    for row in data:
        by_chr[row["chr"]].append(row["cM"])

    for chrom, values in sorted(by_chr.items(), key=lambda x: chr_sort_key(x[0])):
        if not values:
            continue

        cmin = math.floor(min(values) / bin_size) * bin_size
        cmax = math.ceil(max(values) / bin_size) * bin_size

        current = cmin
        while current < cmax:
            left = current
            right = current + bin_size
            count = sum(1 for v in values if left <= v < right or (v == cmax and right == cmax))
            density_rows.append({
                "chr": chrom,
                "bin_start_cM": left,
                "bin_end_cM": right,
                "n_markers": count,
            })
            current += bin_size

    out_tsv = qc_dir / f"{prefix}_marker_density_10cM_bins.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    with open(out_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["chr", "bin_start_cM", "bin_end_cM", "n_markers"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(density_rows)

    out_svg = fig_dir / f"{prefix}_marker_density_10cM.svg"
    make_density_svg(density_rows, out_svg, title)

    return out_tsv, out_svg


def make_density_svg(density_rows, out_path, title):
    chroms = sorted({r["chr"] for r in density_rows}, key=chr_sort_key)
    max_bin = max(r["bin_end_cM"] for r in density_rows)
    max_count = max(r["n_markers"] for r in density_rows) if density_rows else 1

    bin_size = 10
    bins = list(range(0, int(math.ceil(max_bin / bin_size) * bin_size), bin_size))

    cell_w = 92
    cell_h = 28
    left_margin = 150
    top_margin = 90
    title_h = 45

    width = left_margin + len(chroms) * cell_w + 180
    height = top_margin + title_h + len(bins) * cell_h + 90

    body = []
    body.append(f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">{svg_escape(title)}</text>\n')

    by_key = {(r["chr"], int(r["bin_start_cM"])): r["n_markers"] for r in density_rows}

    for i, chrom in enumerate(chroms):
        x = left_margin + i * cell_w
        body.append(f'<text x="{x + cell_w/2}" y="{top_margin}" text-anchor="middle" font-family="Arial" font-size="13">chr{svg_escape(chrom)}</text>\n')

    for j, bin_start in enumerate(bins):
        y = top_margin + title_h + j * cell_h
        body.append(f'<text x="{left_margin - 15}" y="{y + 18}" text-anchor="end" font-family="Arial" font-size="12">{bin_start}-{bin_start+bin_size}</text>\n')

        for i, chrom in enumerate(chroms):
            x = left_margin + i * cell_w
            count = by_key.get((chrom, bin_start), 0)
            shade = 255 - int(200 * (count / max_count)) if max_count else 255
            shade = max(35, min(255, shade))
            fill = f'rgb({shade},{shade},{shade})'

            body.append(f'<rect x="{x}" y="{y}" width="{cell_w-4}" height="{cell_h-3}" fill="{fill}" stroke="white"/>\n')
            if count > 0:
                text_color = "white" if shade < 120 else "black"
                body.append(f'<text x="{x + (cell_w-4)/2}" y="{y + 18}" text-anchor="middle" font-family="Arial" font-size="12" fill="{text_color}">{count}</text>\n')

    body.append(f'<text x="{left_margin - 80}" y="{top_margin + title_h + len(bins)*cell_h/2}" transform="rotate(-90 {left_margin-80},{top_margin + title_h + len(bins)*cell_h/2})" text-anchor="middle" font-family="Arial" font-size="14">cM interval</text>\n')

    legend_x = width - 130
    legend_y = 90
    body.append(f'<text x="{legend_x}" y="{legend_y}" font-family="Arial" font-size="13">n markers</text>\n')
    for k in range(6):
        count = round(max_count * k / 5)
        shade = 255 - int(200 * (count / max_count)) if max_count else 255
        shade = max(35, min(255, shade))
        body.append(f'<rect x="{legend_x}" y="{legend_y + 15 + k*24}" width="25" height="18" fill="rgb({shade},{shade},{shade})" stroke="black" stroke-width="0.5"/>\n')
        body.append(f'<text x="{legend_x + 35}" y="{legend_y + 29 + k*24}" font-family="Arial" font-size="12">{count}</text>\n')

    write_svg(out_path, width, height, "".join(body))


def main():
    for cfg in CONFIGS:
        species_dir = cfg["species_dir"]
        fig_dir = species_dir / "results/figures"
        qc_dir = species_dir / "results/qc"
        fig_dir.mkdir(parents=True, exist_ok=True)
        qc_dir.mkdir(parents=True, exist_ok=True)

        if not cfg["map_path"].exists():
            raise FileNotFoundError(cfg["map_path"])

        data = load_map(cfg["map_path"])

        genetic_svg = fig_dir / f"{cfg['prefix']}_genetic_map_coverage.svg"
        physical_svg = fig_dir / f"{cfg['prefix']}_physical_map_coverage.svg"

        make_rug_map_svg(
            data,
            value_col="cM",
            unit_label="Genetic position, cM",
            out_path=genetic_svg,
            title=f"{cfg['species_label']}: marker coverage by cM",
            panel_label=cfg["map_label"],
        )

        make_rug_map_svg(
            data,
            value_col="pos",
            unit_label="Physical position, Mb",
            out_path=physical_svg,
            title=f"{cfg['species_label']}: marker coverage by physical position",
            panel_label=cfg["map_label"],
        )

        density_tsv, density_svg = make_density_10cm(
            data,
            prefix=cfg["prefix"],
            qc_dir=qc_dir,
            fig_dir=fig_dir,
            title=f"{cfg['species_label']}: marker density in 10 cM bins",
        )

        print(cfg["species_label"])
        print("  Figures:")
        print(f"  {genetic_svg}")
        print(f"  {physical_svg}")
        print(f"  {density_svg}")
        print("  Density table:")
        print(f"  {density_tsv}")


if __name__ == "__main__":
    main()

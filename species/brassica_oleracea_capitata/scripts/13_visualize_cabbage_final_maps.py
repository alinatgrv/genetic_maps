from pathlib import Path
import csv
import math
from collections import defaultdict


SPECIES_DIR = Path("species/brassica_oleracea_capitata")

MAPS = [
    {
        "name": "strict_high_confidence_196",
        "label": "Strict high-confidence, 196 markers",
        "path": SPECIES_DIR / "results/final/cabbage_genetic_map.high_confidence_196.tsv",
    },
    {
        "name": "relaxed_mm1_collinear_222",
        "label": "Relaxed mm1 collinear, 222 markers",
        "path": SPECIES_DIR / "results/final/cabbage_genetic_map.relaxed_mm1_collinear_222.tsv",
    },
]

FIG_DIR = SPECIES_DIR / "results/figures"
QC_DIR = SPECIES_DIR / "results/qc"

FIG_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)


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
    with open(path, "w", encoding="utf-8") as out:
        out.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">\n')
        out.write('<rect width="100%" height="100%" fill="white"/>\n')
        out.write(body)
        out.write("\n</svg>\n")


def make_rug_map_svg(all_data, value_col, unit_label, out_path, title):
    chroms = sorted(
        {row["chr"] for data in all_data.values() for row in data},
        key=chr_sort_key
    )

    width = 1500
    panel_height = 520
    top_margin = 80
    left_margin = 130
    right_margin = 80
    bottom_margin = 70
    panel_gap = 60

    height = top_margin + len(MAPS) * panel_height + (len(MAPS) - 1) * panel_gap + bottom_margin

    plot_width = width - left_margin - right_margin
    chr_step = plot_width / len(chroms)

    all_values = [row[value_col] for data in all_data.values() for row in data]
    global_min = 0
    global_max = max(all_values)
    if value_col == "pos":
        global_max = math.ceil(global_max / 10_000_000) * 10_000_000
    else:
        global_max = math.ceil(global_max / 10) * 10

    def y_scale(value, panel_top):
        plot_top = panel_top + 60
        plot_bottom = panel_top + panel_height - 70
        return plot_top + (value - global_min) / (global_max - global_min) * (plot_bottom - plot_top)

    body = []
    body.append(f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">{svg_escape(title)}</text>\n')

    for panel_i, map_info in enumerate(MAPS):
        map_name = map_info["name"]
        map_label = map_info["label"]
        data = all_data[map_name]

        panel_top = top_margin + panel_i * (panel_height + panel_gap)
        plot_top = panel_top + 60
        plot_bottom = panel_top + panel_height - 70

        body.append(f'<text x="{left_margin}" y="{panel_top + 25}" font-family="Arial" font-size="20" font-weight="bold">{svg_escape(map_label)}</text>\n')

        # y-axis
        body.append(f'<line x1="{left_margin-25}" y1="{plot_top}" x2="{left_margin-25}" y2="{plot_bottom}" stroke="black" stroke-width="1"/>\n')

        n_ticks = 5
        for t in range(n_ticks + 1):
            value = global_min + (global_max - global_min) * t / n_ticks
            y = y_scale(value, panel_top)
            if value_col == "pos":
                label = f"{value/1_000_000:.0f}"
            else:
                label = f"{value:.0f}"

            body.append(f'<line x1="{left_margin-30}" y1="{y}" x2="{left_margin-25}" y2="{y}" stroke="black" stroke-width="1"/>\n')
            body.append(f'<text x="{left_margin-38}" y="{y+4}" text-anchor="end" font-family="Arial" font-size="12">{label}</text>\n')

        body.append(f'<text x="{left_margin-80}" y="{(plot_top+plot_bottom)/2}" transform="rotate(-90 {left_margin-80},{(plot_top+plot_bottom)/2})" text-anchor="middle" font-family="Arial" font-size="14">{svg_escape(unit_label)}</text>\n')

        by_chr = defaultdict(list)
        for row in data:
            by_chr[row["chr"]].append(row)

        for i, chrom in enumerate(chroms):
            x = left_margin + chr_step * (i + 0.5)

            vals = [r[value_col] for r in by_chr[chrom]]
            if vals:
                chr_min = min(vals)
                chr_max = max(vals)
                y1 = y_scale(chr_min, panel_top)
                y2 = y_scale(chr_max, panel_top)
            else:
                y1, y2 = plot_top, plot_bottom

            body.append(f'<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" stroke="black" stroke-width="2"/>\n')

            # marker ticks
            for r in sorted(by_chr[chrom], key=lambda z: z[value_col]):
                y = y_scale(r[value_col], panel_top)
                body.append(f'<line x1="{x-9}" y1="{y}" x2="{x+9}" y2="{y}" stroke="black" stroke-width="1"/>\n')

            body.append(f'<text x="{x}" y="{plot_bottom + 28}" text-anchor="middle" font-family="Arial" font-size="14">chr{svg_escape(chrom)}</text>\n')
            body.append(f'<text x="{x}" y="{plot_bottom + 47}" text-anchor="middle" font-family="Arial" font-size="12">{len(by_chr[chrom])}</text>\n')

        body.append(f'<text x="{width - right_margin}" y="{plot_bottom + 47}" text-anchor="end" font-family="Arial" font-size="12">numbers = markers per chromosome</text>\n')

    write_svg(out_path, width, height, "".join(body))


def make_density_10cm(all_data):
    bin_size = 10.0
    density_rows = []

    for map_info in MAPS:
        map_name = map_info["name"]
        data = all_data[map_name]

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
                    "map_name": map_name,
                    "chr": chrom,
                    "bin_start_cM": left,
                    "bin_end_cM": right,
                    "n_markers": count,
                })
                current += bin_size

    out_tsv = QC_DIR / "cabbage_marker_density_10cM_bins.tsv"
    with open(out_tsv, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["map_name", "chr", "bin_start_cM", "bin_end_cM", "n_markers"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(density_rows)

    make_density_svg(density_rows, FIG_DIR / "cabbage_marker_density_10cM.svg")

    return out_tsv


def make_density_svg(density_rows, out_path):
    chroms = sorted({r["chr"] for r in density_rows}, key=chr_sort_key)
    max_bin = max(r["bin_end_cM"] for r in density_rows)
    max_count = max(r["n_markers"] for r in density_rows) if density_rows else 1

    bin_size = 10
    bins = list(range(0, int(math.ceil(max_bin / bin_size) * bin_size), bin_size))

    cell_w = 95
    cell_h = 28
    left_margin = 150
    top_margin = 80
    panel_gap = 80
    title_h = 40

    width = left_margin + len(chroms) * cell_w + 180
    panel_height = title_h + len(bins) * cell_h + 50
    height = top_margin + len(MAPS) * panel_height + (len(MAPS)-1) * panel_gap + 80

    body = []
    body.append(f'<text x="{width/2}" y="38" text-anchor="middle" font-family="Arial" font-size="26" font-weight="bold">Marker density in 10 cM bins</text>\n')

    by_key = {(r["map_name"], r["chr"], int(r["bin_start_cM"])): r["n_markers"] for r in density_rows}

    for panel_i, map_info in enumerate(MAPS):
        map_name = map_info["name"]
        map_label = map_info["label"]

        panel_top = top_margin + panel_i * (panel_height + panel_gap)

        body.append(f'<text x="{left_margin}" y="{panel_top}" font-family="Arial" font-size="20" font-weight="bold">{svg_escape(map_label)}</text>\n')

        for i, chrom in enumerate(chroms):
            x = left_margin + i * cell_w
            body.append(f'<text x="{x + cell_w/2}" y="{panel_top + 30}" text-anchor="middle" font-family="Arial" font-size="13">chr{svg_escape(chrom)}</text>\n')

        for j, bin_start in enumerate(bins):
            y = panel_top + title_h + j * cell_h
            body.append(f'<text x="{left_margin - 15}" y="{y + 18}" text-anchor="end" font-family="Arial" font-size="12">{bin_start}-{bin_start+bin_size}</text>\n')

            for i, chrom in enumerate(chroms):
                x = left_margin + i * cell_w
                count = by_key.get((map_name, chrom, bin_start), 0)
                shade = 255 - int(200 * (count / max_count)) if max_count else 255
                shade = max(35, min(255, shade))
                fill = f'rgb({shade},{shade},{shade})'

                body.append(f'<rect x="{x}" y="{y}" width="{cell_w-4}" height="{cell_h-3}" fill="{fill}" stroke="white"/>\n')
                if count > 0:
                    text_color = "white" if shade < 120 else "black"
                    body.append(f'<text x="{x + (cell_w-4)/2}" y="{y + 18}" text-anchor="middle" font-family="Arial" font-size="12" fill="{text_color}">{count}</text>\n')

        body.append(f'<text x="{left_margin - 80}" y="{panel_top + title_h + len(bins)*cell_h/2}" transform="rotate(-90 {left_margin-80},{panel_top + title_h + len(bins)*cell_h/2})" text-anchor="middle" font-family="Arial" font-size="14">cM interval</text>\n')

    # legend
    legend_x = width - 130
    legend_y = 80
    body.append(f'<text x="{legend_x}" y="{legend_y}" font-family="Arial" font-size="13">n markers</text>\n')
    for k in range(6):
        count = round(max_count * k / 5)
        shade = 255 - int(200 * (count / max_count)) if max_count else 255
        shade = max(35, min(255, shade))
        body.append(f'<rect x="{legend_x}" y="{legend_y + 15 + k*24}" width="25" height="18" fill="rgb({shade},{shade},{shade})" stroke="black" stroke-width="0.5"/>\n')
        body.append(f'<text x="{legend_x + 35}" y="{legend_y + 29 + k*24}" font-family="Arial" font-size="12">{count}</text>\n')

    write_svg(out_path, width, height, "".join(body))


def main():
    all_data = {}
    for map_info in MAPS:
        if not map_info["path"].exists():
            raise FileNotFoundError(map_info["path"])
        all_data[map_info["name"]] = load_map(map_info["path"])

    make_rug_map_svg(
        all_data,
        value_col="cM",
        unit_label="Genetic position, cM",
        out_path=FIG_DIR / "cabbage_genetic_map_coverage.svg",
        title="Cabbage final genetic maps: marker coverage by cM"
    )

    make_rug_map_svg(
        all_data,
        value_col="pos",
        unit_label="Physical position, Mb",
        out_path=FIG_DIR / "cabbage_physical_map_coverage.svg",
        title="Cabbage final genetic maps: marker coverage by physical position"
    )

    density_tsv = make_density_10cm(all_data)

    print("Done.")
    print("Figures:")
    print(FIG_DIR / "cabbage_genetic_map_coverage.svg")
    print(FIG_DIR / "cabbage_physical_map_coverage.svg")
    print(FIG_DIR / "cabbage_marker_density_10cM.svg")
    print("Density table:")
    print(density_tsv)


if __name__ == "__main__":
    main()

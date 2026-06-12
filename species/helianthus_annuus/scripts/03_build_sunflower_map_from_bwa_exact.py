from pathlib import Path
from collections import Counter, defaultdict
import csv

SPECIES_DIR = Path("species/helianthus_annuus")

sam_file = SPECIES_DIR / "results/intermediate/sunflower_features_25bp.bwa_exact.sam"
metadata_file = SPECIES_DIR / "data/metadata/sunflower_features_metadata.tsv"

out_all = SPECIES_DIR / "results/intermediate/sunflower_genetic_map.bwa_exact_all_unique_chrom_hits.tsv"
out_lg_matched = SPECIES_DIR / "results/intermediate/sunflower_genetic_map.bwa_exact_chr_lg_matched_with_markers.tsv"
out_final = SPECIES_DIR / "results/final/sunflower_genetic_map.bwa_exact_unique.tsv"

out_qc = SPECIES_DIR / "results/qc/sunflower_bwa_exact_map_summary.txt"
out_unmapped = SPECIES_DIR / "results/qc/sunflower_bwa_exact_unmapped_features.txt"
out_multimap = SPECIES_DIR / "results/qc/sunflower_bwa_exact_multimap_or_nonunique_features.tsv"
out_lg_conflicts = SPECIES_DIR / "results/qc/sunflower_bwa_exact_chr_lg_conflicts.tsv"
out_chr_pos_conflicts = SPECIES_DIR / "results/qc/sunflower_bwa_exact_chr_pos_cm_conflicts.tsv"

ref_to_chr = {
    "NC_035433.2": "1",
    "NC_035434.2": "2",
    "NC_035435.2": "3",
    "NC_035436.2": "4",
    "NC_035437.2": "5",
    "NC_035438.2": "6",
    "NC_035439.2": "7",
    "NC_035440.2": "8",
    "NC_035441.2": "9",
    "NC_035442.2": "10",
    "NC_035443.2": "11",
    "NC_035444.2": "12",
    "NC_035445.2": "13",
    "NC_035446.2": "14",
    "NC_035447.2": "15",
    "NC_035448.2": "16",
    "NC_035449.2": "17",
}

def parse_sam_tag(fields, tag):
    prefix = tag + ":"
    for f in fields:
        if f.startswith(prefix):
            return f.split(":")[-1]
    return None

def write_tsv(path, rows, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

metadata = {}
with open(metadata_file, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        metadata[row["feature_id"]] = row

total_sam_rows = 0
unmapped = []
multimap = []
all_unique_chrom_hits = []
lg_conflicts = []

x0_counter = Counter()
rname_counter = Counter()
flag_counter = Counter()

with open(sam_file, encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith("@"):
            continue

        total_sam_rows += 1
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            raise ValueError(f"Malformed SAM line with {len(fields)} fields: {line[:200]}")

        qname = fields[0]
        flag = int(fields[1])
        rname = fields[2]
        pos = int(fields[3]) if fields[3].isdigit() else 0
        mapq = fields[4]
        cigar = fields[5]
        seq = fields[9]

        flag_counter[str(flag)] += 1

        if (flag & 4) or rname == "*":
            unmapped.append(qname)
            continue

        nm = parse_sam_tag(fields[11:], "NM")
        x0 = parse_sam_tag(fields[11:], "X0")
        x1 = parse_sam_tag(fields[11:], "X1")

        x0_counter[str(x0)] += 1
        rname_counter[rname] += 1

        # Strict exact unique 25-bp hit:
        # 25M, NM=0, X0=1, chromosome-level reference only.
        if not (cigar == "25M" and nm == "0" and x0 == "1"):
            multimap.append({
                "feature_id": qname,
                "ref_seq": rname,
                "alignment_start": pos,
                "cigar": cigar,
                "MAPQ": mapq,
                "NM": nm,
                "X0": x0,
                "X1": x1,
            })
            continue

        if rname not in ref_to_chr:
            multimap.append({
                "feature_id": qname,
                "ref_seq": rname,
                "alignment_start": pos,
                "cigar": cigar,
                "MAPQ": mapq,
                "NM": nm,
                "X0": x0,
                "X1": x1,
            })
            continue

        if qname not in metadata:
            raise ValueError(f"{qname} not found in metadata")

        meta = metadata[qname]
        chrom = ref_to_chr[rname]
        marker_pos = pos + 12  # midpoint of 25M alignment, 1-based SAM coordinate

        row = {
            "chr": chrom,
            "pos": str(marker_pos),
            "cM": meta["cM"],
            "feature_id": qname,
            "linkage_group": meta["linkage_group"],
            "unigene": meta["unigene"],
            "plants_scored": meta["plants_scored"],
            "mismatches_to_template": meta["mismatches_to_template"],
            "sequence": meta["sequence"],
            "ref_seq": rname,
            "alignment_start": str(pos),
            "cigar": cigar,
            "MAPQ": mapq,
            "NM": nm,
            "X0": x0,
            "X1": x1,
        }

        all_unique_chrom_hits.append(row)

        # For final map we keep only markers where published linkage group
        # agrees with physical chromosome number.
        if str(meta["linkage_group"]) != chrom:
            lg_conflicts.append(row)

lg_conflict_ids = {r["feature_id"] for r in lg_conflicts}
lg_matched = [r for r in all_unique_chrom_hits if r["feature_id"] not in lg_conflict_ids]

# Exact deduplication by chr-pos-cM
seen_exact = set()
exact_dedup = []
for r in lg_matched:
    key = (r["chr"], r["pos"], r["cM"])
    if key not in seen_exact:
        exact_dedup.append(r)
        seen_exact.add(key)

# Exclude chr-pos positions with conflicting cM values
cm_by_chr_pos = defaultdict(set)
for r in exact_dedup:
    cm_by_chr_pos[(r["chr"], r["pos"])].add(r["cM"])

conflict_chr_pos = {key for key, cms in cm_by_chr_pos.items() if len(cms) > 1}
chr_pos_conflict_rows = [
    r for r in exact_dedup
    if (r["chr"], r["pos"]) in conflict_chr_pos
]

final_rows_with_markers = [
    r for r in exact_dedup
    if (r["chr"], r["pos"]) not in conflict_chr_pos
]

def sort_key(row):
    return (int(row["chr"]), int(row["pos"]), float(row["cM"]), row.get("feature_id", ""))

all_unique_chrom_hits.sort(key=sort_key)
lg_matched.sort(key=sort_key)
lg_conflicts.sort(key=sort_key)
chr_pos_conflict_rows.sort(key=sort_key)
final_rows_with_markers.sort(key=sort_key)

final_simple = [
    {"chr": r["chr"], "pos": r["pos"], "cM": r["cM"]}
    for r in final_rows_with_markers
]

header_all = [
    "chr", "pos", "cM", "feature_id", "linkage_group", "unigene",
    "plants_scored", "mismatches_to_template", "sequence",
    "ref_seq", "alignment_start", "cigar", "MAPQ", "NM", "X0", "X1"
]

write_tsv(out_all, all_unique_chrom_hits, header_all)
write_tsv(out_lg_matched, lg_matched, header_all)
write_tsv(out_final, final_simple, ["chr", "pos", "cM"])
write_tsv(out_multimap, multimap, ["feature_id", "ref_seq", "alignment_start", "cigar", "MAPQ", "NM", "X0", "X1"])
write_tsv(out_lg_conflicts, lg_conflicts, header_all)
write_tsv(out_chr_pos_conflicts, chr_pos_conflict_rows, header_all)

with open(out_unmapped, "w") as out:
    for q in unmapped:
        out.write(q + "\n")

chr_counts = Counter(r["chr"] for r in final_rows_with_markers)

with open(out_qc, "w") as out:
    out.write("Sunflower BWA exact high-confidence map summary\n")
    out.write("==============================================\n\n")
    out.write(f"Input features in metadata: {len(metadata)}\n")
    out.write(f"SAM alignment rows: {total_sam_rows}\n")
    out.write(f"Unmapped features: {len(unmapped)}\n")
    out.write(f"Exact unique chromosome-level hits before LG filtering: {len(all_unique_chrom_hits)}\n")
    out.write(f"LG/chr conflict rows excluded: {len(lg_conflicts)}\n")
    out.write(f"Rows after LG/chr filtering: {len(lg_matched)}\n")
    out.write(f"Rows after exact chr-pos-cM deduplication: {len(exact_dedup)}\n")
    out.write(f"chr-pos coordinates with conflicting cM values: {len(conflict_chr_pos)}\n")
    out.write(f"Rows excluded due to chr-pos cM conflicts: {len(chr_pos_conflict_rows)}\n")
    out.write(f"Final unique chr-pos-cM rows: {len(final_simple)}\n\n")

    out.write("BWA X0 tag distribution:\n")
    for k in sorted(x0_counter, key=lambda x: (-1 if x == "None" else int(x))):
        out.write(f"{k}\t{x0_counter[k]}\n")

    out.write("\nFinal markers by chromosome:\n")
    for chrom in sorted(chr_counts, key=lambda x: int(x)):
        out.write(f"chr{chrom}\t{chr_counts[chrom]}\n")

print("Done.")
print(f"All exact unique chromosome hits: {out_all}")
print(f"LG-matched map with markers: {out_lg_matched}")
print(f"Final map: {out_final}")
print(f"QC summary: {out_qc}")
print(f"Final unique chr-pos-cM rows: {len(final_simple)}")

import pandas as pd
from pathlib import Path

blast_file = Path("results/cabbage/cabbage_primers_vs_ref.full_length_mm1.tsv")
metadata_file = Path("data/metadata/cabbage/cabbage_markers_metadata.tsv")

amplicons_out = Path("results/cabbage/cabbage_amplicons.relaxed_primer_mm1.tsv")
map_out = Path("results/cabbage/cabbage_genetic_map.relaxed_primer_mm1.tsv")
summary_out = Path("results/cabbage/cabbage_amplicons.relaxed_primer_mm1.summary.txt")

MAX_AMPLICON_LEN = 1000

cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen"
]

blast = pd.read_csv(blast_file, sep="\t", names=cols)

# Safety filter:
# full primer length, no gaps, maximum 1 mismatch per primer
blast = blast[
    (blast["length"] == blast["qlen"]) &
    (blast["gapopen"] == 0) &
    (blast["mismatch"] <= 1)
].copy()

parsed = blast["qseqid"].str.extract(r"^(?P<marker_id>.+)__(?P<primer_side>[FR])$")
blast["marker_id"] = parsed["marker_id"]
blast["primer_side"] = parsed["primer_side"]
blast = blast.dropna(subset=["marker_id", "primer_side"])

blast["hit_start"] = blast[["sstart", "send"]].min(axis=1)
blast["hit_end"] = blast[["sstart", "send"]].max(axis=1)
blast["strand"] = blast.apply(lambda row: "+" if row["sstart"] <= row["send"] else "-", axis=1)

amplicons = []

for marker_id, group in blast.groupby("marker_id"):
    f_hits = group[group["primer_side"] == "F"].copy()
    r_hits = group[group["primer_side"] == "R"].copy()

    if f_hits.empty or r_hits.empty:
        continue

    for _, f in f_hits.iterrows():
        for _, r in r_hits.iterrows():

            if f["sseqid"] != r["sseqid"]:
                continue

            orientation_ok = False

            # Case 1: F on plus strand, R on minus strand, F is left of R
            if f["strand"] == "+" and r["strand"] == "-":
                if f["hit_start"] <= r["hit_start"]:
                    orientation_ok = True

            # Case 2: F on minus strand, R on plus strand, R is left of F
            elif f["strand"] == "-" and r["strand"] == "+":
                if r["hit_start"] <= f["hit_start"]:
                    orientation_ok = True

            if not orientation_ok:
                continue

            amp_start = int(min(f["hit_start"], r["hit_start"]))
            amp_end = int(max(f["hit_end"], r["hit_end"]))
            amp_len = amp_end - amp_start + 1

            if amp_len > MAX_AMPLICON_LEN:
                continue

            marker_pos = int(round((amp_start + amp_end) / 2))
            total_mismatch = int(f["mismatch"] + r["mismatch"])

            amplicons.append({
                "marker_id": marker_id,
                "ref_seq": f["sseqid"],
                "amp_start": amp_start,
                "amp_end": amp_end,
                "amp_len": amp_len,
                "marker_pos": marker_pos,
                "forward_mismatch": int(f["mismatch"]),
                "reverse_mismatch": int(r["mismatch"]),
                "total_mismatch": total_mismatch,
                "forward_hit_start": int(f["hit_start"]),
                "forward_hit_end": int(f["hit_end"]),
                "forward_strand": f["strand"],
                "reverse_hit_start": int(r["hit_start"]),
                "reverse_hit_end": int(r["hit_end"]),
                "reverse_strand": r["strand"],
            })

amplicons = pd.DataFrame(amplicons)

if amplicons.empty:
    raise RuntimeError("No relaxed primer-mm1 candidate amplicons found.")

# Keep only markers with exactly one candidate amplicon
amplicons["candidate_amplicons_per_marker"] = (
    amplicons.groupby("marker_id")["marker_id"].transform("count")
)

unique_amplicons = amplicons[
    amplicons["candidate_amplicons_per_marker"] == 1
].copy()

metadata = pd.read_csv(metadata_file, sep="\t")

genetic_map = unique_amplicons.merge(metadata, on="marker_id", how="left")

genetic_map = genetic_map[[
    "ref_seq",
    "marker_pos",
    "cM",
    "marker_id",
    "linkage_group",
    "marker_type",
    "amp_start",
    "amp_end",
    "amp_len",
    "forward_mismatch",
    "reverse_mismatch",
    "total_mismatch"
]].copy()

genetic_map = genetic_map.rename(columns={
    "ref_seq": "chr",
    "marker_pos": "pos"
})

genetic_map = genetic_map.sort_values(["chr", "pos", "cM", "marker_id"])

amplicons.to_csv(amplicons_out, sep="\t", index=False)
genetic_map.to_csv(map_out, sep="\t", index=False)

with open(summary_out, "w") as out:
    out.write("Cabbage relaxed primer-mm1 amplicon summary\n")
    out.write("==========================================\n\n")
    out.write(f"Full-length <=1 mismatch primer hits: {len(blast)}\n")
    out.write(f"Unique primers with full-length <=1 mismatch hits: {blast['qseqid'].nunique()}\n")
    out.write(f"Unique markers with at least one relaxed primer hit: {blast['marker_id'].nunique()}\n")
    out.write(f"Candidate amplicons found: {len(amplicons)}\n")
    out.write(f"Markers with at least one candidate amplicon: {amplicons['marker_id'].nunique()}\n")
    out.write(f"Markers with exactly one candidate amplicon: {unique_amplicons['marker_id'].nunique()}\n")
    out.write(f"Maximum allowed amplicon length: {MAX_AMPLICON_LEN} bp\n\n")

    out.write("Unique amplicons by total mismatch:\n")
    out.write(unique_amplicons["total_mismatch"].value_counts().sort_index().to_string())
    out.write("\n\n")

    out.write("Unique amplicons by marker type:\n")
    out.write(genetic_map["marker_type"].value_counts().to_string())
    out.write("\n\n")

    out.write("Amplicon length summary:\n")
    out.write(unique_amplicons["amp_len"].describe().to_string())
    out.write("\n")

print(f"Relaxed candidate amplicons: {amplicons_out}")
print(f"Relaxed primer-mm1 map: {map_out}")
print(f"Summary: {summary_out}")

print("\nSummary:")
print(f"Full-length <=1 mismatch primer hits: {len(blast)}")
print(f"Unique primers with full-length <=1 mismatch hits: {blast['qseqid'].nunique()}")
print(f"Candidate amplicons found: {len(amplicons)}")
print(f"Markers with candidate amplicons: {amplicons['marker_id'].nunique()}")
print(f"Markers with exactly one candidate amplicon: {unique_amplicons['marker_id'].nunique()}")

print("\nUnique amplicons by total mismatch:")
print(unique_amplicons["total_mismatch"].value_counts().sort_index().to_string())

print("\nUnique amplicons by marker type:")
print(genetic_map["marker_type"].value_counts().to_string())

print("\nAmplicon length summary:")
print(unique_amplicons["amp_len"].describe().to_string())

import pandas as pd
from pathlib import Path

# =========================
# Input / output
# =========================

blast_file = Path("results/cabbage/cabbage_primers_vs_ref.perfect.tsv")
metadata_file = Path("data/metadata/cabbage/cabbage_markers_metadata.tsv")

amplicons_out = Path("results/cabbage/cabbage_amplicons.strict.tsv")
map_out = Path("results/cabbage/cabbage_genetic_map.strict.tsv")
summary_out = Path("results/cabbage/cabbage_amplicons.strict.summary.txt")

# Maximum expected PCR product length.
# This is a conservative first-pass threshold.
MAX_AMPLICON_LEN = 1000

# =========================
# Load BLAST perfect hits
# =========================

cols = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen"
]

blast = pd.read_csv(blast_file, sep="\t", names=cols)

# Parse primer id:
# BoSF1331__F -> marker_id = BoSF1331, primer_side = F
parsed = blast["qseqid"].str.extract(r"^(?P<marker_id>.+)__(?P<primer_side>[FR])$")

blast["marker_id"] = parsed["marker_id"]
blast["primer_side"] = parsed["primer_side"]

# Coordinates independent of alignment strand
blast["hit_start"] = blast[["sstart", "send"]].min(axis=1)
blast["hit_end"] = blast[["sstart", "send"]].max(axis=1)

# BLAST strand relative to reference:
# plus if sstart < send, minus if sstart > send
blast["strand"] = blast.apply(
    lambda row: "+" if row["sstart"] <= row["send"] else "-",
    axis=1
)

# Keep only correctly parsed primer IDs
blast = blast.dropna(subset=["marker_id", "primer_side"])

# =========================
# Pair F and R hits into candidate amplicons
# =========================

amplicons = []

for marker_id, group in blast.groupby("marker_id"):
    f_hits = group[group["primer_side"] == "F"].copy()
    r_hits = group[group["primer_side"] == "R"].copy()

    if f_hits.empty or r_hits.empty:
        continue

    for _, f in f_hits.iterrows():
        for _, r in r_hits.iterrows():

            # Both primers must be on the same reference sequence
            if f["sseqid"] != r["sseqid"]:
                continue

            # Expected PCR orientation:
            # F primer on + strand, R primer on - strand
            # or F primer on - strand, R primer on + strand.
            # In both cases primers should face each other.
            orientation_ok = False

            if f["strand"] == "+" and r["strand"] == "-":
                # F is left, R is right
                if f["hit_start"] <= r["hit_start"]:
                    orientation_ok = True

            elif f["strand"] == "-" and r["strand"] == "+":
                # R is left, F is right
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

            amplicons.append({
                "marker_id": marker_id,
                "ref_seq": f["sseqid"],
                "amp_start": amp_start,
                "amp_end": amp_end,
                "amp_len": amp_len,
                "marker_pos": marker_pos,
                "forward_hit_start": int(f["hit_start"]),
                "forward_hit_end": int(f["hit_end"]),
                "forward_strand": f["strand"],
                "reverse_hit_start": int(r["hit_start"]),
                "reverse_hit_end": int(r["hit_end"]),
                "reverse_strand": r["strand"],
            })

amplicons = pd.DataFrame(amplicons)

if amplicons.empty:
    raise RuntimeError("No candidate amplicons found. Try increasing MAX_AMPLICON_LEN or relaxing filters.")

# Number of candidate amplicons per marker
amplicons["candidate_amplicons_per_marker"] = amplicons.groupby("marker_id")["marker_id"].transform("count")

# Strict: exactly one candidate amplicon per marker
strict = amplicons[amplicons["candidate_amplicons_per_marker"] == 1].copy()

# =========================
# Add marker metadata
# =========================

metadata = pd.read_csv(metadata_file, sep="\t")

strict = strict.merge(metadata, on="marker_id", how="left")

# Build final strict genetic map table
# Keep ref_seq first; later we will rename ref_seq to chromosome names if needed.
genetic_map = strict[[
    "ref_seq",
    "marker_pos",
    "cM",
    "marker_id",
    "linkage_group",
    "marker_type",
    "amp_start",
    "amp_end",
    "amp_len"
]].copy()

genetic_map = genetic_map.rename(columns={
    "ref_seq": "chr",
    "marker_pos": "pos"
})

genetic_map = genetic_map.sort_values(["chr", "pos", "cM"])

# =========================
# Save outputs
# =========================

amplicons_out.parent.mkdir(parents=True, exist_ok=True)

amplicons.to_csv(amplicons_out, sep="\t", index=False)
genetic_map.to_csv(map_out, sep="\t", index=False)

with open(summary_out, "w") as out:
    out.write("Cabbage strict amplicon mapping summary\n")
    out.write("======================================\n\n")
    out.write(f"BLAST perfect hits: {len(blast)}\n")
    out.write(f"Markers in metadata: {metadata['marker_id'].nunique()}\n")
    out.write(f"Markers with at least one perfect primer hit: {blast['marker_id'].nunique()}\n")
    out.write(f"Candidate amplicons found: {len(amplicons)}\n")
    out.write(f"Markers with at least one candidate amplicon: {amplicons['marker_id'].nunique()}\n")
    out.write(f"Strict markers with exactly one candidate amplicon: {strict['marker_id'].nunique()}\n")
    out.write(f"Maximum allowed amplicon length: {MAX_AMPLICON_LEN} bp\n\n")

    out.write("Amplicon length summary for strict markers:\n")
    out.write(strict["amp_len"].describe().to_string())
    out.write("\n\n")

    out.write("Strict markers by linkage group:\n")
    out.write(strict["linkage_group"].value_counts().sort_index().to_string())
    out.write("\n\n")

    out.write("Strict markers by marker type:\n")
    out.write(strict["marker_type"].value_counts().to_string())
    out.write("\n")

print(f"Candidate amplicons written to: {amplicons_out}")
print(f"Strict genetic map written to: {map_out}")
print(f"Summary written to: {summary_out}")

print("\nSummary:")
print(f"Candidate amplicons found: {len(amplicons)}")
print(f"Markers with candidate amplicons: {amplicons['marker_id'].nunique()}")
print(f"Strict markers with exactly one candidate amplicon: {strict['marker_id'].nunique()}")

print("\nStrict markers by linkage group:")
print(strict["linkage_group"].value_counts().sort_index().to_string())

print("\nAmplicon length summary:")
print(strict["amp_len"].describe().to_string())

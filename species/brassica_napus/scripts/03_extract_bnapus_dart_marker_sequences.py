#!/usr/bin/env python3

# Извлечение последовательностей DArT-клонов для маркеров консенсус-карты.
#
# Берём список join-id DArT-маркеров карты (из шага 02) и вытаскиваем
# соответствующие последовательности из общего файла DArT_Brassica.fasta
# (Diversity Arrays Technology). Заголовки FASTA имеют вид ">brPb-657581".
#
# На выходе — FASTA только нужных маркеров (имена приводим к join-id) и
# сводка по покрытию (сколько маркеров карты получили последовательность).
#
# Запускать из корня репозитория.

from pathlib import Path
from collections import Counter
import re

SPECIES_DIR = Path("species/brassica_napus")

dart_ids_file = SPECIES_DIR / "data/metadata/bnapus_consensus_dart_join_ids.txt"
metadata_file = SPECIES_DIR / "data/metadata/bnapus_consensus_map_metadata.tsv"
fasta_in = SPECIES_DIR / "data/markers/DArT_Brassica.fasta"
fasta_out = SPECIES_DIR / "data/markers/bnapus_dart_markers.fasta"
missing_out = SPECIES_DIR / "results/qc/bnapus_dart_markers_without_sequence.txt"
summary_out = SPECIES_DIR / "results/qc/bnapus_marker_sequence_extraction_summary.txt"

fasta_out.parent.mkdir(parents=True, exist_ok=True)
summary_out.parent.mkdir(parents=True, exist_ok=True)


def norm(name):
    return re.sub(r"^x", "", str(name).strip().lower())


# Нужные join-id.
wanted = set()
with open(dart_ids_file, encoding="utf-8") as handle:
    for line in handle:
        jid = line.strip()
        if jid:
            wanted.add(jid)


def iter_fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)


found = {}
seq_lengths = []
for name, seq in iter_fasta(fasta_in):
    jid = norm(name)
    if jid in wanted and jid not in found:
        found[jid] = seq
        seq_lengths.append(len(seq))

# Запись результата (отсортировано для воспроизводимости).
with open(fasta_out, "w", encoding="utf-8") as out:
    for jid in sorted(found):
        out.write(f">{jid}\n{found[jid]}\n")

missing = sorted(wanted - set(found))
with open(missing_out, "w", encoding="utf-8") as out:
    for jid in missing:
        out.write(jid + "\n")

# Покрытие по хромосомам (по метаданным карты).
chr_by_id = {}
with open(metadata_file, encoding="utf-8") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(header)}
    for line in handle:
        f = line.rstrip("\n").split("\t")
        if f[idx["is_dart"]] == "1":
            chr_by_id[f[idx["marker_join_id"]]] = f[idx["chr"]]

chr_order = [f"A{i}" for i in range(1, 11)] + [f"C{i}" for i in range(1, 10)]
found_by_chr = Counter(chr_by_id.get(jid, "?") for jid in found)
wanted_by_chr = Counter(chr_by_id.get(jid, "?") for jid in wanted)

if seq_lengths:
    seq_lengths_sorted = sorted(seq_lengths)
    median_len = seq_lengths_sorted[len(seq_lengths_sorted) // 2]
else:
    median_len = 0

with open(summary_out, "w", encoding="utf-8") as out:
    out.write("Brassica napus DArT marker sequence extraction\n")
    out.write("=" * 50 + "\n\n")
    out.write(f"Source FASTA: {fasta_in}\n")
    out.write(f"Requested DArT markers (in consensus map): {len(wanted)}\n")
    out.write(f"Markers with sequence found: {len(found)}\n")
    out.write(f"Markers without sequence: {len(missing)}\n")
    if seq_lengths:
        out.write(f"Sequence length: min {min(seq_lengths)} / median {median_len} / max {max(seq_lengths)} bp\n")
    out.write("\nCoverage by chromosome (found / requested):\n")
    for ch in chr_order:
        out.write(f"  {ch}\t{found_by_chr.get(ch, 0)}\t{wanted_by_chr.get(ch, 0)}\n")

print("Done.")
print(f"Marker FASTA: {fasta_out}")
print(f"Missing list: {missing_out}")
print(f"Summary:      {summary_out}")
print(f"Found sequences: {len(found)} / {len(wanted)} requested")

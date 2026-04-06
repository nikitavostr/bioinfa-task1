#!/usr/bin/env python3
"""
Инструменты (реализованы через BioPython PairwiseAligner):
  1. Needle  — глобальное  выравнивание (Needleman-Wunsch, NW)
  2. Water   — локальное   выравнивание (Smith-Waterman, SW)
"""

import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

SEQ_DIR   = "sequences"
ALN_DIR   = "alignments"
os.makedirs(ALN_DIR, exist_ok=True)

MATRIX    = substitution_matrices.load("BLOSUM62")
GAP_OPEN  = -10.0     # как в EMBOSS Needle/Water
GAP_EXT   = -0.5

def load_seq(name: str):
    path = os.path.join(SEQ_DIR, f"{name}.fasta")
    return SeqIO.read(path, "fasta")


def make_aligner(mode: str) -> PairwiseAligner:
    aln = PairwiseAligner()
    aln.mode = mode 
    aln.substitution_matrix = MATRIX
    aln.open_gap_score = GAP_OPEN
    aln.extend_gap_score = GAP_EXT
    return aln


def calc_stats(alignment) -> dict:
    fmt_lines = alignment.format().split("\n")
    match_chars = []
    i = 0
    while i < len(fmt_lines):
        line = fmt_lines[i]
        if line.startswith("target"):
            match_line = fmt_lines[i + 1] if i + 1 < len(fmt_lines) else ""
            prefix_len = line.index(line.split()[1]) if len(line.split()) >= 2 else 0
            parts = line.split()
            if len(parts) >= 3:
                seq_start = line.index(parts[2])
                match_seq = match_line[seq_start:seq_start + len(parts[2])]
                match_chars.append(match_seq)
            i += 3
        else:
            i += 1

    match_str = "".join(match_chars)
    aln_len     = len(match_str)
    identities  = match_str.count("|")
    similarities = identities + match_str.count(".")
    gaps        = match_str.count("-")

    return {
        "score":        alignment.score,
        "aln_len":      aln_len,
        "identities":   identities,
        "similarities": similarities,
        "gaps":         gaps,
        "pct_identity": 100.0 * identities   / aln_len if aln_len else 0,
        "pct_similar":  100.0 * similarities / aln_len if aln_len else 0,
        "pct_gaps":     100.0 * gaps          / aln_len if aln_len else 0,
    }


def save_alignment(alignment, path: str):
    """Сохранить выравнивание в файл."""
    with open(path, "w") as f:
        f.write(str(alignment))
    print(f"  Сохранено: {path}")


def print_stats(label: str, stats: dict):
    print(f"\n  {'─'*50}")
    print(f"  {label}")
    print(f"  {'─'*50}")
    print(f"  Score       : {stats['score']:.1f}")
    print(f"  Длина ALN   : {stats['aln_len']}")
    print(f"  Идентичность: {stats['identities']}/{stats['aln_len']}  "
          f"({stats['pct_identity']:.1f}%)")
    print(f"  Схожесть    : {stats['similarities']}/{stats['aln_len']}  "
          f"({stats['pct_similar']:.1f}%)")
    print(f"  Гэпы        : {stats['gaps']}/{stats['aln_len']}  "
          f"({stats['pct_gaps']:.1f}%)")

PAIRS = [
    ("PAH_human",  "PAH_mouse",  "PAH"),
    ("QDPR_human", "QDPR_mouse", "QDPR"),
]

all_results = {}

for seq_human_name, seq_mouse_name, gene in PAIRS:
    print(f"\n{'='*60}")
    print(f"  Ген: {gene}")
    print(f"  Человек : {seq_human_name}")
    print(f"  Мышь    : {seq_mouse_name}")
    print(f"{'='*60}")

    rec_h = load_seq(seq_human_name)
    rec_m = load_seq(seq_mouse_name)

    print(f"  Длина (человек): {len(rec_h.seq)} aa")
    print(f"  Длина (мышь)   : {len(rec_m.seq)} aa")

    gene_results = {}

    for mode, tool_name in [("global", "Needle (NW)"), ("local", "Water (SW)")]:
        aligner = make_aligner(mode)
        alignments = aligner.align(rec_h.seq, rec_m.seq)
        best = next(iter(alignments))

        out_path = os.path.join(ALN_DIR, f"{gene}_{mode}.txt")
        save_alignment(best, out_path)

        stats = calc_stats(best)
        gene_results[mode] = stats
        print_stats(f"{tool_name}  ({mode})", stats)

    all_results[gene] = gene_results

print(f"\n\n{'='*75}")
print(f"{'СВОДНАЯ ТАБЛИЦА ВЫРАВНИВАНИЙ':^75}")
print(f"{'='*75}")
header = f"{'Ген':<6} {'Инструмент':<16} {'Score':>8} {'ALN':>6} {'%ID':>7} {'%SIM':>7} {'%GAP':>7}"
print(header)
print(f"{'─'*75}")

for gene, results in all_results.items():
    for mode, tool_name in [("global", "Needle (NW)"), ("local", "Water (SW)")]:
        s = results[mode]
        print(f"{gene:<6} {tool_name:<16} {s['score']:>8.1f} {s['aln_len']:>6} "
              f"{s['pct_identity']:>7.1f} {s['pct_similar']:>7.1f} {s['pct_gaps']:>7.1f}")

print(f"{'─'*75}")
print("\nПодробные выравнивания сохранены в папке alignments/")
print("Done.")

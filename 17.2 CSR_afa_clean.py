#!/usr/bin/env python3
"""
remove_gap_N_codons.py

Walks through every CSR alignment (.afa) under
  sequences/CDS_CSR/<Genus>/aligned/*.afa
and in each alignment:

  • Finds the shortest sequence length (trimming any overhangs).
  • Ensures length is a multiple of 3 (drops trailing 1–2 nt if needed).
  • Builds a mask of “good” codon positions (no sequence has ‘-’ or ‘N’ in that 3-nt window).
  • Rewrites the alignment in place, keeping only those clean codons across all samples.

Usage:
  python remove_gap_N_codons.py
"""

import os
from pathlib import Path
from Bio import SeqIO

# base directory for your CSR alignments
LOCAL_BASE = Path("sequences/CDS_CSR")

def clean_alignment_file(afa_path: Path):
    # parse all records
    records = list(SeqIO.parse(str(afa_path), "fasta"))
    if not records:
        print(f"[SKIP] {afa_path} — no sequences found.")
        return

    # dict of id -> sequence string
    seqs = {rec.id: str(rec.seq) for rec in records}

    # 1) uniform length: truncate to shortest
    min_len = min(len(s) for s in seqs.values())
    for k in seqs:
        seqs[k] = seqs[k][:min_len]

    # 2) make length a multiple of 3
    remainder = min_len % 3
    if remainder:
        min_len -= remainder
        for k in seqs:
            seqs[k] = seqs[k][:min_len]

    if min_len == 0:
        print(f"[WARN] {afa_path} truncated to zero length → clearing file.")
        afa_path.write_text("")  # wipe out
        return

    num_codons = min_len // 3

    # 3) build mask: True if codon OK in ALL sequences
    codon_mask = [True] * num_codons
    for i in range(num_codons):
        start, end = 3*i, 3*i + 3
        for s in seqs.values():
            codon = s[start:end]
            if "-" in codon or "N" in codon:
                codon_mask[i] = False
                break

    # 4) reconstruct each sequence using only passing codons
    cleaned = {}
    for id_, full in seqs.items():
        parts = []
        for i in range(num_codons):
            if codon_mask[i]:
                parts.append(full[3*i:3*i+3])
        cleaned[id_] = "".join(parts)

    # 5) write back in FASTA format (overwrite)
    with afa_path.open("w") as out:
        for id_, seq in cleaned.items():
            out.write(f">{id_}\n{seq}\n")

    kept = sum(codon_mask)
    dropped = num_codons - kept
    print(f"[CLEANED] {afa_path.name}: {kept}/{num_codons} codons kept ({dropped} dropped)")

def main():
    for genus_dir in LOCAL_BASE.iterdir():
        if not genus_dir.is_dir():
            continue
        aligned_dir = genus_dir / "aligned"
        if not aligned_dir.is_dir():
            continue

        for afa_file in aligned_dir.glob("*.afa"):
            clean_alignment_file(afa_file)

if __name__ == "__main__":
    main()

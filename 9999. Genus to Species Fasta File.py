#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Split genus-level concatenated FASTAs into per-species FASTAs
and print per-species / overall sequence counts.

Examples of inputs (read-only):
  sequences/CDS_Genus/Balaena/CDS_nucleotide_gapped/Balaena_concatenated_gapped_sequences.fasta
  sequences/CDS_Genus/Balaenoptera/CDS_nucleotide_gapped/Balaenoptera_concatenated_gapped_sequences.fasta

Outputs created:
  sequences/CDS_Species/<Genus_species>/<Genus_species>_concatenated_gapped_sequences.fasta

Usage (from your repo root):
  python 9999_genus_to_species_fasta.py
  # with options
  python 9999_genus_to_species_fasta.py --genus-root sequences/CDS_Genus \
                                        --out-root sequences/CDS_Species \
                                        --csv-summary sequences/CDS_Species/_species_counts.csv \
                                        --show-skipped-examples 5
"""

from pathlib import Path
import re
import sys
import csv
import argparse
from collections import defaultdict

# ---------------- FASTA I/O ----------------
try:
    from Bio import SeqIO  # type: ignore
    BIOPYTHON = True
except Exception:
    BIOPYTHON = False


def iter_fasta_records(fpath: Path):
    """Yield (header, sequence) tuples from a FASTA file."""
    if BIOPYTHON:
        for rec in SeqIO.parse(str(fpath), "fasta"):
            # Prefer full description; strip leading '>'
            header = rec.description.strip()
            if header.startswith(">"):
                header = header[1:]
            seq = str(rec.seq).replace("\n", "").strip()
            yield header, seq
    else:
        header = None
        chunks = []
        with open(fpath, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(chunks)
                    header = line[1:].strip()
                    chunks = []
                else:
                    chunks.append(line.strip())
            if header is not None:
                yield header, "".join(chunks)


def wrap_fasta_seq(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# ---------------- Species extraction ----------------
# NOTE: No \b (word boundaries) because underscores are "word" chars in Python's regex.
# Match anywhere: "Genus_species", "Genus species", or "Genus-species"
BINOMIAL_RE = re.compile(r"([A-Z][a-z]+)[ _-]([a-z]+)")


def GENUS_FIRST_RE(genus: str):
    """Regex that prefers the folder's genus first (case-sensitive)."""
    return re.compile(rf"{re.escape(genus)}[ _-]([a-z]+)")


def extract_species_name(header: str, genus_hint: str | None = None) -> str | None:
    """
    Extract 'Genus_species' from a FASTA header.
    Prefer a match whose Genus == genus_hint (if provided).
    """
    h = header.strip()

    if genus_hint:
        m = GENUS_FIRST_RE(genus_hint).search(h)
        if m:
            return f"{genus_hint}_{m.group(1)}"

    m = BINOMIAL_RE.search(h)
    if m:
        return f"{m.group(1)}_{m.group(2)}"

    return None


# ---------------- File ops ----------------
def find_genus_concat_fastas(root: Path):
    """Yield paths to *_concatenated_gapped_sequences.fasta under **/CDS_nucleotide_gapped/."""
    pattern = "**/CDS_nucleotide_gapped/*_concatenated_gapped_sequences.fasta"
    return sorted(root.glob(pattern))


def write_species_fasta(out_root: Path, species: str, records: list[tuple[str, str]], overwrite: bool):
    """
    Write a multi-FASTA for a species to:
      {out_root}/{species}/{species}_concatenated_gapped_sequences.fasta
    """
    out_dir = out_root / species
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fa = out_dir / f"{species}_concatenated_gapped_sequences.fasta"

    mode = "w" if overwrite else ("a" if out_fa.exists() else "w")
    with open(out_fa, mode, encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            fh.write(wrap_fasta_seq(seq) + "\n")
    return out_fa


def process_genus_fasta(fpath: Path, show_skipped_examples: int = 0):
    """
    Split one genus-level FASTA into per-species collections.
    Returns (species_bins, unknown_count, genus_hint, skipped_examples)
    """
    try:
        # parent=CDS_nucleotide_gapped, parent.parent=<Genus>
        genus_hint = fpath.parent.parent.name
    except Exception:
        genus_hint = None

    species_bins: dict[str, list[tuple[str, str]]] = defaultdict(list)
    unknown = 0
    skipped_examples: list[str] = []

    for header, seq in iter_fasta_records(fpath):
        sp = extract_species_name(header, genus_hint=genus_hint)
        if sp is None:
            unknown += 1
            if show_skipped_examples and len(skipped_examples) < show_skipped_examples:
                skipped_examples.append(header[:200])
            continue
        species_bins[sp].append((header, seq))

    return species_bins, unknown, genus_hint, skipped_examples


def print_summary(species_counts: dict[str, int], total_unknown: int):
    print("\n[SUMMARY] Per-species sequence counts")
    if not species_counts:
        print("  (No sequences assigned to any species.)")
    else:
        items = sorted(species_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        name_w = max(len(k) for k in species_counts.keys())
        count_w = max(len(str(v)) for v in species_counts.values())
        for sp, n in items:
            print(f"  {sp:<{name_w}} : {n:>{count_w}}")

        total_assigned = sum(species_counts.values())
        print("\n[TOTALS]")
        print(f"  Unique species                  : {len(species_counts)}")
        print(f"  Total sequences (assigned)      : {total_assigned}")
        print(f"  Total records skipped (unknown) : {total_unknown}")


def write_csv_summary(csv_path: Path, species_counts: dict[str, int], total_unknown: int):
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    rows = [["Species", "Count"]] + sorted(species_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(rows)
        writer.writerow([])
        writer.writerow(["Total unique species", len(species_counts)])
        writer.writerow(["Total sequences (assigned)", sum(species_counts.values())])
        writer.writerow(["Total records skipped (unknown)", total_unknown])
    print(f"[INFO] Wrote CSV summary → {csv_path}")


# ---------------- Main ----------------
def parse_args():
    p = argparse.ArgumentParser(description="Split genus-level concatenated FASTAs into per-species FASTAs.")
    p.add_argument("--genus-root", type=Path, default=Path("sequences/CDS_Genus"),
                   help="Input root containing <Genus>/CDS_nucleotide_gapped/*.fasta (default: sequences/CDS_Genus)")
    p.add_argument("--out-root", type=Path, default=Path("sequences/CDS_Species"),
                   help="Output root for per-species FASTAs (default: sequences/CDS_Species)")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing per-species files (default if neither overwrite nor append is set).")
    p.add_argument("--append", action="store_true",
                   help="Append to existing per-species files instead of overwriting.")
    p.add_argument("--show-skipped-examples", type=int, default=0,
                   help="Print up to N example headers that were skipped (no species detected).")
    p.add_argument("--csv-summary", type=Path, default=None,
                   help="Optional path to write a CSV summary of species counts.")
    return p.parse_args()


def main():
    args = parse_args()

    if args.overwrite and args.append:
        print("[ERROR] --overwrite and --append are mutually exclusive.", file=sys.stderr)
        sys.exit(2)

    overwrite = True
    if args.append:
        overwrite = False
    elif args.overwrite:
        overwrite = True

    genus_root: Path = args.genus_root
    out_root: Path = args.out_root

    if not genus_root.exists():
        print(f"[ERROR] Input root not found: {genus_root.resolve()}", file=sys.stderr)
        sys.exit(1)

    out_root.mkdir(parents=True, exist_ok=True)

    genus_fastas = find_genus_concat_fastas(genus_root)
    if not genus_fastas:
        print(f"[WARN] No genus-level concatenated FASTAs found under {genus_root}/**/CDS_nucleotide_gapped/")
        sys.exit(0)

    total_unknown = 0
    species_counts: dict[str, int] = defaultdict(int)
    total_species_files_written = 0

    for fpath in genus_fastas:
        species_bins, unknown, genus_hint, skipped_examples = process_genus_fasta(
            fpath, show_skipped_examples=args.show_skipped_examples
        )
        total_unknown += unknown

        print(f"\n[INFO] Processing: {fpath}")
        if genus_hint:
            print(f"       Genus hint: {genus_hint}")
        if unknown:
            print(f"       [WARN] {unknown} record(s) had no detectable species binomial and were skipped.")
            if skipped_examples:
                print("       Examples of skipped headers:")
                for ex in skipped_examples:
                    print(f"         - {ex}")

        for species, recs in sorted(species_bins.items()):
            out_fa = write_species_fasta(out_root, species, recs, overwrite=overwrite)
            total_species_files_written += 1
            species_counts[species] += len(recs)
            print(f"       Wrote {len(recs):>3} seq(s) → {out_fa}")

    print_summary(species_counts, total_unknown)
    print(f"\n[DONE] Created/updated {total_species_files_written} species FASTA file(s) under {out_root}.")

    if args.csv_summary:
        write_csv_summary(args.csv_summary, species_counts, total_unknown)


if __name__ == "__main__":
    main()

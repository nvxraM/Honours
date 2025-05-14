#!/usr/bin/env python3
"""
compute_zscores.py

1) Read your closest-relative map.
2) Load the genus-level FASTA for each species’ genus.
3) From that FASTA, pull only the records whose header endswith _Genus_species.
4) Count mismatches vs. the reference sequence of the closest relative,
   SKIPPING any records whose length doesn’t match.
5) Average those counts and write out an Excel.

Usage:
  pip install biopython pandas openpyxl
  python compute_zscores.py
"""

import os, sys
from Bio import SeqIO
import pandas as pd

# 1) species → closest relative (all 33)
closest_relative = {
    "Balaena mysticetus":    "Eubalaena japonica",
    "Balaenoptera musculus": "Balaenoptera physalus",
    "Balaenoptera physalus": "Balaenoptera musculus",
    "Delphinapterus leucas": "Peponocephala electra",
    "Delphinus delphis":     "Tursiops aduncus",
    "Eubalaena australis":   "Eubalaena glacialis",
    "Eubalaena glacialis":   "Eubalaena australis",
    "Eubalaena japonica":    "Eubalaena australis",
    "Globicephala macrorhynchus": "Pseudorca crassidens",
    "Hyperoodon ampullatus":       "Ziphius cavirostris",
    "Mesoplodon densirostris":     "Mesoplodon grayi",
    "Mesoplodon europaeus":        "Mesoplodon mirus",
    "Mesoplodon grayi":            "Mesoplodon densirostris",
    "Mesoplodon mirus":            "Mesoplodon europaeus",
    "Monodon monoceros":           "Delphinapterus leucas",
    "Neophocaena asiaeorientalis": "Neophocaena phocaenoides",
    "Neophocaena phocaenoides":    "Neophocaena asiaeorientalis",
    "Orcaella brevirostris":       "Orcinus orca",
    "Orcinus orca":                "Orcaella brevirostris",
    "Peponocephala electra":       "Delphinapterus leucas",
    "Phocoena phocoena":           "Phocoenoides dalli",
    "Phocoena sinus":              "Phocoena spinipinnis",
    "Phocoena spinipinnis":        "Phocoena sinus",
    "Phocoenoides dalli":          "Phocoena phocoena",
    "Physeter macrocephalus":      "Platanista gangetica",
    "Platanista gangetica":        "Physeter macrocephalus",
    "Pseudorca crassidens":        "Globicephala macrorhynchus",
    "Stenella attenuata":          "Stenella longirostris",
    "Stenella longirostris":       "Stenella attenuata",
    "Tursiops aduncus":            "Delphinus delphis",
    "Tursiops australis":          "Tursiops truncatus",
    "Tursiops truncatus":          "Tursiops australis",
    "Ziphius cavirostris":         "Hyperoodon ampullatus",
}

# 2) Paths
REF_FASTA  = "sequences/CDS_Reference/CDS_nucleotide_gapped/" \
             "CDS_Reference_concatenated_gapped_sequences.fasta"
GENUS_DIR  = "sequences/CDS_Genus"
OUTPUT_XLS = "sequences/zscores/scores.xlsx"


def load_reference_sequences(path):
    """
    Parse the reference FASTA into a dict: "Genus species" → sequence.
    Expects headers like >NC_XXXXX_Genus_species
    """
    refs = {}
    for rec in SeqIO.parse(path, "fasta"):
        parts = rec.id.split("_")
        if len(parts) >= 4 and parts[0] == "NC":
            species = f"{parts[2]} {parts[3]}"
        else:
            raise ValueError(f"Can't parse species from ref header: {rec.id}")
        refs[species] = str(rec.seq)
    return refs


def compute_average_mismatches_for_species(species, ref_seq, genus_fasta):
    """
    From genus_fasta, pick only records ending in _Genus_species,
    skip any whose length != len(ref_seq),
    tally mismatches vs. ref_seq, and return the average.
    """
    # suffix to identify only our focal species
    sp_und = species.replace(" ", "_")
    want_suffix = f"_{sp_und}"

    all_recs = list(SeqIO.parse(genus_fasta, "fasta"))
    # filter to only those that match species suffix
    recs = [r for r in all_recs if r.id.endswith(want_suffix)]
    if not recs:
        raise RuntimeError(f"No records for {species} in {genus_fasta}")

    L = len(ref_seq)
    scores = []
    skipped = 0

    for rec in recs:
        seq = str(rec.seq)
        if len(seq) != L:
            skipped += 1
            print(f"  [WARN] skipping {rec.id} (len {len(seq)} ≠ {L})")
            continue
        diff = sum(1 for a, b in zip(seq, ref_seq) if a != b)
        scores.append(diff)

    if not scores:
        raise RuntimeError(f"All {len(recs)} records for {species} were wrong length.")
    if skipped:
        print(f"  [INFO] {skipped} of {len(recs)} records skipped for {species}")

    return sum(scores) / len(scores)


def main():
    # load references
    if not os.path.exists(REF_FASTA):
        sys.exit(f"Ref FASTA missing: {REF_FASTA}")
    ref_dict = load_reference_sequences(REF_FASTA)

    # ensure output dir
    os.makedirs(os.path.dirname(OUTPUT_XLS), exist_ok=True)

    # compute per species
    results = []
    for sp, rel in closest_relative.items():
        if rel not in ref_dict:
            sys.exit(f"Ref for '{rel}' not in {REF_FASTA}")

        genus = sp.split()[0]
        genus_fasta = os.path.join(
            GENUS_DIR, genus,
            "CDS_nucleotide_gapped",
            f"{genus}_concatenated_gapped_sequences.fasta"
        )
        if not os.path.exists(genus_fasta):
            sys.exit(f"Genus FASTA not found: {genus_fasta}")

        print(f"\nProcessing {sp}  (ref = {rel})")
        avg_score = compute_average_mismatches_for_species(
            sp, ref_dict[rel], genus_fasta
        )
        print(f" → avg mismatches = {avg_score:.1f}")
        results.append({
            "Species": sp,
            "Score":   avg_score,
            "Closest_relative": rel
        })

    # write Excel
    df = pd.DataFrame(results)
    df.to_excel(OUTPUT_XLS, index=False)
    print(f"\nDone! Written to {OUTPUT_XLS}")


if __name__ == "__main__":
    main()

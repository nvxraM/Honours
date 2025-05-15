#!/usr/bin/env python3
"""
generate_csr_fastas.py

Script to generate CSR FASTA files combining each species' multiple gapped sequences
with the reference sequence of its closest relative.

Directory structure assumed under `sequences/`:

sequences/
  CDS_Genus/
    <Genus>/
      CDS_nucleotide_gapped/
        <Genus>_concatenated_gapped_sequences.fasta
  CDS_Reference/
    CDS_nucleotide_gapped/
      CDS_Reference_concatenated_gapped_sequences.fasta
  CDS_CSR/  # output will be placed here

Requires Biopython.
"""
import os
from pathlib import Path
from Bio import SeqIO

# ---------------------- Configuration ----------------------
BASE_DIR = Path("sequences")
GENUS_DIR = BASE_DIR / "CDS_Genus"
REF_FASTA = BASE_DIR / "CDS_Reference" / "CDS_nucleotide_gapped" / "CDS_Reference_concatenated_gapped_sequences.fasta"
CSR_DIR = BASE_DIR / "CDS_CSR"

# Mapping from species to its closest relative (reference)
closest_relatives = {
    "Balaena_mysticetus": "Eubalaena_japonica",
    "Balaenoptera_musculus": "Balaenoptera_physalus",
    "Balaenoptera_physalus": "Balaenoptera_musculus",
    "Delphinapterus_leucas": "Peponocephala_electra",
    "Delphinus_delphis": "Tursiops_aduncus",
    "Eubalaena_australis": "Eubalaena_glacialis",
    "Eubalaena_glacialis": "Eubalaena_australis",
    "Eubalaena_japonica": "Eubalaena_australis",
    "Globicephala_macrorhynchus": "Pseudorca_crassidens",
    "Hyperoodon_ampullatus": "Ziphius_cavirostris",
    "Mesoplodon_densirostris": "Mesoplodon_grayi",
    "Mesoplodon_europaeus": "Mesoplodon_mirus",
    "Mesoplodon_grayi": "Mesoplodon_densirostris",
    "Mesoplodon_mirus": "Mesoplodon_europaeus",
    "Monodon_monoceros": "Delphinapterus_leucas",
    "Neophocaena_asiaeorientalis": "Neophocaena_phocaenoides",
    "Neophocaena_phocaenoides": "Neophocaena_asiaeorientalis",
    "Orcaella_brevirostris": "Orcinus_orca",
    "Orcinus_orca": "Orcaella_brevirostris",
    "Peponocephala_electra": "Delphinapterus_leucas",
    "Phocoena_phocoena": "Phocoenoides_dalli",
    "Phocoena_sinus": "Phocoena_spinipinnis",
    "Phocoena_spinipinnis": "Phocoena_sinus",
    "Phocoenoides_dalli": "Phocoena_phocoena",
    "Physeter_macrocephalus": "Platanista_gangetica",
    "Platanista_gangetica": "Physeter_macrocephalus",
    "Pseudorca_crassidens": "Globicephala_macrorhynchus",
    "Stenella_attenuata": "Stenella_longirostris",
    "Stenella_longirostris": "Stenella_attenuata",
    "Tursiops_aduncus": "Delphinus_delphis",
    "Tursiops_australis": "Tursiops_truncatus",
    "Tursiops_truncatus": "Tursiops_australis",
    "Ziphius_cavirostris": "Hyperoodon_ampullatus",
}
# -----------------------------------------------------------

def load_ref_by_species(ref_path):
    """
    Parse the reference FASTA and map by species key (last two name fields).
    """
    ref_by_species = {}
    for rec in SeqIO.parse(ref_path, "fasta"):
        parts = rec.id.split("_")
        key = "_".join(parts[-2:])
        ref_by_species[key] = rec
    return ref_by_species


def main():
    # Ensure output folder exists
    CSR_DIR.mkdir(parents=True, exist_ok=True)

    # Load reference sequences
    if not REF_FASTA.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {REF_FASTA}")
    ref_by_species = load_ref_by_species(str(REF_FASTA))

    total_individuals = 0
    # Iterate each genus
    for genus_path in GENUS_DIR.iterdir():
        if not genus_path.is_dir():
            continue
        genus = genus_path.name
        gapped_fasta = genus_path / "CDS_nucleotide_gapped" / f"{genus}_concatenated_gapped_sequences.fasta"
        if not gapped_fasta.exists():
            print(f"[WARN] Missing gapped file for genus {genus}: {gapped_fasta}")
            continue

        # Read all records for this genus
        records = list(SeqIO.parse(str(gapped_fasta), "fasta"))

        # Group by species
        species_groups = {}
        for rec in records:
            parts = rec.id.split("_")
            species = "_".join(parts[-2:])
            species_groups.setdefault(species, []).append(rec)

        # Create CSR for each species
        out_dir = CSR_DIR / genus
        out_dir.mkdir(parents=True, exist_ok=True)

        for species, seq_recs in species_groups.items():
            closest = closest_relatives.get(species)
            if not closest:
                print(f"[WARN] No mapping for {species}")
                continue
            ref_rec = ref_by_species.get(closest)
            if not ref_rec:
                print(f"[ERROR] Reference for {closest} not found")
                continue

            out_file = out_dir / f"{genus}_{species}.fasta"
            with open(out_file, "w") as fh:
                # Write full reference header first
                SeqIO.write(ref_rec, fh, "fasta")
                # Then all sample sequences
                SeqIO.write(seq_recs, fh, "fasta")
            # Log with count and CSR info
            count = len(seq_recs)
            total_individuals += count + 1
            print(f"[INFO] Written CSR for {species}: {out_file} (Number of individuals: {count}, CSR is {closest})")

    # After processing all species
    print(f"[INFO] Total individuals across all species: {total_individuals}")

if __name__ == "__main__":
    main()

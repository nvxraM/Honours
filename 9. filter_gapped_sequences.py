import os
from Bio import SeqIO


def enforce_minimum_sequences(base_dir, min_sequences=5):
    """
    For each species directory under the given base directory:
    1. Identify the concatenated gapped sequences FASTA file.
    2. Count the number of sequences in that file.
    3. If there are fewer than `min_sequences` sequences, move all sequences to
       the excluded FASTA file.

    Directory structure assumed:
    sequences/
      ├─ CDS/
      │  ├─ [Species]/
      │  │  ├─ CDS_nucleotide_gapped/
      │  │  │  └─ [Species]_concatenated_gapped_sequences.fasta
      ├─ excluded/
      │  └─ CDS/
      │     └─ [Species]/
      │        └─ excluded_[Species]_concatenated_gapped_sequences.fasta

    Args:
        base_dir (str): The base directory containing 'CDS' and 'excluded' subdirectories.
        min_sequences (int): The minimum number of sequences required to remain in the file.
                             If fewer, they are all moved to the excluded directory.
    """
    # Path to the base CDS directory
    cds_path = os.path.join(base_dir, "CDS")

    # List all directories in the CDS directory that represent species
    species_list = [
        species for species in os.listdir(cds_path)
        if os.path.isdir(os.path.join(cds_path, species))
    ]

    for species in species_list:
        # Paths for input and output
        gapped_dir = os.path.join(cds_path, species, "CDS_nucleotide_gapped")
        input_file = os.path.join(gapped_dir, f"{species}_concatenated_gapped_sequences.fasta")

        # If the input file does not exist, skip this species
        if not os.path.exists(input_file):
            print(f"No concatenated FASTA found for {species}, skipping.")
            continue

        # Read all sequences
        records = list(SeqIO.parse(input_file, "fasta"))
        num_sequences = len(records)

        # If fewer than the minimum required sequences, exclude all of them
        if num_sequences < min_sequences:
            excluded_dir = os.path.join(base_dir, "excluded", "CDS", species)
            os.makedirs(excluded_dir, exist_ok=True)
            excluded_file = os.path.join(excluded_dir, f"excluded_{species}_concatenated_gapped_sequences.fasta")

            # Move all sequences to the excluded file
            with open(excluded_file, "w") as out_excluded:
                SeqIO.write(records, out_excluded, "fasta")

            # Remove the original file or replace it with an empty file
            # Since the user might want to know that we had a file, we'll just empty it.
            with open(input_file, "w") as out:
                pass  # Create/overwrite with an empty file

            print(f"Fewer than {min_sequences} sequences found for {species}, moved all to excluded.")
        else:
            print(f"{species} has {num_sequences} sequences, which meets the minimum requirement.")


# Base directory containing the sequences
base_directory = "sequences"

# Enforce minimum number of sequences (5)
enforce_minimum_sequences(base_directory, min_sequences=5)

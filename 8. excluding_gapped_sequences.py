import os
from Bio import SeqIO


def exclude_gapped_sequences(base_dir):
    """
    For each species directory under the given base directory:
    1. Identify the concatenated gapped sequences FASTA file.
    2. Check each sequence in that file for gaps ('-').
    3. Write sequences containing gaps to an excluded FASTA file.
    4. Write sequences without gaps back into the original concatenated FASTA file.

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

        # Prepare the excluded directory and file
        excluded_dir = os.path.join(base_dir, "excluded", "CDS", species)
        os.makedirs(excluded_dir, exist_ok=True)
        excluded_file = os.path.join(excluded_dir, f"excluded_{species}_concatenated_gapped_sequences.fasta")

        # Read all sequences and classify them into gapped or gapless
        records = list(SeqIO.parse(input_file, "fasta"))
        gapless_records = []
        gapped_records = []

        for record in records:
            if "-" in str(record.seq):
                gapped_records.append(record)
            else:
                gapless_records.append(record)

        # Write gapless sequences back to the main concatenated file
        with open(input_file, "w") as out:
            SeqIO.write(gapless_records, out, "fasta")

        # Write gapped sequences to the excluded file
        if gapped_records:
            with open(excluded_file, "w") as out_excluded:
                SeqIO.write(gapped_records, out_excluded, "fasta")
            print(f"Excluded {len(gapped_records)} gapped sequence(s) for {species}.")
        else:
            print(f"No gapped sequences found for {species}.")


# Base directory containing the sequences
base_directory = "sequences"

# Run the exclusion process
exclude_gapped_sequences(base_directory)

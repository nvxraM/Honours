import os
from Bio import SeqIO

# Original gene order including ND6
GENE_ORDER = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"
]

# Set this to True to exclude ND6 from the concatenation.
# Change this to False if you want to include ND6 again.
EXCLUDE_ND6 = True


def get_filtered_gene_order():
    """
    Return the gene order, potentially excluding ND6 if EXCLUDE_ND6 is set to True.
    If EXCLUDE_ND6 is False, ND6 will be included.

    Returns:
        list: A list of gene names in the order they should be concatenated.
    """
    if EXCLUDE_ND6:
        return [gene for gene in GENE_ORDER if gene != "ND6"]
    else:
        return GENE_ORDER


def concatenate_gapped_sequences(base_dir, species):
    """
    Concatenate the nucleotide sequences of gapped genes for the given species
    into a single FASTA file. ND6 is excluded if EXCLUDE_ND6 is True.

    This function:
    1. Determines the order of genes to be included (with or without ND6).
    2. Iterates over each gene FASTA file within the species directory.
    3. Concatenates sequences for each accession across all provided genes.
    4. Writes the concatenated sequences to a single FASTA output file.

    Args:
        base_dir (str): The base directory containing species subdirectories.
        species (str): The name of the species subdirectory.
    """
    # Input and output directories for the given species
    input_dir = os.path.join(base_dir, species, "CDS_nucleotide_gapped")
    output_dir = os.path.join(base_dir, species, "CDS_nucleotide_gapped")
    os.makedirs(output_dir, exist_ok=True)

    # Output file path for concatenated nucleotide sequences
    output_file = os.path.join(output_dir, f"{species}_concatenated_gapped_sequences.fasta")

    # Dictionary to hold concatenated sequences keyed by accession
    concatenated_sequences = {}

    # Get the filtered gene order list
    filtered_genes = get_filtered_gene_order()

    # Process each gene in the filtered gene set
    for gene in filtered_genes:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")

        # Warn if the gene file is missing and skip it
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue

        # Parse the FASTA file for the current gene
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id

                # Remove the appended gene name from the accession if present
                # For example, if the accession is "ABC123_ND2", remove "_ND2"
                for gene_name in GENE_ORDER:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break

                # Initialize the accession key if it's the first time we see it
                if accession not in concatenated_sequences:
                    concatenated_sequences[accession] = ""

                # Append the current gene sequence to the accession's sequence
                concatenated_sequences[accession] += str(record.seq)

    # Write out the concatenated sequences to the output file
    with open(output_file, "w") as output:
        for accession, sequence in concatenated_sequences.items():
            output.write(f">{accession}\n{sequence}\n")


# Base directory containing the species subdirectories
base_directory = "sequences/CDS"

# Identify all directories in the base directory that represent species
species_list = [
    species for species in os.listdir(base_directory)
    if os.path.isdir(os.path.join(base_directory, species))
]

# Concatenate gapped sequences for each species
for species in species_list:
    concatenate_gapped_sequences(base_directory, species)

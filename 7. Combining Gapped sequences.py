import os
from Bio import SeqIO

# Original gene order including ND6
gene_order = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
              "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"]

# Set this to True to exclude ND6 from the concatenation.
# Change this to False if you want to include ND6 again.
exclude_nd6 = True


def get_filtered_gene_order():
    """
    Return the gene order, potentially excluding ND6 if exclude_nd6 is set to True.
    If exclude_nd6 is False, ND6 will be included.
    """
    if exclude_nd6:
        return [gene for gene in gene_order if gene != "ND6"]
    else:
        return gene_order


def concatenate_sequences(base_dir, species):
    """
    Concatenate the nucleotide sequences of genes for the given species
    into a single FASTA file. ND6 is excluded if exclude_nd6 is True.
    """
    # Input and output directories
    input_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    output_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    os.makedirs(output_dir, exist_ok=True)

    # Output file path for concatenated nucleotide sequences
    output_file = os.path.join(output_dir, f"{species}_concatenated_sequences.fasta")

    # Dictionary to hold concatenated sequences keyed by accession
    concatenated_sequences = {}

    # Get the filtered gene order list (without ND6 if exclude_nd6=True)
    filtered_genes = get_filtered_gene_order()

    for gene in filtered_genes:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")

        # If the gene file does not exist, print a warning and skip it
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue

        # Parse each gene FASTA file
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id

                # Remove the appended gene name from the accession if present
                for gene_name in gene_order:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break

                # If this accession hasn't been encountered, initialize it
                if accession not in concatenated_sequences:
                    concatenated_sequences[accession] = ""

                # Append the current gene sequence to the accession's sequence
                concatenated_sequences[accession] += str(record.seq)

    # Write out the concatenated sequences to the output file
    with open(output_file, "w") as output:
        for accession, sequence in concatenated_sequences.items():
            output.write(f">{accession}\n{sequence}\n")


def concatenate_proteins(base_dir, species):
    """
    Concatenate the protein sequences of genes for the given species
    into a single FASTA file. ND6 is excluded if exclude_nd6 is True.
    """
    # Input and output directories
    input_dir = os.path.join(base_dir, species, "CDS_protein")
    output_dir = os.path.join(base_dir, species, "CDS_protein")
    os.makedirs(output_dir, exist_ok=True)

    # Output file path for concatenated protein sequences
    output_file = os.path.join(output_dir, f"{species}_concatenated_proteins.fasta")

    # Dictionary to hold concatenated protein sequences keyed by accession
    concatenated_proteins = {}

    # Get the filtered gene order list (without ND6 if exclude_nd6=True)
    filtered_genes = get_filtered_gene_order()

    for gene in filtered_genes:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")

        # If the gene file does not exist, print a warning and skip it
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue

        # Parse each gene FASTA file
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id

                # Remove the appended gene name from the accession if present
                for gene_name in gene_order:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break

                # If this accession hasn't been encountered, initialize it
                if accession not in concatenated_proteins:
                    concatenated_proteins[accession] = ""

                # Append the current protein sequence to the accession's sequence
                concatenated_proteins[accession] += str(record.seq)

    # Write out the concatenated protein sequences to the output file
    with open(output_file, "w") as output:
        for accession, sequence in concatenated_proteins.items():
            output.write(f">{accession}\n{sequence}\n")


# Base directory containing the species subdirectories
base_directory = "sequences/CDS"

# List all directories in the base directory that represent species
species_list = [
    species for species in os.listdir(base_directory)
    if os.path.isdir(os.path.join(base_directory, species))
]

# Concatenate sequences and proteins for each species
for species in species_list:
    concatenate_sequences(base_directory, species)
    concatenate_proteins(base_directory, species)

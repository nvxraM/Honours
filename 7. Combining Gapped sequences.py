import os
from Bio import SeqIO

gene_order = ["ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"]

def concatenate_sequences(base_dir, species):
    input_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    output_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{species}_concatenated_sequences.fasta")

    concatenated_sequences = {}

    for gene in gene_order:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id
                for gene_name in gene_order:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break
                if accession not in concatenated_sequences:
                    concatenated_sequences[accession] = ""
                concatenated_sequences[accession] += str(record.seq)

    with open(output_file, "w") as output:
        for accession, sequence in concatenated_sequences.items():
            output.write(f">{accession}\n{sequence}\n")

def concatenate_proteins(base_dir, species):
    input_dir = os.path.join(base_dir, species, "CDS_protein")
    output_dir = os.path.join(base_dir, species, "CDS_protein")
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{species}_concatenated_proteins.fasta")

    concatenated_proteins = {}

    for gene in gene_order:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id
                for gene_name in gene_order:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break
                if accession not in concatenated_proteins:
                    concatenated_proteins[accession] = ""
                concatenated_proteins[accession] += str(record.seq)

    with open(output_file, "w") as output:
        for accession, sequence in concatenated_proteins.items():
            output.write(f">{accession}\n{sequence}\n")

base_directory = "sequences/CDS"
species_list = [species for species in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, species))]

for species in species_list:
    concatenate_sequences(base_directory, species)
    concatenate_proteins(base_directory, species)

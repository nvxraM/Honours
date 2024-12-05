import os
from Bio import SeqIO

gene_order = ["ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6"]

def concatenate_sequences(base_dir, species):
    input_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    nucleotide_output_file = os.path.join(input_dir, f"{species}_concatenated_sequences_nucleotide.fasta")

    concatenated_nucleotide_sequences = {}

    # Process nucleotide sequences
    for gene in gene_order:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id
                if accession not in concatenated_nucleotide_sequences:
                    concatenated_nucleotide_sequences[accession] = ""
                concatenated_nucleotide_sequences[accession] += str(record.seq)

    with open(nucleotide_output_file, "w") as output:
        for accession, sequence in concatenated_nucleotide_sequences.items():
            output.write(f">{accession}\n{sequence}\n")

base_directory = "sequences/CDS"
species_list = ["Phocoena"]

for species in species_list:
    concatenate_sequences(base_directory, species)

import os
from Bio import SeqIO

gene_order = ["ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6"]

def concatenate_sequences(base_dir, species):
    aligned_input_dir = os.path.join(base_dir, species, "CDS_nucleotide", "aligned")
    os.makedirs(aligned_input_dir, exist_ok=True)

    nucleotide_output_file = os.path.join(aligned_input_dir, f"{species}_concatenated_sequences_nucleotide_aligned.fasta")

    concatenated_nucleotide_sequences = {}

    # Process aligned nucleotide sequences
    for gene in gene_order:
        gene_file = os.path.join(aligned_input_dir, f"{gene}.afa")
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

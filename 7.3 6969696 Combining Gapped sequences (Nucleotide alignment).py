import os
from Bio import SeqIO

gene_order = ["ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","CYTB","ND6"]

def concatenate_sequences(base_dir, species):
    nucleotide_input_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    gapped_input_dir = os.path.join(base_dir, species, "CDS_nucleotide_gapped")
    nucleotide_output_dir = os.path.join(base_dir, species, "CDS_nucleotide")
    gapped_output_dir = os.path.join(base_dir, species, "CDS_nucleotide_gapped")
    os.makedirs(nucleotide_output_dir, exist_ok=True)
    os.makedirs(gapped_output_dir, exist_ok=True)

    nucleotide_output_file = os.path.join(nucleotide_output_dir, f"{species}_concatenated_sequences_nucleotide.fasta")
    gapped_output_file = os.path.join(gapped_output_dir, f"{species}_concatenated_sequences_gapped.fasta")

    concatenated_nucleotide_sequences = {}
    concatenated_gapped_sequences = {}

    # Process nucleotide sequences
    for gene in gene_order:
        gene_file = os.path.join(nucleotide_input_dir, f"{gene}.fasta")
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

    # Process gapped sequences
    for gene in gene_order:
        gene_file = os.path.join(gapped_input_dir, f"{gene}.fasta")
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
                if accession not in concatenated_gapped_sequences:
                    concatenated_gapped_sequences[accession] = ""
                concatenated_gapped_sequences[accession] += str(record.seq)

    with open(gapped_output_file, "w") as output:
        for accession, sequence in concatenated_gapped_sequences.items():
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
species_list = ["Phocoena"]

for species in species_list:
    concatenate_sequences(base_directory, species)
    concatenate_proteins(base_directory, species)

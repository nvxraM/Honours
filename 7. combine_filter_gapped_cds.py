import os
from Bio import SeqIO

# Define the original gene order
GENE_ORDER = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"
]

# Set this to True if you wish to exclude ND6 from the concatenation
EXCLUDE_ND6 = True

def get_filtered_gene_order():
    """
    Return the gene order, excluding ND6 if EXCLUDE_ND6 is True.
    Otherwise, return the full gene order.
    """
    if EXCLUDE_ND6:
        return [gene for gene in GENE_ORDER if gene != "ND6"]
    else:
        return GENE_ORDER

def process_gene_sequences(gene_sequences):
    """
    Given a dictionary of gene sequences {accession: seq_string}, apply the following steps:
    1. Remove trailing gaps ('-') from all sequences.
    2. Determine the shortest sequence length across all sequences and truncate all to that length.
    3. Ensure the final length is divisible by 3. If not, trim off extra nucleotides at the end.
    4. Identify any codon positions (every 3 nucleotides) that contain 'N' in any sequence. Remove
       those codons from all sequences.
    5. After removing codons with 'N', ensure all sequences remain the same length and are aligned.

    Returns:
        dict: {accession: processed_sequence_str}
    """
    if not gene_sequences:
        return gene_sequences

    # Step 1: Remove trailing gaps
    # Using rstrip('-') to remove any trailing '-' characters
    for acc in gene_sequences:
        gene_sequences[acc] = gene_sequences[acc].rstrip('-')

    # Step 2: Truncate all sequences to the shortest length
    min_length = min(len(seq) for seq in gene_sequences.values())
    for acc in gene_sequences:
        gene_sequences[acc] = gene_sequences[acc][:min_length]

    # Step 3: Ensure length is divisible by 3
    remainder = min_length % 3
    if remainder != 0:
        min_length -= remainder
        for acc in gene_sequences:
            gene_sequences[acc] = gene_sequences[acc][:min_length]

    # If no sequence length remains, just return empty sequences
    if min_length == 0:
        return {acc: "" for acc in gene_sequences}

    # Step 4: Remove codons containing 'N'
    # We have min_length which is divisible by 3 now
    num_codons = min_length // 3

    # Create a mask that indicates which codon positions are retained
    codon_mask = [True] * num_codons

    # Check each codon position in all sequences
    for i in range(num_codons):
        start = i * 3
        end = start + 3
        # If any sequence at this codon position has 'N', mark codon as False
        for seq in gene_sequences.values():
            if 'N' in seq[start:end]:
                codon_mask[i] = False
                break

    # Build new sequences excluding codons with 'N'
    for acc in gene_sequences:
        seq = gene_sequences[acc]
        new_seq_parts = [seq[i*3:(i*3)+3] for i in range(num_codons) if codon_mask[i]]
        gene_sequences[acc] = "".join(new_seq_parts)

    # After removing codons, ensure all sequences align at the same length again
    if gene_sequences:
        final_min_length = min(len(s) for s in gene_sequences.values())
        for acc in gene_sequences:
            gene_sequences[acc] = gene_sequences[acc][:final_min_length]

    return gene_sequences

def concatenate_gapped_sequences(base_dir, genus):
    """
    For a given genus directory, this function:
    1. Identifies the FASTA files in CDS_nucleotide_gapped corresponding to the genes.
    2. Processes each gene's sequences by applying the rules described in process_gene_sequences().
    3. Concatenates all processed gene sequences for each accession in the order of the filtered gene set.
    4. Writes the concatenated sequences to a single FASTA file.

    Args:
        base_dir (str): Base directory containing genus subdirectories.
        genus (str): Name of the genus subdirectory (e.g., 'Balaena').
    """
    # Input and output directories for this genus
    input_dir = os.path.join(base_dir, genus, "CDS_nucleotide_gapped")
    output_dir = os.path.join(base_dir, genus, "CDS_nucleotide_gapped")
    os.makedirs(output_dir, exist_ok=True)

    # Output file path for concatenated sequences
    output_file = os.path.join(output_dir, f"{genus}_concatenated_gapped_sequences.fasta")

    # Dictionary to hold concatenated sequences {accession: combined_sequence}
    concatenated_sequences = {}

    # Get the gene order (excluding ND6 if EXCLUDE_ND6 is True)
    filtered_genes = get_filtered_gene_order()

    # Process each gene
    for gene in filtered_genes:
        gene_file = os.path.join(input_dir, f"{gene}.fasta")

        # If the gene file does not exist, skip it
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping.")
            continue

        # Parse the FASTA file for the current gene
        gene_sequences = {}
        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id
                # Remove any appended gene name from the accession
                # Example: ABC123_ND2 -> ABC123
                for gene_name in GENE_ORDER:
                    if f"_{gene_name}" in accession:
                        accession = accession.replace(f"_{gene_name}", "")
                        break
                gene_sequences[accession] = str(record.seq)

        # Process the gene sequences as per the rules
        processed_sequences = process_gene_sequences(gene_sequences)

        # Concatenate the processed sequences to our main dictionary
        for acc, seq in processed_sequences.items():
            if acc not in concatenated_sequences:
                concatenated_sequences[acc] = ""
            concatenated_sequences[acc] += seq

    # Write the final concatenated sequences to the output file
    with open(output_file, "w") as output:
        for accession, sequence in concatenated_sequences.items():
            output.write(f">{accession}\n{sequence}\n")

# Base directory containing the genus subdirectories
base_directory = "sequences/CDS"

# Identify all directories in the base directory that represent genera
genus_list = [
    g for g in os.listdir(base_directory)
    if os.path.isdir(os.path.join(base_directory, g))
]

# Process each genus directory
for genus in genus_list:
    concatenate_gapped_sequences(base_directory, genus)

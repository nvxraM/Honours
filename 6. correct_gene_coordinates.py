import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def add_gaps_to_all_species(base_dir):
    """
    For each species directory within the specified base directory:
    1. Locate aligned protein sequences (.afa files).
    2. Identify the corresponding nucleotide sequences for the same genes.
    3. Insert codon-based gaps ('---') into the nucleotide sequences wherever the aligned protein sequence has gaps ('-').
    4. Save the resulting gapped nucleotide sequences into a separate output directory.

    Parameters
    ----------
    base_dir : str
        The path to the base directory containing species subdirectories. Each species directory should contain:
        - A 'CDS_protein/aligned' subdirectory with aligned protein sequences (.afa).
        - A 'CDS_nucleotide' subdirectory with the corresponding gene sequences (.fasta).
    """

    # Check if the base directory exists
    if not os.path.exists(base_dir):
        print(f"Base directory {base_dir} does not exist.")
        return

    # Identify all species directories within the base directory
    species_dirs = [
        d for d in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, d))
    ]

    # Process each species directory
    for species in species_dirs:
        # Construct path to the aligned protein directory for the current species
        aligned_protein_dir = os.path.join(base_dir, species, 'CDS_protein', 'aligned')

        # Skip species if aligned protein directory doesn't exist
        if not os.path.exists(aligned_protein_dir):
            print(f"Aligned protein directory {aligned_protein_dir} does not exist for species {species}")
            continue

        # Iterate over all files in the aligned protein directory
        for filename in os.listdir(aligned_protein_dir):
            # Only process alignment files with '.afa' extension
            if not filename.endswith('.afa'):
                continue

            # Full path to the aligned protein file
            protein_path = os.path.join(aligned_protein_dir, filename)

            # Extract the gene name from the file name.
            # Assumes filenames are structured as "GENENAME.afa"
            gene_name = os.path.splitext(filename)[0]

            # Path to the corresponding nucleotide gene file
            gene_path = os.path.join(base_dir, species, 'CDS_nucleotide', f"{gene_name}.fasta")

            # Check if the corresponding gene file exists
            if not os.path.exists(gene_path):
                print(f"Gene file {gene_path} not found for protein alignment {filename} in species {species}")
                continue

            # Parse the aligned protein sequences
            protein_records = list(SeqIO.parse(protein_path, 'fasta'))
            # Parse the corresponding gene sequences
            gene_records = list(SeqIO.parse(gene_path, 'fasta'))

            # Create a dictionary to quickly access gene records by their ID.
            # The keys match the record IDs in the protein alignments.
            gene_dict = {record.id: record for record in gene_records}

            modified_gene_records = []

            # For each protein record, insert gaps into the corresponding gene sequence
            for protein_record in protein_records:
                # Extract accession number or unique identifier
                # Assumes ID structure is something like "ACCNUM_xxx".
                accession_number = protein_record.id.split('_')[0]

                # Retrieve the corresponding gene sequence record by the full ID
                gene_record = gene_dict.get(protein_record.id)
                if not gene_record:
                    print(f"No corresponding gene found for {protein_record.id} in species {species}")
                    continue

                # We'll build a new gene sequence with gaps inserted
                new_gene_seq = []
                gene_index = 0  # Tracks the position within the nucleotide sequence

                # Iterate over each amino acid in the protein alignment
                for aa in protein_record.seq:
                    if aa == '-':
                        # Protein gap -> insert a codon-sized gap of '---' in the nucleotide sequence
                        new_gene_seq.append('---')
                    else:
                        # For a valid amino acid, extract the next codon (3 nucleotides) from the gene
                        if gene_index < len(gene_record.seq):
                            codon = str(gene_record.seq[gene_index:gene_index + 3])
                            new_gene_seq.append(codon)
                            gene_index += 3
                        else:
                            # If we run out of gene sequence, this indicates a mismatch or incomplete data
                            # In a perfect scenario this shouldn't happen, but we don't halt processing
                            print(
                                f"Ran out of nucleotides in {gene_record.id} for species {species} at amino acid {aa}.")
                            break

                # Combine the new sequence into a single string
                new_gene_seq_str = ''.join(new_gene_seq)

                # Create a new SeqRecord for the modified gene sequence with gaps
                new_gene_record = SeqRecord(
                    Seq(new_gene_seq_str),
                    id=gene_record.id,
                    description=gene_record.description
                )
                modified_gene_records.append(new_gene_record)

            # Create the output directory for gapped nucleotide sequences if it doesn't already exist
            species_output_dir = os.path.join(base_dir, species, 'CDS_nucleotide_gapped')
            os.makedirs(species_output_dir, exist_ok=True)

            # Write the modified (gapped) gene sequences to a FASTA file
            output_path = os.path.join(species_output_dir, f"{gene_name}.fasta")
            SeqIO.write(modified_gene_records, output_path, 'fasta')
            print(f"Gapped gene sequences for {gene_name} in species {species} saved to {output_path}")

# Define the main directory containing species-specific subdirectories
main_dir = r'sequences/CDS'

# Run the function to process all species
add_gaps_to_all_species(main_dir)

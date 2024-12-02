import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to add gaps to gene sequences based on aligned protein sequences
def add_gaps_to_all_species(base_dir):
    # Check if the base directory exists before proceeding
    if not os.path.exists(base_dir):
        print(f"Base directory {base_dir} does not exist.")
        return

    # Iterate over each species folder to find corresponding nucleotide and protein sequences
    species_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    for species in species_dirs:
        aligned_protein_dir = os.path.join(base_dir, species, 'CDS_protein', 'aligned')
        if not os.path.exists(aligned_protein_dir):
            print(f"Aligned protein directory {aligned_protein_dir} does not exist for species {species}")
            continue

        # Iterate over all files in the aligned protein directory
        for filename in os.listdir(aligned_protein_dir):
            if filename.endswith('.afa'):  # Process only alignment files with '.afa' extension
                protein_path = os.path.join(aligned_protein_dir, filename)

                # Extract gene name from filename (e.g., ATP6, COX1, etc.)
                gene_name = filename.split('_')[0]

                gene_path = os.path.join(base_dir, species, 'CDS_nucleotide', f"{gene_name}.fasta")

                # Check if the corresponding gene file exists before proceeding
                if not os.path.exists(gene_path):
                    print(f"Gene file {gene_path} not found for protein alignment {filename} in species {species}")
                    continue

                # Parse the aligned protein sequences and corresponding gene sequences using iterators
                protein_records = SeqIO.parse(protein_path, 'fasta')
                gene_records = SeqIO.parse(gene_path, 'fasta')

                # Create a dictionary to quickly access gene records by their ID
                gene_dict = {record.id: record for record in gene_records}

                modified_gene_records = []

                # Iterate over each protein record to modify the corresponding gene sequence
                for protein_record in protein_records:
                    accession_number = protein_record.id.split('_')[0]  # Extract accession number
                    gene_record = gene_dict.get(protein_record.id)
                    if not gene_record:
                        # Skip if no corresponding gene sequence is found for the protein
                        print(f"No corresponding gene found for {protein_record.id} in species {species}")
                        continue

                    new_gene_seq = []  # List to store the modified gene sequence with gaps
                    gene_index = 0  # Index to track position in the gene sequence

                    # Iterate over each amino acid in the protein sequence
                    for aa in protein_record.seq:
                        if aa == '-':  # Add a gap ('---') for each gap in the protein sequence
                            new_gene_seq.append('---')
                        else:
                            # Add the corresponding codon (3 nucleotides) from the gene sequence
                            if gene_index < len(gene_record.seq):
                                new_gene_seq.append(str(gene_record.seq[gene_index:gene_index+3]))
                                gene_index += 3

                    # Combine the modified gene sequence into a single string
                    new_gene_seq = ''.join(new_gene_seq)
                    # Create a new SeqRecord for the modified gene sequence
                    new_gene_record = SeqRecord(Seq(new_gene_seq), id=gene_record.id, description=gene_record.description)
                    modified_gene_records.append(new_gene_record)

                # Determine the output directory for the specific species
                species_output_dir = os.path.join(base_dir, species, 'CDS_nucleotide_gapped')
                os.makedirs(species_output_dir, exist_ok=True)

                # Write the modified gene sequences to the output file
                output_path = os.path.join(species_output_dir, f"{gene_name}.fasta")
                SeqIO.write(modified_gene_records, output_path, 'fasta')
                print(f"Gapped gene sequences for {gene_name} in species {species} saved to {output_path}")


# Define the main directory containing species-specific subdirectories
main_dir = r'sequences/CDS'

# Run the function to process all species
add_gaps_to_all_species(main_dir)

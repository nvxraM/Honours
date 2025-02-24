import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def process_genbank_files(base_directory):
    """
    Process GenBank files from a base directory.

    This function does the following:
    - Iterates over a directory structure containing GenBank files.
    - Attempts to find corresponding FASTA files with matching IDs.
    - Extracts CDS_Genus (coding sequence) information (both nucleotide and protein).
    - Logs progress and any issues encountered.
    - Adjusts sequences to remove stop codons and partial stop codons so that the sequences are fully divisible by 3.
    """
    format_type = 'gb'
    input_directory = os.path.join(base_directory, format_type)
    output_directory = os.path.join(base_directory, "CDS_Genus")
    log_file_path = os.path.join(base_directory, "CDS_Acquisition_Log.txt")

    # Check if the input directory exists
    if not os.path.exists(input_directory):
        log_message(log_file_path, f"Input directory {input_directory} does not exist.")
        return

    # Loop through each species folder
    for species_folder in os.listdir(input_directory):
        species_path = os.path.join(input_directory, species_folder)

        if os.path.isdir(species_path):
            # Extract genus name from the folder name (assumes structure like Genus_species)
            genus_name = species_folder.split('_')[0]

            # --------------------------------------------------------------
            # Removed the merging logic here. Originally we had:
            # if genus_name == "Phocoenoides":
            #     genus_name = "Phocoena"
            # --------------------------------------------------------------

            # Create output directories for this genus
            genus_cds_path = os.path.join(output_directory, genus_name)
            gene_dir = os.path.join(genus_cds_path, 'CDS_nucleotide')
            protein_dir = os.path.join(genus_cds_path, 'CDS_protein')
            os.makedirs(gene_dir, exist_ok=True)
            os.makedirs(protein_dir, exist_ok=True)

            # Find the corresponding FASTA file for this species
            fasta_file_path = find_fasta_file(base_directory, species_folder)
            if not fasta_file_path:
                log_message(log_file_path, f"FASTA file for species {species_folder} not found.")
                continue

            # Parse the FASTA file to get a dictionary of SeqRecords by accession ID
            with open(fasta_file_path, 'r') as fasta_handle:
                fasta_records = {record.id: record for record in SeqIO.parse(fasta_handle, "fasta")}

            # Now process all GenBank files in the species path
            for genbank_file in os.listdir(species_path):
                if genbank_file.endswith('.gb'):
                    file_path = os.path.join(species_path, genbank_file)
                    with open(file_path, 'r') as gb_handle:
                        for record in SeqIO.parse(gb_handle, format_type):
                            # Extract and save CDS_Genus sequences to files
                            extract_and_save(
                                [record],
                                gene_dir,
                                protein_dir,
                                fasta_records,
                                species_folder,
                                log_file_path
                            )

            # Log that we processed this species
            log_message(log_file_path, f"Processed CDS_Genus for species: {species_folder}")

    # After processing all species, recheck all sequences for proper length
    recheck_sequence_lengths(output_directory, log_file_path)


def find_fasta_file(base_directory, species_folder):
    """
    Attempts to find a FASTA file corresponding to a given species_folder.
    It first checks a "fasta/<species_folder>/" directory.
    If not found, it searches recursively in the "fasta" base directory.
    """
    species_fasta_dir = os.path.join(base_directory, "fasta", species_folder)
    if os.path.exists(species_fasta_dir):
        # Look for a .fasta file in species-specific directory
        for file in os.listdir(species_fasta_dir):
            if file.endswith('.fasta'):
                return os.path.join(species_fasta_dir, file)

    # If not found, search recursively in the main "fasta" directory
    fasta_base_dir = os.path.join(base_directory, "fasta")
    for root, _, files in os.walk(fasta_base_dir):
        for file in files:
            if file.endswith('.fasta'):
                fasta_path = os.path.join(root, file)
                # Check if the species_folder appears in the FASTA record IDs
                with open(fasta_path, 'r') as fasta_handle:
                    for record in SeqIO.parse(fasta_handle, "fasta"):
                        if species_folder in record.id:
                            return fasta_path
    return None


def extract_and_save(records, gene_directory, protein_directory, fasta_records, species_name, log_file_path):
    """
    Extracts CDS_Genus nucleotide and protein sequences from a list of SeqRecords.
    Writes them into appropriate FASTA files grouped by gene.
    """
    gene_sequences = {}
    protein_sequences = {}

    for record in records:
        # Identify CDS_Genus features (features of type "CDS_Genus" with a "translation" qualifier)
        cds_features = [
            feature for feature in record.features
            if feature.type == "CDS" and "translation" in feature.qualifiers
        ]

        # Process each CDS_Genus feature
        for feature in cds_features:
            # Determine gene name from feature qualifiers
            if "gene" in feature.qualifiers:
                gene_name = normalize_gene_name(feature.qualifiers['gene'][0])
            elif "product" in feature.qualifiers:
                # If no gene qualifier, try to deduce from product (specifically for NADH genes)
                product_name = feature.qualifiers['product'][0].lower()
                if product_name.startswith("nadh dehydrogenase subunit"):
                    subunit_number = product_name.split()[-1]
                    if subunit_number.isdigit():
                        gene_name = f"ND{subunit_number}"
                    else:
                        # If unable to parse subunit number properly, skip
                        continue
                else:
                    # Not a known NADH product, skip
                    continue
            else:
                # No identifiable gene name or recognizable product
                continue

            # Extract the appropriate sequence (check if record has sequence or use FASTA)
            accession_id = record.id.split()[0]
            location = feature.location

            if len(record.seq) == 0:
                # If the GenBank record doesn't have sequence, try to get from FASTA
                if accession_id in fasta_records:
                    seq_to_use = fasta_records[accession_id].seq
                    sequence_source = 'FASTA record'
                else:
                    log_message(log_file_path, f"Sequence not found for accession {accession_id}, gene {gene_name}.")
                    continue
            else:
                seq_to_use = record.seq
                sequence_source = 'GenBank record'

            # Extract the nucleotide sequence for this CDS_Genus feature
            nucleotide_sequence = str(feature.extract(seq_to_use))
            log_message(log_file_path,
                        f"Extracted sequence for accession {accession_id}, gene {gene_name} from {sequence_source}. Strand: {location.strand}")

            # Remove stop codon or partial stop codons so that length is divisible by 3
            nucleotide_sequence = remove_stop_codon(nucleotide_sequence, accession_id,
                                                    record.annotations.get('organism', 'Unknown'), log_file_path)

            # Extract the protein sequence from the feature translation qualifier
            protein_sequence = feature.qualifiers['translation'][0]

            # Create a header for the sequence
            species_gene_header = f">{record.id}_{species_name}"

            # Store sequences by gene name (may be multiple sequences per gene)
            if gene_name not in gene_sequences:
                gene_sequences[gene_name] = []
            gene_sequences[gene_name].append(f"{species_gene_header}\n{nucleotide_sequence}\n")

            if gene_name not in protein_sequences:
                protein_sequences[gene_name] = []
            protein_sequences[gene_name].append(f"{species_gene_header}\n{protein_sequence}\n")

    # Write out collected gene and protein sequences
    for gene_name, sequences in gene_sequences.items():
        combined_gene_name = normalize_combined_gene_name(gene_name)
        gene_file = os.path.join(gene_directory, f"{combined_gene_name}.fasta")
        with open(gene_file, 'a') as gf:
            gf.writelines(sequences)

    for gene_name, sequences in protein_sequences.items():
        combined_gene_name = normalize_combined_gene_name(gene_name)
        protein_file = os.path.join(protein_directory, f"{combined_gene_name}.fasta")
        with open(protein_file, 'a') as pf:
            pf.writelines(sequences)


def remove_stop_codon(sequence, accession_id, species_name, log_file_path):
    """
    Removes full or partial stop codons at the end of a sequence so that the final sequence length
    is divisible by 3. Logs when stop codons are found and removed.
    """
    stop_codons = ["TAA", "TAG", "AGA", "AGG"]

    # Check for a full stop codon (3 nucleotides) at the end
    if len(sequence) >= 3 and sequence[-3:] in stop_codons:
        log_message(log_file_path,
                    f"{species_name} {accession_id} [{sequence[-3:]}] Full stop codon found and removed.")
        sequence = sequence[:-3]
        if len(sequence) % 3 == 0:
            return sequence

    # If not divisible by 3, look for partial stop codons (less than 3 nucleotides at the end)
    while len(sequence) > 0 and len(sequence) % 3 != 0:
        if len(sequence) >= 2 and sequence[-2:] in ["TA", "TG", "AG"] and (len(sequence) - 2) % 3 == 0:
            log_message(log_file_path,
                        f"{species_name} {accession_id} [{sequence[-2:]}] Partial stop codon (2 nt) found and removed.")
            sequence = sequence[:-2]
        elif len(sequence) >= 1 and sequence[-1:] in ["T", "A", "G"] and (len(sequence) - 1) % 3 == 0:
            log_message(log_file_path,
                        f"{species_name} {accession_id} [{sequence[-1:]}] Partial stop codon (1 nt) found and removed.")
            sequence = sequence[:-1]
        else:
            # If removing partial codons doesn't yield a divisible-by-3 sequence, break
            break
    return sequence


def normalize_gene_name(gene_name):
    """
    Normalizes gene names to a standard format (e.g., NADH genes to ND#).
    """
    gene_name = gene_name.upper()
    if gene_name == "NADH1":
        return "ND1"
    elif gene_name == "NADH2":
        return "ND2"
    elif gene_name == "NADH3":
        return "ND3"
    elif gene_name == "NADH4":
        return "ND4"
    elif gene_name == "NADH4L":
        return "ND4L"
    elif gene_name == "NADH5":
        return "ND5"
    elif gene_name == "NADH6":
        return "ND6"
    elif gene_name in ["ATPASE 8", "ATPASE8"]:
        return "ATP8"
    elif gene_name in ["ATPASE 6", "ATPASE6", "ATP"]:
        return "ATP6"
    elif gene_name in ["COI", "COXI", "CO1"]:
        return "COX1"
    elif gene_name in ["COII", "COXII", "CO2"]:
        return "COX2"
    elif gene_name in ["COIII", "COXIII", "CO3"]:
        return "COX3"
    else:
        return gene_name


def normalize_combined_gene_name(gene_name):
    """
    Ensures consistency in gene naming across different gene variants.
    """
    if gene_name in ["ND1", "NAD1"]:
        return "ND1"
    elif gene_name in ["ND2", "NAD2"]:
        return "ND2"
    elif gene_name in ["ND3", "NAD3"]:
        return "ND3"
    elif gene_name in ["ND4", "NAD4"]:
        return "ND4"
    elif gene_name in ["ND4L", "NAD4L"]:
        return "ND4L"
    elif gene_name in ["ND5", "NAD5"]:
        return "ND5"
    elif gene_name in ["ND6", "NAD6"]:
        return "ND6"
    elif gene_name in ["COX1", "CO1"]:
        return "COX1"
    elif gene_name in ["COX2", "CO2"]:
        return "COX2"
    elif gene_name in ["COX3", "CO3"]:
        return "COX3"
    elif gene_name in ["ATP6", "ATPASE6"]:
        return "ATP6"
    elif gene_name in ["ATP8", "ATPASE8"]:
        return "ATP8"
    elif gene_name in ["COB"]:
        return "CYTB"
    else:
        return gene_name


def recheck_sequence_lengths(directory, log_file_path):
    """
    After all sequences have been initially processed and saved,
    re-check them to ensure they are divisible by 3 and remove any
    leftover partial stop codons if needed.
    """
    for genus_folder in os.listdir(directory):
        genus_path = os.path.join(directory, genus_folder, 'CDS_nucleotide')
        if os.path.isdir(genus_path):
            for gene_file in os.listdir(genus_path):
                if gene_file.endswith('.fasta'):
                    file_path = os.path.join(genus_path, gene_file)
                    adjusted_sequences = []
                    # Parse each sequence and remove stop codons if necessary
                    with open(file_path, 'r') as f:
                        for record in SeqIO.parse(f, "fasta"):
                            adjusted_sequence = remove_stop_codon(str(record.seq), record.id, genus_folder,
                                                                  log_file_path)
                            adjusted_sequences.append((record.id, adjusted_sequence))
                    # Rewrite the file with adjusted sequences
                    with open(file_path, 'w') as f:
                        for record_id, sequence in adjusted_sequences:
                            f.write(f">{record_id}\n{sequence}\n")


def log_message(log_file_path, message):
    """
    Logs a message to the specified log file.
    """
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + "\n")


if __name__ == "__main__":
    # Specify the base directory containing 'gb' and 'fasta' subdirectories
    base_directory = "sequences/"
    process_genbank_files(base_directory)

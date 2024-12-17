import os
import shutil
import re

# Directory paths
fasta_dir = 'sequences/fasta'  # Directory containing FASTA files organized by species
gb_dir = 'sequences/gb'  # Directory containing GenBank (GB) files organized by species
excluded_fasta_dir = 'sequences/excluded/fasta'  # Directory to move excluded FASTA files/species
excluded_gb_dir = 'sequences/excluded/gb'  # Directory to move excluded GB files/species

# Ensure the excluded directories exist
os.makedirs(excluded_fasta_dir, exist_ok=True)
os.makedirs(excluded_gb_dir, exist_ok=True)


def count_fasta_headers(fasta_file):
    """
    Count the number of headers (lines starting with '>') in a given FASTA file.

    Parameters:
        fasta_file (str): Path to the FASTA file.

    Returns:
        int: The number of FASTA headers in the file.
    """
    with open(fasta_file, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))


def check_cds_count(gb_file_path):
    """
    Count the number of CDS_Genus features in a given GenBank file.

    Parameters:
        gb_file_path (str): Path to the GenBank file.

    Returns:
        int: The number of CDS_Genus entries found in the GenBank file.
    """
    with open(gb_file_path, 'r') as file:
        content = file.read()
        # Matches lines indicating a CDS_Genus entry. CDS_Genus often start with a line containing "CDS_Genus" preceded by whitespace.
        cds_count = len(re.findall(r'\n\s+CDS_Genus\s', content))
        return cds_count


def exclude_species_based_on_fasta(fasta_directory, gb_directory, header_threshold=4):
    """
    Exclude species whose FASTA files have a number of headers less than or equal to a given threshold.
    Moves both FASTA and corresponding GB directories to the excluded directories.

    Parameters:
        fasta_directory (str): Directory containing species subdirectories with FASTA files.
        gb_directory (str): Directory containing species subdirectories with GB files.
        header_threshold (int): The minimum number of headers required to keep the species.
                                Species with headers <= this value are moved to excluded directories.
    """
    with os.scandir(fasta_directory) as entries:
        for entry in entries:
            if entry.is_dir():
                # Each subdirectory under 'fasta_directory' should correspond to a species
                species = entry.name
                species_fasta_path = os.path.join(entry.path, f"{species}.fasta")

                if os.path.isfile(species_fasta_path):
                    # Count how many headers the FASTA file has
                    header_count = count_fasta_headers(species_fasta_path)

                    # If headers are below or equal to threshold, exclude this species
                    if header_count <= header_threshold:
                        excluded_fasta_path = os.path.join(excluded_fasta_dir, species)

                        # Remove existing species folder in excluded if present, then move
                        if os.path.exists(excluded_fasta_path):
                            shutil.rmtree(excluded_fasta_path)
                        shutil.move(entry.path, excluded_fasta_path)
                        print(f"Moved {species} FASTA to excluded due to only {header_count} headers.")

                        # Also move the corresponding GB folder
                        species_gb_path = os.path.join(gb_directory, species)
                        if os.path.exists(species_gb_path):
                            excluded_gb_path = os.path.join(excluded_gb_dir, species)
                            if os.path.exists(excluded_gb_path):
                                shutil.rmtree(excluded_gb_path)
                            shutil.move(species_gb_path, excluded_gb_path)
                            print(f"Moved {species} GB to excluded due to corresponding FASTA headers.")


def check_and_move_species_based_on_cds(fasta_directory, gb_directory):
    """
    For each species in the GB directory:
    - Checks each GenBank file for the number of CDS_Genus features.
    - If a file does not have the expected number of CDS_Genus features (e.g., 13), the file is moved to the excluded directory.
    - The corresponding entries in the FASTA file are also moved to the excluded directory.
    - The original FASTA file is updated to remove the excluded sequences.

    Parameters:
        fasta_directory (str): Directory containing species subdirectories with FASTA files.
        gb_directory (str): Directory containing species subdirectories with GB files.
    """
    with os.scandir(gb_directory) as entries:
        for entry in entries:
            if entry.is_dir():
                species = entry.name
                species_gb_path = entry.path
                gb_files = [f for f in os.listdir(species_gb_path) if f.endswith('.gb')]

                # Keep track of accessions to exclude
                excluded_accessions = []

                # Check each GB file for the correct number of CDS_Genus features
                for gb_file in gb_files:
                    gb_file_path = os.path.join(species_gb_path, gb_file)
                    accession = gb_file.split('.')[0]

                    if os.path.isfile(gb_file_path):
                        cds_count = check_cds_count(gb_file_path)
                        # If the number of CDS_Genus entries is not as expected (13 in this example),
                        # move the file and track the accession.
                        if cds_count != 13:
                            excluded_gb_path = os.path.join(excluded_gb_dir, species)
                            os.makedirs(excluded_gb_path, exist_ok=True)

                            # Move the GB file to the excluded directory
                            shutil.move(gb_file_path, os.path.join(excluded_gb_path, gb_file))
                            print(f"Moved {gb_file} to excluded/gb due to {cds_count} CDS_Genus entries.")
                            excluded_accessions.append(accession)

                # If we excluded any GB files, we must also update the corresponding FASTA files
                if excluded_accessions:
                    species_fasta_path = os.path.join(fasta_directory, species, f"{species}.fasta")
                    if os.path.exists(species_fasta_path):
                        excluded_fasta_path = os.path.join(excluded_fasta_dir, species)
                        os.makedirs(excluded_fasta_path, exist_ok=True)
                        new_fasta_path = os.path.join(excluded_fasta_path, f"{species}.fasta")

                        # Write excluded sequences to a new FASTA file in the excluded directory
                        with open(species_fasta_path, 'r') as infile, open(new_fasta_path, 'w') as outfile:
                            write = False
                            for line in infile:
                                if line.startswith('>'):
                                    # Parse the accession number from header line
                                    acc = line.strip().split()[0][1:].split('.')[0]
                                    write = (acc in excluded_accessions)
                                if write:
                                    outfile.write(line)

                        # Create a temporary file to rewrite the original FASTA file without excluded sequences
                        temp_fasta_path = species_fasta_path + ".temp"
                        with open(species_fasta_path, 'r') as infile, open(temp_fasta_path, 'w') as outfile:
                            write = True
                            for line in infile:
                                if line.startswith('>'):
                                    acc = line.strip().split()[0][1:].split('.')[0]
                                    write = (acc not in excluded_accessions)
                                if write:
                                    outfile.write(line)

                        # Replace the original FASTA file with the updated one
                        os.replace(temp_fasta_path, species_fasta_path)
                        print(
                            f"Updated {species}.fasta to exclude sequences with accessions {', '.join(excluded_accessions)}.")


def run_additional_header_check(fasta_directory, gb_directory, header_threshold=4):
    """
    Run an additional check for species based on the FASTA header count and exclude them if they don't meet the threshold.

    Parameters:
        fasta_directory (str): Directory containing FASTA files.
        gb_directory (str): Directory containing GB files.
        header_threshold (int): Threshold for minimum number of headers.
    """
    exclude_species_based_on_fasta(fasta_directory, gb_directory, header_threshold)


# Run the exclusion processes
exclude_species_based_on_fasta(fasta_dir, gb_dir)
check_and_move_species_based_on_cds(fasta_dir, gb_dir)
run_additional_header_check(fasta_dir, gb_dir)

import os
import shutil
import re

# Paths to the fasta and GB folders
fasta_dir = 'sequences/fasta'  # Directory containing FASTA files
gb_dir = 'sequences/gb'  # Directory containing GB files
excluded_fasta_dir = 'sequences/excluded/fasta'  # Directory to store excluded FASTA files
excluded_gb_dir = 'sequences/excluded/gb'  # Directory to store excluded GB files

# Ensure the excluded directories exist
os.makedirs(excluded_fasta_dir, exist_ok=True)  # Create excluded FASTA directory if it doesn't exist
os.makedirs(excluded_gb_dir, exist_ok=True)  # Create excluded GB directory if it doesn't exist

def count_fasta_headers(fasta_file):
    with open(fasta_file, 'r') as file:
        return sum(1 for line in file if line.startswith('>'))

def check_cds_count(gb_file_path):
    with open(gb_file_path, 'r') as file:
        content = file.read()
        cds_count = len(re.findall(r'\n\s+CDS\s', content))
        return cds_count

def exclude_species_based_on_fasta(fasta_directory, gb_directory, header_threshold=4):
    with os.scandir(fasta_directory) as entries:
        for entry in entries:
            if entry.is_dir():
                species = entry.name
                species_fasta_path = os.path.join(entry.path, f"{species}.fasta")

                if os.path.isfile(species_fasta_path):
                    header_count = count_fasta_headers(species_fasta_path)
                    if header_count <= header_threshold:
                        excluded_fasta_path = os.path.join(excluded_fasta_dir, species)
                        if os.path.exists(excluded_fasta_path):
                            shutil.rmtree(excluded_fasta_path)
                        shutil.move(entry.path, excluded_fasta_path)
                        print(f"Moved {species} to excluded/fasta due to {header_count} headers.")

                        species_gb_path = os.path.join(gb_directory, species)
                        if os.path.exists(species_gb_path):
                            excluded_gb_path = os.path.join(excluded_gb_dir, species)
                            if os.path.exists(excluded_gb_path):
                                shutil.rmtree(excluded_gb_path)
                            shutil.move(species_gb_path, excluded_gb_path)
                            print(f"Moved {species} to excluded/gb due to corresponding FASTA file with {header_count} headers.")

def check_and_move_species_based_on_cds(fasta_directory, gb_directory):
    with os.scandir(gb_directory) as entries:
        for entry in entries:
            if entry.is_dir():
                species = entry.name
                species_gb_path = entry.path
                gb_files = [f for f in os.listdir(species_gb_path) if f.endswith('.gb')]

                excluded_accessions = []
                for gb_file in gb_files:
                    gb_file_path = os.path.join(species_gb_path, gb_file)
                    accession = gb_file.split('.')[0]

                    if os.path.isfile(gb_file_path):
                        cds_count = check_cds_count(gb_file_path)
                        if cds_count != 13:
                            excluded_gb_path = os.path.join(excluded_gb_dir, species)
                            if not os.path.exists(excluded_gb_path):
                                os.makedirs(excluded_gb_path)
                            shutil.move(gb_file_path, os.path.join(excluded_gb_path, gb_file))
                            print(f"Moved {gb_file} to excluded/gb due to {cds_count} CDS entries.")
                            excluded_accessions.append(accession)

                if excluded_accessions:
                    species_fasta_path = os.path.join(fasta_directory, species, f"{species}.fasta")
                    if os.path.exists(species_fasta_path):
                        excluded_fasta_path = os.path.join(excluded_fasta_dir, species)
                        if not os.path.exists(excluded_fasta_path):
                            os.makedirs(excluded_fasta_path)
                        new_fasta_path = os.path.join(excluded_fasta_path, f"{species}.fasta")

                        with open(species_fasta_path, 'r') as infile, open(new_fasta_path, 'w') as outfile:
                            write = False
                            for line in infile:
                                if line.startswith('>'):
                                    accession = line.strip().split()[0][1:].split('.')[0]
                                    if accession in excluded_accessions:
                                        write = True
                                    else:
                                        write = False
                                if write:
                                    outfile.write(line)

                        temp_fasta_path = species_fasta_path + ".temp"
                        with open(species_fasta_path, 'r') as infile, open(temp_fasta_path, 'w') as outfile:
                            write = True
                            for line in infile:
                                if line.startswith('>'):
                                    accession = line.strip().split()[0][1:].split('.')[0]
                                    if accession in excluded_accessions:
                                        write = False
                                    else:
                                        write = True
                                if write:
                                    outfile.write(line)
                        os.replace(temp_fasta_path, species_fasta_path)
                        print(f"Updated {species}.fasta to exclude sequences with accessions {', '.join(excluded_accessions)} due to no corresponding CDS entries.")

def run_additional_header_check(fasta_directory, gb_directory, header_threshold=4):
    exclude_species_based_on_fasta(fasta_directory, gb_directory, header_threshold)

exclude_species_based_on_fasta(fasta_dir, gb_dir)
check_and_move_species_based_on_cds(fasta_dir, gb_dir)
run_additional_header_check(fasta_dir, gb_dir)

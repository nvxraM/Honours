import os
import shutil
from Bio import SeqIO


def combine_subspecies_fasta(base_directory):
    """
    Combine FASTA files from subspecies folders into their corresponding base species folders.

    Process:
    1. For each folder inside `base_directory`, determine if it represents a base species or a subspecies.
    2. If it's identified as a subspecies folder (based on folder name):
       - Determine the corresponding base species name.
       - Create the base species folder if it doesn't exist.
       - Combine all .fasta files from the subspecies folder into a single .fasta file in the base species folder.
       - Move any other non-FASTA files into the base species folder.
       - Remove the now-empty subspecies folder.
    """
    print(f"Combining subspecies folders and FASTA files in {base_directory}")

    # Iterate over each species (or subspecies) folder in the base directory.
    for species_folder in os.listdir(base_directory):
        species_folder_path = os.path.join(base_directory, species_folder)

        # Check if the current item is a directory before processing.
        if os.path.isdir(species_folder_path):
            # Special case: "Delphinus_bairdii" should be treated as "Delphinus_delphis".
            if species_folder == "Delphinus_bairdii":
                base_species_name = "Delphinus_delphis"
            else:
                # Generally, the base species name is derived by taking the first two parts of the folder name.
                species_name_parts = species_folder.split('_')
                base_species_name = '_'.join(species_name_parts[:2])

            # If the folder name does not match the base species name, it is considered a subspecies folder.
            if base_species_name != species_folder:
                # Construct the base species folder path and ensure it exists.
                base_species_folder_path = os.path.join(base_directory, base_species_name)
                os.makedirs(base_species_folder_path, exist_ok=True)

                # Combine FASTA files from this subspecies folder into the corresponding base species folder.
                combine_fasta_files(species_folder_path, base_species_folder_path)

                # Move non-FASTA files (if any) to the base species folder.
                for file_name in os.listdir(species_folder_path):
                    source_file = os.path.join(species_folder_path, file_name)
                    target_file = os.path.join(base_species_folder_path, file_name)

                    # Move all non-FASTA files (e.g., metadata, txt files) to the base folder, if they don't already exist there.
                    if os.path.isfile(source_file) and not file_name.endswith('.fasta'):
                        if not os.path.exists(target_file):
                            shutil.move(source_file, target_file)
                            print(f"Moved {file_name} from {species_folder} to {base_species_name}")
                        else:
                            print(f"Skipping {file_name}, already exists in {base_species_name}")

                # Remove the now-empty subspecies folder.
                shutil.rmtree(species_folder_path)
                print(f"Removed folder {species_folder}")


def combine_fasta_files(subspecies_folder, base_species_folder):
    """
    Combine all FASTA files from a subspecies folder into a single FASTA file in the corresponding base species folder.

    Process:
    1. Identify all .fasta files in the subspecies folder.
    2. If a base species FASTA file already exists, append the subspecies sequences to it.
    3. If it does not exist, move the subspecies FASTA file to serve as the base species FASTA file.
    """
    # Construct the path for the combined base species FASTA file.
    base_species_fasta_path = os.path.join(base_species_folder, f"{os.path.basename(base_species_folder)}.fasta")

    # Loop through all FASTA files in the subspecies folder.
    for file_name in os.listdir(subspecies_folder):
        if file_name.endswith(".fasta"):
            subspecies_fasta_path = os.path.join(subspecies_folder, file_name)

            # Parse sequences from the subspecies FASTA file.
            subspecies_sequences = list(SeqIO.parse(subspecies_fasta_path, "fasta"))

            # If the base species FASTA already exists, append to it; otherwise, move the subspecies file as the starter.
            if os.path.exists(base_species_fasta_path):
                print(f"Combining FASTA files from {subspecies_fasta_path} into {base_species_fasta_path}")
                with open(base_species_fasta_path, "a") as base_fasta_file:
                    SeqIO.write(subspecies_sequences, base_fasta_file, "fasta")
            else:
                print(f"Creating new FASTA file {base_species_fasta_path} from {subspecies_fasta_path}")
                shutil.move(subspecies_fasta_path, base_species_fasta_path)

            print(
                f"Combined {len(subspecies_sequences)} sequences from {subspecies_fasta_path} into {base_species_fasta_path}")


def move_gb_files(base_directory):
    """
    Move .gb files from subspecies folders into their corresponding base species folder.

    Process:
    1. For each folder in `base_directory`, determine if it represents a base species or a subspecies.
    2. If subspecies, move all .gb files into the corresponding base species folder.
    3. Remove the empty subspecies folder afterward.
    """
    print(f"Moving GB files from subspecies to their base species folder in {base_directory}")

    # Iterate over all directories in the base directory.
    for species_folder in os.listdir(base_directory):
        species_folder_path = os.path.join(base_directory, species_folder)

        # Check if it is a directory.
        if os.path.isdir(species_folder_path):
            # Special case handling similar to the FASTA combination step.
            if species_folder == "Delphinus_bairdii":
                base_species_name = "Delphinus_delphis"
            else:
                species_name_parts = species_folder.split('_')
                base_species_name = '_'.join(species_name_parts[:2])

            # If it doesn't match the base species name, it's considered a subspecies folder.
            if base_species_name != species_folder:
                # Ensure the base species directory exists.
                base_species_folder_path = os.path.join(base_directory, base_species_name)
                os.makedirs(base_species_folder_path, exist_ok=True)

                # Move all .gb files from the subspecies folder to the base species folder.
                for file_name in os.listdir(species_folder_path):
                    if file_name.endswith(".gb"):
                        source_file = os.path.join(species_folder_path, file_name)
                        target_file = os.path.join(base_species_folder_path, file_name)

                        # Move the GB file if it doesn't exist in the target.
                        if not os.path.exists(target_file):
                            shutil.move(source_file, target_file)
                            print(f"Moved {file_name} from {species_folder} to {base_species_folder_path}")
                        else:
                            print(f"Skipping {file_name}, already exists in {base_species_folder_path}")

                # Remove the now-empty subspecies folder.
                shutil.rmtree(species_folder_path)
                print(f"Removed folder {species_folder}")


if __name__ == "__main__":
    # Directories where FASTA and GB sequences are stored.
    base_directory_fasta = "sequences/fasta"
    base_directory_gb = "sequences/gb"

    # First, combine subspecies FASTA files into their respective base species folders.
    combine_subspecies_fasta(base_directory_fasta)

    # Then, move GB files from subspecies folders to base species folders.
    move_gb_files(base_directory_gb)

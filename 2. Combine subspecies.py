import os
import shutil
from Bio import SeqIO


def combine_subspecies_fasta(base_directory):
    """
    Combines the FASTA files in the subspecies folder into the corresponding base species folder.
    """
    print(f"Combining subspecies folders and FASTA files in {base_directory}")

    # Loop through all the directories in the base folder
    for species_folder in os.listdir(base_directory):
        species_folder_path = os.path.join(base_directory, species_folder)

        # Check if it's a directory
        if os.path.isdir(species_folder_path):
            # Split the folder name into species and subspecies components
            if species_folder == "Delphinus_bairdii":
                base_species_name = "Delphinus_delphis"
            else:
                species_name_parts = species_folder.split('_')
                base_species_name = '_'.join(species_name_parts[:2])  # Take the first two parts as the base species name

            # If the folder is a subspecies folder (not the base species folder)
            if base_species_name != species_folder:
                # Determine the path to the base species folder
                base_species_folder_path = os.path.join(base_directory, base_species_name)
                # Create the base species folder if it doesn't exist
                os.makedirs(base_species_folder_path, exist_ok=True)

                # Combine FASTA files from the subspecies folder into the base species folder
                combine_fasta_files(species_folder_path, base_species_folder_path)

                # Move any other files from subspecies to base species folder
                for file_name in os.listdir(species_folder_path):
                    source_file = os.path.join(species_folder_path, file_name)
                    target_file = os.path.join(base_species_folder_path, file_name)

                    # Move non-FASTA files to the base species folder
                    if os.path.isfile(source_file) and not file_name.endswith('.fasta'):
                        if not os.path.exists(target_file):
                            shutil.move(source_file, target_file)
                            print(f"Moved {file_name} from {species_folder} to {base_species_name}")
                        else:
                            print(f"Skipping {file_name}, already exists in {base_species_name}")

                # Remove the subspecies folder after all files are moved and FASTA files are combined
                shutil.rmtree(species_folder_path)
                print(f"Removed folder {species_folder}")


def combine_fasta_files(subspecies_folder, base_species_folder):
    """
    Combines FASTA files from the subspecies folder into the corresponding FASTA file in the base species folder.
    """
    for file_name in os.listdir(subspecies_folder):
        if file_name.endswith(".fasta"):
            subspecies_fasta_path = os.path.join(subspecies_folder, file_name)
            base_species_fasta_path = os.path.join(base_species_folder,
                                                   f"{os.path.basename(base_species_folder)}.fasta")

            # Read sequences from subspecies FASTA file
            subspecies_sequences = list(SeqIO.parse(subspecies_fasta_path, "fasta"))

            # If base species FASTA file exists, append sequences
            if os.path.exists(base_species_fasta_path):
                print(f"Combining FASTA files from {subspecies_fasta_path} into {base_species_fasta_path}")
                with open(base_species_fasta_path, "a") as base_fasta_file:
                    # Append the sequences from the subspecies FASTA file to the base species FASTA file
                    SeqIO.write(subspecies_sequences, base_fasta_file, "fasta")
            else:
                # If base species FASTA does not exist, move the subspecies FASTA to the base species folder
                print(f"Creating new FASTA file {base_species_fasta_path} from {subspecies_fasta_path}")
                shutil.move(subspecies_fasta_path, base_species_fasta_path)
            print(
                f"Combined {len(subspecies_sequences)} sequences from {subspecies_fasta_path} into {base_species_fasta_path}")


def move_gb_files(base_directory):
    """
    Moves GB files from all subspecies folders into their respective base species folder (without combining contents).
    """
    print(f"Moving GB files from subspecies to their base species folder in {base_directory}")

    # Loop through all the directories in the base folder
    for species_folder in os.listdir(base_directory):
        species_folder_path = os.path.join(base_directory, species_folder)

        # Check if it's a directory
        if os.path.isdir(species_folder_path):
            # Split the folder name into species and subspecies components
            if species_folder == "Delphinus_bairdii":
                base_species_name = "Delphinus_delphis"
            else:
                species_name_parts = species_folder.split('_')
                base_species_name = '_'.join(species_name_parts[:2])  # Take the first two parts as the base species name

            # If the folder is a subspecies folder (not the base species folder)
            if base_species_name != species_folder:
                # Determine the path to the base species folder
                base_species_folder_path = os.path.join(base_directory, base_species_name)
                # Create the base species folder if it doesn't exist
                os.makedirs(base_species_folder_path, exist_ok=True)

                # Move all GB files from subspecies folders to the base species folder
                for file_name in os.listdir(species_folder_path):
                    if file_name.endswith(".gb"):
                        source_file = os.path.join(species_folder_path, file_name)
                        target_file = os.path.join(base_species_folder_path, file_name)

                        # Move the GB file if it doesn't already exist in the target folder
                        if not os.path.exists(target_file):
                            shutil.move(source_file, target_file)
                            print(f"Moved {file_name} from {species_folder} to {base_species_folder_path}")
                        else:
                            print(f"Skipping {file_name}, already exists in {base_species_folder_path}")

                # Remove the subspecies folder after all GB files are moved
                shutil.rmtree(species_folder_path)
                print(f"Removed folder {species_folder}")


if __name__ == "__main__":
    base_directory_fasta = "sequences/fasta"
    base_directory_gb = "sequences/gb"

    # Combine FASTA files from subspecies folders into base species folders
    combine_subspecies_fasta(base_directory_fasta)

    # Move GB files into their respective base species folders
    move_gb_files(base_directory_gb)

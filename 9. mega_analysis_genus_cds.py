import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import AlignIO


# The following function for converting AFA to FASTA is commented out since it's already done elsewhere
# def convert_afa_to_fasta(afa_path, fasta_path):
#     """
#     Convert alignment files in AFA format to FASTA format.
#     (Not used here because conversion is handled elsewhere)
#     """
#     # alignment = AlignIO.read(afa_path, "fasta")
#     # AlignIO.write(alignment, fasta_path, "fasta")

def create_grp_files(base_dir):
    """
    Create .grp files for each genus based on existing FASTA files.

    For each genus directory found in base_dir:
    - Finds the FASTA file named "*_concatenated_gapped_sequences.fasta".
    - For each sequence header in that FASTA file, it extracts the full species name
      from the last two underscore-separated parts.

    Examples:
    FASTA header: >MN145937.1_Balaena_mysticetus
    parts = ["MN145937.1", "Balaena", "mysticetus"]
    species_name = "Balaena_mysticetus"
    GRP line: MN145937.1_Balaena_mysticetus=Balaena_mysticetus(Balaena)

    Special case:
    FASTA header: >NC_005268.1_Balaena_mysticetus
    parts = ["NC", "005268.1", "Balaena", "mysticetus"]
    last two parts: ["Balaena", "mysticetus"]
    species_name = "Balaena_mysticetus"
    GRP line: NC_005268.1_Balaena_mysticetus=Balaena_mysticetus(Balaena)
    """
    for genus_folder in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus_folder, "CDS_nucleotide_gapped")

        # Only process if the path is a directory
        if not os.path.isdir(genus_dir):
            continue

        # We look for the main concatenated FASTA file
        for fasta_file in os.listdir(genus_dir):
            if fasta_file.endswith("_concatenated_gapped_sequences.fasta"):
                fasta_file_path = os.path.join(genus_dir, fasta_file)

                # Extract the genus name from the folder name
                genus = genus_folder

                # Create output directory for grp files
                genus_output_dir = os.path.join(base_dir, genus, "MEGA_Grouped_by_Species", "grp")
                os.makedirs(genus_output_dir, exist_ok=True)

                grp_file_path = os.path.join(genus_output_dir, f"{genus}.grp")

                with open(fasta_file_path, "r") as fasta, open(grp_file_path, "w") as grp:
                    for line in fasta:
                        if line.startswith(">"):
                            sequence_id = line[1:].strip()

                            # Split by underscore
                            parts = sequence_id.split("_")

                            # The last two parts should represent the species name: Genus_species
                            if len(parts) >= 2:
                                species_name = "_".join(parts[-2:])
                            else:
                                # Fallback if there's not enough parts (very unusual case)
                                species_name = parts[-1]

                            # Write the grp line: sequence_id=species_name(genus)
                            grp.write(f"{sequence_id}={species_name}({genus})\n")

                print(f"Created group file: {grp_file_path}")


def run_mega_analysis():
    """
    Run MEGA_Grouped_by_Species analyses using multiple settings files on a set of FASTA files organized by genus.

    Steps:
    - Checks for the MEGA_Grouped_by_Species executable.
    - Creates .grp files for each genus.
    - For each genus, identifies the concatenated FASTA file and its corresponding .grp file.
    - For each analysis type (defined by a .mao file), runs MEGA in parallel.
    - Prints the status of each analysis.
    """
    # Path to the MEGA command-line executable
    megacc_executable = r"C:\Program Files\MEGA11\megacc.exe"

    # Verify MEGA executable exists
    if not os.path.isfile(megacc_executable):
        print(f"Error: MEGA_Grouped_by_Species executable not found at {megacc_executable}")
        return

    # List of .mao setting files for different analyses
    settings_files = [
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_nucleotide.mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only).mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only).mao"
    ]

    # Base directory containing all genus folders
    base_dir = r"sequences/CDS_Genus"

    # Create .grp files for each genus
    create_grp_files(base_dir)

    tasks = []

    # Iterate over each genus directory
    for genus_folder in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus_folder, "CDS_nucleotide_gapped")

        if not os.path.isdir(genus_dir):
            continue

        # Identify concatenated FASTA file
        for fasta_file in os.listdir(genus_dir):
            if fasta_file.endswith("_concatenated_gapped_sequences.fasta"):
                fasta_file_path = os.path.join(genus_dir, fasta_file)

                genus = genus_folder
                genus_output_dir = os.path.join(base_dir, genus, "MEGA_Grouped_by_Species")
                os.makedirs(genus_output_dir, exist_ok=True)

                # Corresponding grp file
                group_file_path = os.path.join(genus_output_dir, "grp", f"{genus}.grp")

                # If no grp file, skip
                if not os.path.isfile(group_file_path):
                    print(f"Warning: Group file not found for genus '{genus}'")
                    continue

                # For each analysis type (mao file), create a task
                for mao_file in settings_files:
                    if not os.path.isfile(mao_file):
                        print(f"Warning: Settings file not found: {mao_file}")
                        continue

                    # Determine analysis type based on .mao filename
                    analysis_type = os.path.basename(mao_file).replace(".mao", "")
                    if analysis_type == "distance_estimation_within_grp_avg_nucleotide":
                        analysis_type = "d_value"
                    elif analysis_type == "distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only)":
                        analysis_type = "nonsynonymous_only"
                    elif analysis_type == "distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only)":
                        analysis_type = "Synonymous_only"

                    # Output file for this analysis
                    output_file = os.path.join(genus_output_dir, f"{genus}_{analysis_type}.meg")

                    # Build MEGA command
                    command = [
                        megacc_executable,
                        "-a", mao_file,  # Analysis settings
                        "-d", fasta_file_path,  # Input FASTA
                        "-g", group_file_path,  # Group file
                        "-o", output_file  # Output result
                    ]

                    tasks.append((genus, analysis_type, command))

    # If no tasks, nothing to run
    if not tasks:
        print("No tasks to run. Check input files and directories.")
        return

    # Run tasks in parallel
    with ThreadPoolExecutor(max_workers=12) as executor:
        future_to_task = {
            executor.submit(subprocess.run, task[2], check=True): task
            for task in tasks
        }

        for future in as_completed(future_to_task):
            genus, analysis_type, _ = future_to_task[future]
            try:
                future.result()  # Raises CalledProcessError if non-zero exit code
                print(f"Analysis '{analysis_type}' completed for genus '{genus}'.")
            except subprocess.CalledProcessError as e:
                print(f"Error running analysis '{analysis_type}' for genus '{genus}': {e}")


if __name__ == "__main__":
    run_mega_analysis()

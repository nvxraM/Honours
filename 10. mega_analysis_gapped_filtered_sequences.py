import subprocess
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio import AlignIO


# The following function for converting AFA to FASTA is commented out since it's already done elsewhere
# def convert_afa_to_fasta(afa_path, fasta_path):
#     """
#     Convert alignment files in AFA format to FASTA format.
#     """
#     alignment = AlignIO.read(afa_path, "fasta")
#     AlignIO.write(alignment, fasta_path, "fasta")

def run_mega_analysis():
    """
    Run MEGA analyses using multiple settings files on a set of FASTA files organized by genus.

    This function:
    - Checks for the MEGA executable.
    - Iterates through genus-specific directories containing gapped CDS FASTA files.
    - For each genus and each analysis setting (MAO file), it prepares a MEGA analysis command.
    - Executes these analyses in parallel using a ThreadPoolExecutor.
    - Prints the status and any errors encountered during execution.
    """
    # Define the path to the MEGA command-line executable
    megacc_executable = r"C:\Program Files\MEGA11\megacc.exe"

    # Check if the MEGA executable exists at the given path
    if not os.path.isfile(megacc_executable):
        print(f"Error: MEGA executable not found at {megacc_executable}")
        return

    # List of .mao setting files for different analyses
    # Each corresponds to a different type of distance estimation
    settings_files = [
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_nucleotide.mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only).mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only).mao"
    ]

    # Base directory containing all genus folders
    base_dir = r"sequences\CDS"

    # This list will hold tuples of (genus, analysis_type, command)
    # which will be used to run MEGA analyses in parallel
    tasks = []

    # Iterate over each genus directory within the base directory
    for genus_folder in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus_folder, "CDS_nucleotide_gapped")

        # Only process if the path is indeed a directory
        if not os.path.isdir(genus_dir):
            continue

        # We are looking specifically for files ending with "_concatenated_gapped_sequences.fasta"
        # as these represent aligned, concatenated sequences for the genus
        for fasta_file in os.listdir(genus_dir):
            if fasta_file.endswith("_concatenated_gapped_sequences.fasta"):
                fasta_file_path = os.path.join(genus_dir, fasta_file)

                # Extract the genus name from the folder name
                genus = genus_folder

                # Create an output directory for the genus if it doesn't exist
                genus_output_dir = os.path.join(base_dir, genus, "MEGA")
                os.makedirs(genus_output_dir, exist_ok=True)

                # The group file contains the grouping information for the sequences
                # It is assumed to be located in a directory named "A_Mega_NoGroups_grp"
                group_file_path = os.path.join(base_dir, "A_Mega_NoGroups_grp", f"{genus}_mds_groups.grp")

                # Check if the group file exists before proceeding
                if not os.path.isfile(group_file_path):
                    print(f"Warning: Group file not found for {genus}")
                    continue

                # Add a task for each of the settings (MAO) files
                for mao_file in settings_files:
                    # Ensure the settings file exists
                    if not os.path.isfile(mao_file):
                        print(f"Warning: Settings file not found: {mao_file}")
                        continue

                    # Determine the analysis type based on the name of the MAO file
                    analysis_type = os.path.basename(mao_file).replace(".mao", "")
                    if analysis_type == "distance_estimation_within_grp_avg_nucleotide":
                        analysis_type = "d_value"
                    elif analysis_type == "distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only)":
                        analysis_type = "nonsynonymous_only"
                    elif analysis_type == "distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only)":
                        analysis_type = "Synonymous_only"

                    # The output file will be a .meg file named according to the genus and analysis type
                    output_file = os.path.join(genus_output_dir, f"{genus}_{analysis_type}.meg")

                    # Build the command list for subprocess
                    command = [
                        megacc_executable,
                        "-a", mao_file,  # MEGA analysis settings file
                        "-d", fasta_file_path,  # Input FASTA file
                        "-g", group_file_path,  # Group file with grouping info
                        "-o", output_file  # Output file to store analysis results
                    ]

                    # Append the command and its associated metadata as a task
                    tasks.append((genus, analysis_type, command))

    # If no tasks were gathered, it means no suitable input was found
    if not tasks:
        print("No tasks to run. Please check the input files and directories.")
        return

    # Execute the tasks in parallel using up to 12 threads
    with ThreadPoolExecutor(max_workers=12) as executor:
        future_to_task = {
            executor.submit(subprocess.run, task[2], check=True): task
            for task in tasks
        }

        # Process the futures as they complete
        for future in as_completed(future_to_task):
            genus, analysis_type, _ = future_to_task[future]
            try:
                # future.result() will raise an exception if the command failed
                future.result()
                print(f"Analysis '{analysis_type}' completed for genus '{genus}'.")
            except subprocess.CalledProcessError as e:
                print(f"Error running MEGA analysis '{analysis_type}' for genus '{genus}': {e}")


# Execute the function if the script is run as the main program
if __name__ == "__main__":
    run_mega_analysis()

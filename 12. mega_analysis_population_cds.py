import subprocess
import os
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed

def create_grp_file_from_csv(csv_path):
    """
    Given a path to a CSV file that has "IID" and "Group" columns,
    create a .grp file if there is more than one unique group.
    Returns the path to the .grp file if created, otherwise None.
    """
    # Parse species and genus from the path:
    # Path structure: sequences/Interactive_group_editor/[Genus]/[Species]/[Species]_clustered.csv
    # We can split by os.sep and deduce genus and species
    parts = csv_path.split(os.sep)
    # parts[-3] should be the genus folder, parts[-2] the species folder, parts[-1] the csv file
    genus = parts[-3]
    species = parts[-2]

    # Read the CSV
    iids = []
    groups = []
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            iid = row["IID"]
            grp = row["Group"]
            iids.append(iid)
            groups.append(grp)

    unique_groups = set(groups)
    # If only one group, do not create .grp file
    if len(unique_groups) <= 1:
        return None

    # Create MEGA_Populations folder next to the CSV
    species_dir = os.path.dirname(csv_path)
    mega_pop_dir = os.path.join(species_dir, "MEGA_Populations")
    os.makedirs(mega_pop_dir, exist_ok=True)

    grp_filename = f"{species}_population.grp"
    grp_path = os.path.join(mega_pop_dir, grp_filename)

    # Write .grp file: each line "IID=species_Population_GroupNumber"
    with open(grp_path, "w", encoding='utf-8') as grp_file:
        for iid, grp in zip(iids, groups):
            grp_file.write(f"{iid}={species}_Population_{grp}\n")

    return grp_path


def find_csv_and_create_grp_files(base_dir):
    """
    Traverse the directory structure:
    sequences/Interactive_group_editor/[Genus]/[Species]/[Species]_clustered.csv

    For each such CSV, try to create a corresponding .grp file.
    Return a list of tuples: (genus, species, grp_file_path, fasta_file_path)
    for those species that have multiple groups and hence got a .grp file.
    """
    results = []
    for genus in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus)
        if not os.path.isdir(genus_dir):
            continue

        # Each species directory
        for species in os.listdir(genus_dir):
            species_dir = os.path.join(genus_dir, species)
            if not os.path.isdir(species_dir):
                continue

            # CSV file: [Species]_clustered.csv
            csv_file = os.path.join(species_dir, f"{species}_clustered.csv")
            if not os.path.isfile(csv_file):
                continue

            # FASTA file: [Species]_concatenated_gapped_sequences.fasta
            fasta_file = os.path.join(species_dir, f"{species}_concatenated_gapped_sequences.fasta")
            if not os.path.isfile(fasta_file):
                continue

            grp_file = create_grp_file_from_csv(csv_file)
            # Only append to results if grp_file was created (i.e., multiple groups)
            if grp_file is not None:
                results.append((genus, species, grp_file, fasta_file))
    return results


def run_mega_analysis(megacc_executable, settings_files, tasks):
    """
    Run MEGA analyses given a list of tasks.
    Each task: (genus, species, grp_file, fasta_file)
    For each .mao in settings_files, run analysis.
    If the output file already exists, skip to avoid re-running.
    """
    if not os.path.isfile(megacc_executable):
        print(f"Error: MEGA executable not found at {megacc_executable}")
        return

    futures = []
    with ThreadPoolExecutor(max_workers=12) as executor:
        for (genus, species, grp_file, fasta_file) in tasks:
            # Output directory (MEGA_Populations) is where grp_file is
            out_dir = os.path.dirname(grp_file)

            for mao_file in settings_files:
                if not os.path.isfile(mao_file):
                    print(f"Warning: Settings file not found: {mao_file}")
                    continue

                mao_basename = os.path.basename(mao_file).lower()

                # Determine analysis type and output file name based on .mao file name
                if "nucleotide" in mao_basename:
                    analysis_type = "d_value"
                    output_file = os.path.join(out_dir, f"{species}_{analysis_type}.meg")
                elif "nonsynonymous_only" in mao_basename and "between" not in mao_basename:
                    # For within-group nonsynonymous_only
                    analysis_type = "nonsynonymous_only"
                    output_file = os.path.join(out_dir, f"{species}_{analysis_type}.meg")
                elif "synonymous_only" in mao_basename and "between" not in mao_basename:
                    # For within-group Synonymous_only
                    analysis_type = "Synonymous_only"
                    output_file = os.path.join(out_dir, f"{species}_{analysis_type}.meg")
                elif "synonymous_only" in mao_basename and "between" in mao_basename:
                    # New case: Between-group Synonymous_only
                    analysis_type = "Synonymous_only_between"
                    # As requested, output as a .txt file
                    output_file = os.path.join(out_dir, f"{species}_Synonymous_only_between.txt")
                else:
                    # Generic fallback
                    analysis_type = mao_basename.replace(".mao", "")
                    output_file = os.path.join(out_dir, f"{species}_{analysis_type}.meg")

                # Check if output file already exists, skip if it does
                if os.path.exists(output_file):
                    print(f"Skipping analysis '{analysis_type}' for {genus}/{species} - output already exists.")
                    continue

                command = [
                    megacc_executable,
                    "-a", mao_file,
                    "-d", fasta_file,
                    "-g", grp_file,
                    "-o", output_file
                ]

                future = executor.submit(subprocess.run, command, check=True)
                futures.append((genus, species, analysis_type, future))

        # Collect results
        for genus, species, analysis_type, future in futures:
            try:
                future.result()
                print(f"Analysis '{analysis_type}' completed for {genus}/{species}.")
            except subprocess.CalledProcessError as e:
                print(f"Error running analysis '{analysis_type}' for {genus}/{species}: {e}")


if __name__ == "__main__":
    # Base directory with the structure:
    # sequences/Interactive_group_editor/[Genus]/[Species]/[Species]_clustered.csv
    base_dir = "sequences/Interactive_group_editor"

    # Path to the MEGA command-line executable
    megacc_executable = r"C:\Program Files\MEGA11\megacc.exe"

    # List of .mao setting files, including the new "between"-group file
    settings_files = [
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_nucleotide.mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only).mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only).mao",
        # New settings file for between-group Synonymous-only distances
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_between_grp_avg_syn-nonsynonymous(Synonymous_only).mao"
    ]

    # Create .grp files and gather tasks
    tasks = find_csv_and_create_grp_files(base_dir)

    # Run MEGA analysis for all tasks that have a grp file
    if tasks:
        run_mega_analysis(megacc_executable, settings_files, tasks)
    else:
        print("No .grp files created (only single-group species found or no CSV/FASTA pairs). No tasks to run.")

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
    parts = csv_path.split(os.sep)
    # parts[-3]: Genus folder, parts[-2]: Species folder, parts[-1]: CSV file
    genus = parts[-3]
    species = parts[-2]

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
    # If only one group, no .grp file needed
    if len(unique_groups) <= 1:
        return None

    # Create the MEGA_Populations folder
    species_dir = os.path.dirname(csv_path)
    mega_pop_dir = os.path.join(species_dir, "MEGA_Populations")
    os.makedirs(mega_pop_dir, exist_ok=True)

    grp_filename = f"{species}_population.grp"
    grp_path = os.path.join(mega_pop_dir, grp_filename)

    # Write .grp file
    with open(grp_path, "w", encoding='utf-8') as grp_file:
        for iid, grp in zip(iids, groups):
            grp_file.write(f"{iid}={species}_Population_{grp}\n")

    return grp_path


def find_csv_and_create_grp_files(base_dir):
    """
    Traverse the directory structure:
      sequences/Interactive_group_editor/[Genus]/[Species]/[Species]_clustered.csv
    For each CSV, create a .grp file if multiple groups exist.
    Return a list of tuples: (genus, species, grp_file_path, fasta_file_path)
    where grp_file is non-None.
    """
    results = []
    for genus in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus)
        if not os.path.isdir(genus_dir):
            continue

        for species in os.listdir(genus_dir):
            species_dir = os.path.join(genus_dir, species)
            if not os.path.isdir(species_dir):
                continue

            csv_file = os.path.join(species_dir, f"{species}_clustered.csv")
            if not os.path.isfile(csv_file):
                continue

            fasta_file = os.path.join(species_dir, f"{species}_concatenated_gapped_sequences.fasta")
            if not os.path.isfile(fasta_file):
                continue

            grp_file = create_grp_file_from_csv(csv_file)
            if grp_file is not None:
                results.append((genus, species, grp_file, fasta_file))
    return results


def run_mega_analysis(megacc_executable, settings_files, tasks):
    """
    Run MEGA-CC analysis for each task and each settings file.
    Skip analysis if output file already exists.
    """
    if not os.path.isfile(megacc_executable):
        print(f"Error: MEGA executable not found at {megacc_executable}")
        return

    futures = []
    with ThreadPoolExecutor(max_workers=12) as executor:
        for (genus, species, grp_file, fasta_file) in tasks:
            out_dir = os.path.dirname(grp_file)

            for mao_file in settings_files:
                if not os.path.isfile(mao_file):
                    print(f"Warning: Settings file not found: {mao_file}")
                    continue

                mao_basename = os.path.basename(mao_file).lower()

                # We only have one .mao: it is the between-group NonSynonymous-only
                # But let's keep the logic to parse the analysis type
                if "nonsynonymous_only" in mao_basename and "between" in mao_basename:
                    analysis_type = "NonSynonymous_only_between"
                    output_file = os.path.join(out_dir, f"{species}_NonSynonymous_only_between.meg")
                else:
                    # Fallback, though we'll never hit this if we only have the single .mao
                    analysis_type = mao_basename.replace(".mao", "")
                    output_file = os.path.join(out_dir, f"{species}_{analysis_type}.meg")

                # Skip if output already exists
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

        # Gather results
        for genus, species, analysis_type, future in futures:
            try:
                future.result()
                print(f"Analysis '{analysis_type}' completed for {genus}/{species}.")
            except subprocess.CalledProcessError as e:
                print(f"Error running analysis '{analysis_type}' for {genus}/{species}: {e}")


if __name__ == "__main__":
    # Set base directory where your [Genus]/[Species] folders live
    base_dir = "sequences/Interactive_group_editor"

    # Path to MEGA-CC executable
    megacc_executable = r"C:\Program Files\MEGA11\megacc.exe"

    # Only one settings file: between-group NonSynonymous-only
    settings_files = [
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_between_grp_avg_syn-nonsynonymous(nonSynonymous_only).mao"
    ]

    # Identify which species have multiple groups => create tasks
    tasks = find_csv_and_create_grp_files(base_dir)

    # If no multi-group species found, skip
    if not tasks:
        print("No .grp files created (only single-group species found or no CSV/FASTA pairs). No tasks to run.")
    else:
        # Run only the single requested analysis
        run_mega_analysis(megacc_executable, settings_files, tasks)

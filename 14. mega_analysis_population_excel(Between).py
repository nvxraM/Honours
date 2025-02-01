import os
import re
import sys
import subprocess

# 1) Ensure xlsxwriter is installed (to avoid "ModuleNotFoundError: No module named 'xlsxwriter'")
try:
    import xlsxwriter
except ImportError:
    print("Module 'xlsxwriter' not found; installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "xlsxwriter"])
    import xlsxwriter

import pandas as pd

BASE_DIR = 'sequences/Interactive_group_editor'
# Original output
OUTPUT_XLSX_SYN = 'sequences/MEGA/All_Synonymous_only_dS_matrices_stacked.xlsx'
# New output
OUTPUT_XLSX_NON_SYN = 'sequences/MEGA/All_nonSynonymous_only_dN_matrices_stacked.xlsx'


def parse_mega_between_group(filepath):
    """
    Parses a MEGA 'between-group' file in lower-left format like:

        [1] #Balaena_mysticetus_Population_1
        [2] #Balaena_mysticetus_Population_2
        [3] #Balaena_mysticetus_Population_3

        [               1             2             3 ]
        [1]                [ 0.00125518 ][ 0.00133804 ]
        [2]    0.00802042                [ 0.00161139 ]
        [3]    0.00941509    0.01092384

    Returns an N×N DataFrame capturing EXACTLY the original matrix values,
    preserving any asymmetry (i.e. not forcing symmetry).
    Rows/columns are labeled with the population names.
    """
    label_pattern = re.compile(r'^\[(\d+)\]\s*#(.+)$')
    matrix_line_pattern = re.compile(r'^\[(\d+)\]')

    with open(filepath, 'r', encoding='utf-8') as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    # Extract bracket-index => label, e.g. {1: "Balaena_mysticetus_Population_1", ...}
    pop_labels = {}
    for line in lines:
        m = label_pattern.match(line)
        if m:
            idx_str, label_str = m.groups()
            idx = int(idx_str)
            pop_labels[idx] = label_str

    n = len(pop_labels)
    if n == 0:
        return pd.DataFrame()

    dist_matrix = [[0.0] * n for _ in range(n)]

    for line in lines:
        mat_match = matrix_line_pattern.match(line)
        if mat_match:
            row_idx = int(mat_match.group(1))  # bracket index, 1-based
            remainder = line[line.index(']') + 1:].strip()

            # Find all numeric floats
            numbers = re.findall(r'([\d.]+)', remainder)

            left_count = row_idx - 1
            right_count = n - row_idx

            left_vals = numbers[:left_count]
            right_vals = numbers[left_count:left_count + right_count]

            for col_i, val_str in enumerate(left_vals, start=1):
                dist_matrix[row_idx - 1][col_i - 1] = float(val_str)

            for offset, val_str in enumerate(right_vals, start=1):
                col_idx = row_idx + offset
                dist_matrix[row_idx - 1][col_idx - 1] = float(val_str)

    sorted_indices = sorted(pop_labels.keys())
    label_list = [pop_labels[i] for i in sorted_indices]
    df = pd.DataFrame(dist_matrix, index=label_list, columns=label_list)
    return df


def main():
    # We'll gather all species data for synonyms and non-synonyms separately
    all_results_syn = []
    all_results_non_syn = []

    for genus in os.listdir(BASE_DIR):
        genus_path = os.path.join(BASE_DIR, genus)
        if not os.path.isdir(genus_path):
            continue

        for species_dir_name in os.listdir(genus_path):
            species_path = os.path.join(genus_path, species_dir_name)
            if not os.path.isdir(species_path):
                continue

            mega_pop_path = os.path.join(species_path, "MEGA_Populations")
            if not os.path.isdir(mega_pop_path):
                continue

            # 1) Look for Synonymous
            syn_file = os.path.join(mega_pop_path, f"{species_dir_name}_Synonymous_only_between.meg")
            if os.path.exists(syn_file):
                df_syn = parse_mega_between_group(syn_file)
                if not df_syn.empty:
                    all_results_syn.append((genus, species_dir_name, df_syn))

            # 2) Look for NonSynonymous
            non_syn_file = os.path.join(mega_pop_path, f"{species_dir_name}_NonSynonymous_only_between.meg")
            if os.path.exists(non_syn_file):
                df_non_syn = parse_mega_between_group(non_syn_file)
                if not df_non_syn.empty:
                    all_results_non_syn.append((genus, species_dir_name, df_non_syn))

    # Ensure output folder (sequences/MEGA) exists
    os.makedirs(os.path.dirname(OUTPUT_XLSX_SYN), exist_ok=True)

    # -------------------------------------------------------------------
    # Write the Synonymous stacked Excel (always overwrite, as original)
    # -------------------------------------------------------------------
    if all_results_syn:
        with pd.ExcelWriter(OUTPUT_XLSX_SYN, engine='xlsxwriter') as writer:
            workbook = writer.book
            sheet_name = 'Synonymous_only'
            worksheet = workbook.add_worksheet(sheet_name)
            writer.sheets[sheet_name] = worksheet

            row_offset = 0

            for (genus, species, df) in all_results_syn:
                # Write a header row
                worksheet.write(row_offset, 0, f"{genus} - {species}")
                row_offset += 1

                # Write the NxN table
                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    startrow=row_offset,
                    startcol=0,
                    index=True,
                    header=True
                )
                # Add blank gap of 3 rows after
                row_offset += len(df.index) + 1 + 3

        print(f"Done! Stacked Excel of all Synonymous data: {OUTPUT_XLSX_SYN}")
    else:
        print("No '_Synonymous_only_between.meg' files found or no valid data parsed.")

    # -------------------------------------------------------------------
    # Write the NonSynonymous stacked Excel
    # but CHECK if the file already exists first
    # -------------------------------------------------------------------
    if os.path.exists(OUTPUT_XLSX_NON_SYN):
        print(f"Output '{OUTPUT_XLSX_NON_SYN}' already exists. Skipping creation.")
    else:
        # Only create if we have data
        if all_results_non_syn:
            with pd.ExcelWriter(OUTPUT_XLSX_NON_SYN, engine='xlsxwriter') as writer:
                workbook = writer.book
                sheet_name = 'NonSynonymous_only'
                worksheet = workbook.add_worksheet(sheet_name)
                writer.sheets[sheet_name] = worksheet

                row_offset = 0

                for (genus, species, df) in all_results_non_syn:
                    worksheet.write(row_offset, 0, f"{genus} - {species}")
                    row_offset += 1

                    df.to_excel(
                        writer,
                        sheet_name=sheet_name,
                        startrow=row_offset,
                        startcol=0,
                        index=True,
                        header=True
                    )
                    row_offset += len(df.index) + 1 + 3

            print(f"Done! Stacked Excel of all NonSynonymous data: {OUTPUT_XLSX_NON_SYN}")
        else:
            print("No '_NonSynonymous_only_between.meg' files found or no valid data parsed.")


if __name__ == "__main__":
    main()

import os
import re
import pandas as pd

def parse_fst_matrix(filepath):
    """Parse the Fst matrix from a DnaSP RAD.PW.out file."""
    with open(filepath, encoding='utf-8') as f:
        lines = f.readlines()
    # Find the Fst section
    fst_start = None
    for i, line in enumerate(lines):
        if line.strip() == "Fst":
            fst_start = i
            break
    if fst_start is None:
        raise ValueError(f"Fst section not found in {filepath}")
    # Get the header (population names)
    header_line = lines[fst_start + 1].strip()
    while not header_line:  # skip blanks
        fst_start += 1
        header_line = lines[fst_start + 1].strip()
    col_names = re.split(r'\t+', header_line)
    col_names = [x for x in col_names if x]
    data = []
    row_names = []
    row_i = fst_start + 2
    while row_i < len(lines):
        line = lines[row_i].strip()
        if not line or not any(char.isdigit() for char in line):
            break
        parts = re.split(r'\t+', line)
        row_names.append(parts[0])
        row_data = parts[1:]
        data.append(row_data)
        row_i += 1
    # Build lower-triangle DataFrame
    n = len(col_names)
    df = pd.DataFrame('', index=row_names, columns=col_names)
    for i, row in enumerate(data):
        for j, val in enumerate(row):
            col_idx = i + j  # Lower triangle: shift right by row index
            if val.strip():
                df.iloc[i, col_idx] = float(val)
    return df

base_dir = "sequences/Species_POP_Mega_Haploid"
output_dir = "sequences/Species_POP_Mega_CSVOutput"
os.makedirs(output_dir, exist_ok=True)
species_results = []
long_records = []

for species_folder in sorted(os.listdir(base_dir)):
    result_dir = os.path.join(base_dir, species_folder, "Results")
    if not os.path.isdir(result_dir):
        continue
    for fname in os.listdir(result_dir):
        if fname.endswith('.RAD.PW.out'):
            path = os.path.join(result_dir, fname)
            try:
                fst_df = parse_fst_matrix(path)
                species_results.append((species_folder, fst_df))
                # --- Collect for long-format
                for i, pop1 in enumerate(fst_df.index):
                    for j, pop2 in enumerate(fst_df.columns):
                        val = fst_df.iloc[i, j]
                        if isinstance(val, float):  # Only keep actual Fst values
                            long_records.append({
                                "Species": species_folder,
                                "Population 1": pop1,
                                "Population 2": pop2,
                                "Fst": val
                            })
            except Exception as e:
                print(f"Could not parse {path}: {e}")

# ----- First Sheet: Stacked tables -----
stacked_rows = []
for species, df in species_results:
    stacked_rows.append([species])  # species title as a separator row
    stacked_rows.append([''] + list(df.columns))
    for idx, row in df.iterrows():
        stacked_rows.append([idx] + list(row.values))
    stacked_rows.append([''])  # empty row between species

stacked_df = pd.DataFrame(stacked_rows)

# ----- Second Sheet: Long format -----
long_df = pd.DataFrame(long_records)

# ----- Write to Excel -----
excel_path = os.path.join(output_dir, "Fst_results.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    stacked_df.to_excel(writer, sheet_name="All Fst tables", header=False, index=False)
    long_df.to_excel(writer, sheet_name="Long format", index=False)

print(f"Results written to {excel_path}")

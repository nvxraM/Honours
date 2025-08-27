import os
import re
import pandas as pd

# ---------- Helpers ----------

# Common-name map (extend anytime)
COMMON = {
    "Balaena_mysticetus": "Bowhead whale",
    "Balaenoptera_musculus": "Blue whale",
    "Balaenoptera_physalus": "Fin whale",
    "Delphinapterus_leucas": "Beluga whale",
    "Eubalaena_japonica": "North Pacific right whale",
    "Globicephala_macrorhynchus": "Short-finned pilot whale",
    "Hyperoodon_ampullatus": "Northern bottlenose whale",
    "Mesoplodon_grayi": "Gray's beaked whale",
    "Monodon_monoceros": "Narwhal",
    "Orcaella_brevirostris": "Irrawaddy dolphin",
    "Orcinus_orca": "Killer whale",
    "Peponocephala_electra": "Melon-headed whale",
    "Phocoena_phocoena": "Harbour porpoise",
    "Phocoena_sinus": "Vaquita",
    "Physeter_macrocephalus": "Sperm whale",
    "Pseudorca_crassidens": "False killer whale",
    "Stenella_longirostris": "Spinner dolphin",
    "Tursiops_aduncus": "Indo-Pacific bottlenose dolphin",
    "Tursiops_truncatus": "Common bottlenose dolphin",
    "Ziphius_cavirostris": "Cuvier's beaked whale",
}

ORDINAL_TO_NUM = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10
}

def abbrev_scientific(species_key: str) -> str:
    """'Balaena_mysticetus' -> 'B. mysticetus'"""
    genus, epithet = species_key.split("_", 1)
    return f"{genus[0]}. {epithet}"

def species_label(species_key: str) -> str:
    """'Balaena_mysticetus' -> 'Bowhead whale (B. mysticetus)'"""
    common = COMMON.get(
        species_key, species_key.replace("_", " ").title()
    )
    return f"{common} ({abbrev_scientific(species_key)})"

def pop_to_gp(pop: str) -> tuple[str, int]:
    """
    Extract group label and numeric index.
    '..._Population_two' -> ('Gp2', 2)
    """
    m = re.search(r"Population_(one|two|three|four|five|six|seven|eight|nine|ten)", pop, re.I)
    if not m:
        return ("Gp?", 0)
    num = ORDINAL_TO_NUM[m.group(1).lower()]
    return (f"Gp{num}", num)

def comparison_label(p1: str, p2: str) -> tuple[str, int, int]:
    """
    Build 'GpX & GpY' with larger group first (e.g., 'Gp3 & Gp1').
    Also return numeric parts for sorting.
    """
    gp1, n1 = pop_to_gp(p1)
    gp2, n2 = pop_to_gp(p2)
    if n1 >= n2:
        return f"{gp1} & {gp2}", n1, n2
    else:
        return f"{gp2} & {gp1}", n2, n1

# ---------- Fst parsing ----------

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

    # Header (population names)
    header_line = lines[fst_start + 1].strip()
    while not header_line:
        fst_start += 1
        header_line = lines[fst_start + 1].strip()
    col_names = [x for x in re.split(r'\t+', header_line) if x]

    # Rows (lower-triangle numbers)
    data, row_names = [], []
    row_i = fst_start + 2
    while row_i < len(lines):
        line = lines[row_i].strip()
        if not line or not any(ch.isdigit() for ch in line):
            break
        parts = re.split(r'\t+', line)
        row_names.append(parts[0])
        data.append(parts[1:])
        row_i += 1

    n = len(col_names)
    df = pd.DataFrame('', index=row_names, columns=col_names)
    for i, row in enumerate(data):
        for j, val in enumerate(row):
            col_idx = i + j  # lower triangle
            if val.strip():
                df.iloc[i, col_idx] = float(val)
    return df

# ---------- Main ----------

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
                for i, pop1 in enumerate(fst_df.index):
                    for j, pop2 in enumerate(fst_df.columns):
                        val = fst_df.iloc[i, j]
                        if isinstance(val, float):
                            long_records.append({
                                "Species": species_folder,
                                "Population 1": pop1,
                                "Population 2": pop2,
                                "Fst": val
                            })
            except Exception as e:
                print(f"Could not parse {path}: {e}")

# ----- Sheet 1: Stacked tables (unchanged) -----
stacked_rows = []
for species, df in species_results:
    stacked_rows.append([species])
    stacked_rows.append([''] + list(df.columns))
    for idx, row in df.iterrows():
        stacked_rows.append([idx] + list(row.values))
    stacked_rows.append([''])
stacked_df = pd.DataFrame(stacked_rows)

# ----- Sheet 2: Long format (unchanged) -----
long_df = pd.DataFrame(long_records)

# ----- Sheet 3: Simple view -----
simple_rows = []
for rec in long_records:
    sp_key = rec["Species"]
    comp_lbl, hi, lo = comparison_label(rec["Population 1"], rec["Population 2"])
    simple_rows.append({
        "Combined Common": species_label(sp_key),
        "Comparisons between populations": comp_lbl,
        "FST": rec["Fst"],
        "_sort_hi": hi,     # helper for ordering
        "_sort_lo": lo
    })

simple_df = pd.DataFrame(simple_rows)
# Sort by species, then by (higher GP desc, lower GP asc) to get: Gp2&Gp1, Gp3&Gp1, Gp3&Gp2, ...
simple_df = simple_df.sort_values(
    by=["Combined Common", "_sort_hi", "_sort_lo"],
    ascending=[True, False, True]
).drop(columns=["_sort_hi", "_sort_lo"])

# ----- Write to Excel -----
excel_path = os.path.join(output_dir, "Fst_results.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    stacked_df.to_excel(writer, sheet_name="All Fst tables", header=False, index=False)
    long_df.to_excel(writer, sheet_name="Long format", index=False)
    simple_df.to_excel(writer, sheet_name="Simple view", index=False)

print(f"Results written to {excel_path}")

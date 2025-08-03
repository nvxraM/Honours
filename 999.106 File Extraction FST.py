import re
import pandas as pd
from pathlib import Path

def extract_fst_table(text):
    """Extract Fst table as a pandas DataFrame from DnaSP output text."""
    # Find the Fst table section
    match = re.search(r"Fst\s*\n((?:[^\n]*\n)+)", text)
    if not match:
        return None

    lines = match.group(1).strip().split('\n')
    if not lines:
        return None

    # The first row is header
    headers = [x.strip() for x in lines[0].split('\t') if x.strip()]
    n_pops = len(headers)

    data = []
    row_names = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        row_name = parts[0].strip()
        values = parts[1:]
        row_names.append(row_name)
        data.append(values)

    # Build a DataFrame (lower triangle)
    df = pd.DataFrame(data, columns=headers, index=row_names)
    df = df.replace("", float("nan")).astype(float)
    return df

def melt_fst_table(df, species):
    """Turn a lower triangle Fst table DataFrame into long form (species, pop1, pop2, fst)."""
    melted = []
    for i, pop1 in enumerate(df.index):
        for j, pop2 in enumerate(df.columns):
            # Only keep lower triangle (i < j, as DnaSP reports)
            if j > i:
                continue
            if pd.notnull(df.iloc[i, j]):
                melted.append({
                    "Species": species,
                    "Population_1": pop1,
                    "Population_2": pop2,
                    "Fst": df.iloc[i, j]
                })
    return melted

root = Path("sequences/Species_POP_Mega_Haploid")
output_dir = Path("sequences/Species_POP_Mega_CSVOutput")
output_dir.mkdir(exist_ok=True, parents=True)
all_results = []

# For every species
for species_folder in root.iterdir():
    if not species_folder.is_dir():
        continue
    species = species_folder.name
    results_dir = species_folder / "Results"
    for file in results_dir.glob("*.RAD.PW.out"):
        with open(file) as f:
            text = f.read()
        fst_df = extract_fst_table(text)
        if fst_df is not None:
            all_results.extend(melt_fst_table(fst_df, species))

# Final long dataframe
df = pd.DataFrame(all_results)
df = df[["Species", "Population_1", "Population_2", "Fst"]]

excel_out = output_dir / "All_species_pairwise_Fst.xlsx"
df.to_excel(excel_out, index=False)
print(f"All Fst results saved to {excel_out}")

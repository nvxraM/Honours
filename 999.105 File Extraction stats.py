import re
import pandas as pd
from pathlib import Path

def extract_fields_from_dnasp_out(text, filename, species):
    # Population: prefer from file content, fallback to filename if not found
    pop_match = re.search(r'Population used:\s*(\S+)', text)
    if pop_match:
        population = pop_match.group(1)
    else:
        m = re.search(r'Population_(\d+)', filename)
        population = f"Population_{m.group(1)}" if m else ""
    # Patterns for each field
    patterns = {
        "Number of Sequences used": r"Number of sequences used:\s*(\d+)",
        "Variable sites (VS)": r"Number of polymorphic \(segregating\) sites, S:\s*(\d+)",
        "Haplotype number (HN)": r"Number of Haplotypes, h:\s*(\d+)",
        "Haplotype diversity (Hd)": r"Haplotype \(gene\) diversity, Hd:\s*([0-9.]+)",
        "Haplotype diversity Standard Deviation (Hd.SD)": r"Standard Deviation of Haplotype diversity:\s*([0-9.]+)",
        "Nucleotide diversity (Pi)": r"Nucleotide diversity, Pi:\s*([0-9.]+)",
        "Nucleotide diversity standard deviation (Pi.SD)": r"Standard deviation of Pi:\s*([0-9.]+)",
    }
    row = {
        "Species": species,
        "Population": population
    }
    for field, pat in patterns.items():
        m = re.search(pat, text)
        row[field] = m.group(1) if m else ""
    return row

# Input and output directories
root_dir = Path("sequences/Species_POP_Mega_Haploid_split")
output_dir = Path("sequences/Species_POP_Mega_CSVOutput")
output_dir.mkdir(parents=True, exist_ok=True)

rows = []
for species_folder in root_dir.iterdir():
    if not species_folder.is_dir():
        continue
    species = species_folder.name
    for out_file in species_folder.glob("*.out"):
        with open(out_file) as f:
            text = f.read()
        row = extract_fields_from_dnasp_out(text, out_file.name, species)
        rows.append(row)

# Use your required order for columns
columns = [
    "Species",
    "Population",
    "Number of Sequences used",
    "Variable sites (VS)",
    "Haplotype number (HN)",
    "Haplotype diversity (Hd)",
    "Haplotype diversity Standard Deviation (Hd.SD)",
    "Nucleotide diversity (Pi)",
    "Nucleotide diversity standard deviation (Pi.SD)"
]

df = pd.DataFrame(rows)[columns]
excel_out = output_dir / "All_species_population_summary.xlsx"
df.to_excel(excel_out, index=False)
print(f"Done! Output written to {excel_out}")

import re
import pandas as pd
from pathlib import Path

# ---------- Config ----------
ROOT_DIR = Path("sequences/Species_POP_Mega_Haploid_split")
OUTPUT_DIR = Path("sequences/Species_POP_Mega_CSVOutput")
OUTPUT_XLSX = OUTPUT_DIR / "All_species_population_summary.xlsx"

COLUMNS = [
    "Species",
    "Common name",
    "Family",
    "Population",
    "Combined Common (Pop)",  # <-- new column
    "Number of Sequences used",
    "Variable sites (VS)",
    "Haplotype number (HN)",
    "Haplotype diversity (Hd)",
    "Haplotype diversity Standard Deviation (Hd.SD)",
    "Nucleotide diversity (Pi)",
    "Nucleotide diversity standard deviation (Pi.SD)",
]

# Keep integers as integers; preserve decimal fields as TEXT (exact strings)
INT_COLS = ["Number of Sequences used", "Variable sites (VS)", "Haplotype number (HN)"]
TEXT_DECIMAL_COLS = [
    "Haplotype diversity (Hd)",
    "Haplotype diversity Standard Deviation (Hd.SD)",
    "Nucleotide diversity (Pi)",
    "Nucleotide diversity standard deviation (Pi.SD)",
]

# Decimal pattern that captures everything incl. scientific notation
DEC = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"

COMMON_NAME = {
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
    "Orcinus_orca": "Killer whale (Orca)",
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

FAMILY_BY_SPECIES = {
    "Balaena_mysticetus": "Balaenidae",
    "Eubalaena_japonica": "Balaenidae",
    "Balaenoptera_musculus": "Balaenopteridae",
    "Balaenoptera_physalus": "Balaenopteridae",
    "Delphinapterus_leucas": "Monodontidae",
    "Monodon_monoceros": "Monodontidae",
    "Globicephala_macrorhynchus": "Delphinidae",
    "Orcaella_brevirostris": "Delphinidae",
    "Orcinus_orca": "Delphinidae",
    "Peponocephala_electra": "Delphinidae",
    "Pseudorca_crassidens": "Delphinidae",
    "Stenella_longirostris": "Delphinidae",
    "Tursiops_aduncus": "Delphinidae",
    "Tursiops_truncatus": "Delphinidae",
    "Phocoena_phocoena": "Phocoenidae",
    "Phocoena_sinus": "Phocoenidae",
    "Physeter_macrocephalus": "Physeteridae",
    "Hyperoodon_ampullatus": "Ziphiidae",
    "Mesoplodon_grayi": "Ziphiidae",
    "Ziphius_cavirostris": "Ziphiidae",
}

FAMILY_BY_GENUS = {
    "Balaena": "Balaenidae", "Eubalaena": "Balaenidae",
    "Balaenoptera": "Balaenopteridae",
    "Delphinapterus": "Monodontidae", "Monodon": "Monodontidae",
    "Globicephala": "Delphinidae", "Orcaella": "Delphinidae",
    "Orcinus": "Delphinidae", "Peponocephala": "Delphinidae",
    "Pseudorca": "Delphinidae", "Stenella": "Delphinidae",
    "Tursiops": "Delphinidae",
    "Phocoena": "Phocoenidae",
    "Physeter": "Physeteridae",
    "Hyperoodon": "Ziphiidae", "Mesoplodon": "Ziphiidae", "Ziphius": "Ziphiidae",
}

WORD_TO_NUM = {
    "one": "1", "two": "2", "three": "3", "four": "4", "five": "5",
    "six": "6", "seven": "7", "eight": "8", "nine": "9", "ten": "10"
}

def map_population_to_group(pop_text: str) -> str:
    if not pop_text:
        return ""
    m = re.search(r"\bPopulation_(\d+|one|two|three|four|five|six|seven|eight|nine|ten)\b",
                  pop_text, flags=re.IGNORECASE)
    if not m:
        return pop_text
    token = m.group(1)
    num_str = WORD_TO_NUM.get(token.lower(), token)
    return f"Gp{num_str}"

def common_name_for(scientific: str) -> str:
    return COMMON_NAME.get(scientific, scientific.replace("_", " ").title())

def family_for(scientific: str) -> str:
    if scientific in FAMILY_BY_SPECIES:
        return FAMILY_BY_SPECIES[scientific]
    genus = scientific.split("_", 1)[0]
    return FAMILY_BY_GENUS.get(genus, "Unknown")

def extract_fields_from_dnasp_out(text: str, filename: str, species: str) -> dict:
    # Population from content or filename
    pop_match = re.search(r"Population used:\s*(\S+)", text)
    if pop_match:
        population_raw = pop_match.group(1)
    else:
        mfile = re.search(r"(Population_(\d+|one|two|three|four|five|six|seven|eight|nine|ten))",
                          filename, re.IGNORECASE)
        population_raw = mfile.group(1) if mfile else ""
    population = map_population_to_group(population_raw)

    # Capture EXACT strings for decimals (no conversion!)
    patterns = {
        "Number of Sequences used": r"Number of sequences used:\s*(\d+)",
        "Variable sites (VS)": r"Number of polymorphic \(segregating\) sites, S:\s*(\d+)",
        "Haplotype number (HN)": r"Number of Haplotypes, h:\s*(\d+)",
        "Haplotype diversity (Hd)": rf"Haplotype \(gene\) diversity, Hd:\s*{DEC}",
        "Haplotype diversity Standard Deviation (Hd.SD)": rf"Standard Deviation of Haplotype diversity:\s*{DEC}",
        "Nucleotide diversity (Pi)": rf"Nucleotide diversity, Pi:\s*{DEC}",
        "Nucleotide diversity standard deviation (Pi.SD)": rf"Standard deviation of Pi:\s*{DEC}",
    }

    row = {
        "Species": species,
        "Common name": common_name_for(species),
        "Family": family_for(species),
        "Population": population,
    }
    for field, pat in patterns.items():
        m2 = re.search(pat, text)
        row[field] = m2.group(1) if m2 else ""
    return row

def _combine_common_pop(common: str, pop: str) -> str:
    common = ("" if pd.isna(common) else str(common).strip())
    pop = ("" if pd.isna(pop) else str(pop).strip())
    if common and pop:
        return f"{common} ({pop})"
    return common or pop

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for species_folder in ROOT_DIR.iterdir():
        if not species_folder.is_dir():
            continue
        species = species_folder.name
        for out_file in species_folder.glob("*.out"):
            with open(out_file, "r", encoding="utf-8", errors="ignore") as f:
                text = f.read()
            rows.append(extract_fields_from_dnasp_out(text, out_file.name, species))

    df = pd.DataFrame(rows)

    # Build the combined label BEFORE reindexing to COLUMNS
    df["Combined Common (Pop)"] = df.apply(
        lambda r: _combine_common_pop(r.get("Common name", ""), r.get("Population", "")),
        axis=1
    )

    # Ensure all expected columns exist, then order them
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df = df[COLUMNS]

    # Integers as integers; leave the decimal columns as TEXT (exact strings)
    if INT_COLS:
        df[INT_COLS] = df[INT_COLS].apply(pd.to_numeric, errors="coerce").astype("Int64")

    # Write Excel and force decimal columns + combined column to text so Excel shows all digits verbatim
    try:
        with pd.ExcelWriter(
            OUTPUT_XLSX,
            engine="xlsxwriter",
            engine_kwargs={"options": {"strings_to_numbers": False}},  # don't auto-convert
        ) as writer:
            df.to_excel(writer, index=False, sheet_name="Summary")
            wb = writer.book
            ws = writer.sheets["Summary"]

            # Autosize columns
            for i, col in enumerate(df.columns):
                ws.set_column(i, i, max(14, min(48, len(col) + 2)))

            # Integer columns as 0-decimal numbers
            int_fmt = wb.add_format({"num_format": "0"})
            for col_name in INT_COLS:
                col_idx = df.columns.get_loc(col_name)
                ws.set_column(col_idx, col_idx, None, int_fmt)

            # Text format for decimal columns + the new combined column
            text_fmt = wb.add_format({"num_format": "@"})
            for col_name in TEXT_DECIMAL_COLS + ["Combined Common (Pop)"]:
                col_idx = df.columns.get_loc(col_name)
                ws.set_column(col_idx, col_idx, None, text_fmt)

    except Exception as e:
        # Fallback: still writes strings; Excel may display 'General'
        df.to_excel(OUTPUT_XLSX, index=False)
        print(f"Note: wrote without column formats ({e})")

    print(f"Done! Output written to {OUTPUT_XLSX}")

if __name__ == "__main__":
    main()

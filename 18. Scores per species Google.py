import os
import glob
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

# --- Configuration ---
# Path to the reference sequence FASTA file (Humpback Whale)
REFERENCE_FASTA_PATH = "sequences/CDS_Reference/CDS_nucleotide_gapped/CDS_Reference_concatenated_gapped_sequences.fasta"
# Directory containing FASTA files for other species (grouped by Genus)
COMPARISON_SEQUENCES_DIR = "sequences/CDS_Genus/CDS_nucleotide_gapped/"
# Output Excel file path
OUTPUT_EXCEL_PATH = "sequences/scores/scores.xlsx"  # Using .xlsx for modern Excel

# --- Mappings and Data ---

# Scientific name of the reference species (Humpback Whale)
# This is used if FASTA headers for reference provide scientific name
REFERENCE_SPECIES_COMMON_NAME = "Hump back whale"  # As per your text
# (Megaptera novaeangliae is the scientific name, but not strictly needed if ref file is clearly identified)

# Map scientific names (from your table) to common names (for FASTA parsing and output)
# This is crucial and needs to be accurate for your dataset.
# Add all species you expect to encounter.
SCIENTIFIC_TO_COMMON_MAP = {
    "Balaena mysticetus": "Bowhead Whale",
    "Balaenoptera musculus": "Blue Whale",
    "Balaenoptera physalus": "Fin Whale",
    "Delphinapterus leucas": "Beluga Whale",
    "Delphinus delphis": "Short-beaked Common Dolphin",  # Assuming common species
    "Eubalaena australis": "Southern Right Whale",
    "Eubalaena glacialis": "North Atlantic Right Whale",
    "Eubalaena japonica": "North Pacific Right Whale",
    "Globicephala macrorhynchus": "Short-finned Pilot Whale",
    "Hyperoodon ampullatus": "Northern Bottlenose Whale",
    "Mesoplodon densirostris": "Blainville's Beaked Whale",
    "Mesoplodon europaeus": "Gervais' Beaked Whale",
    "Mesoplodon grayi": "Gray's Beaked Whale",
    "Mesoplodon mirus": "True's Beaked Whale",
    "Monodon monoceros": "Narwhal",
    "Neophocaena asiaeorientalis": "Yangtze Finless Porpoise",
    "Neophocaena phocaenoides": "Indo-Pacific Finless Porpoise",
    "Orcaella brevirostris": "Irrawaddy Dolphin",
    "Orcinus orca": "Orca",
    "Peponocephala electra": "Melon-headed Whale",
    "Phocoena phocoena": "Harbour Porpoise",
    "Phocoena sinus": "Vaquita",
    "Phocoena spinipinnis": "Burmeister's Porpoise",
    "Phocoenoides dalli": "Dall's Porpoise",
    "Physeter macrocephalus": "Sperm Whale",
    "Platanista gangetica": "Ganges River Dolphin",
    "Pseudorca crassidens": "False Killer Whale",
    "Stenella attenuata": "Pantropical Spotted Dolphin",
    "Stenella longirostris": "Spinner Dolphin",
    "Tursiops aduncus": "Indo-Pacific Bottlenose Dolphin",
    "Tursiops australis": "Burrunan Dolphin",
    "Tursiops truncatus": "Common Bottlenose Dolphin",
    "Ziphius cavirostris": "Cuvier's Beaked Whale",
    "Megaptera novaeangliae": REFERENCE_SPECIES_COMMON_NAME  # Scientific name for Humpback
}
COMMON_TO_SCIENTIFIC_MAP = {v: k for k, v in SCIENTIFIC_TO_COMMON_MAP.items()}

# Data from your "List of the closest relatives" image (Image 2)
# Keys are scientific names of species, values are scientific names of their closest relatives.
CLOSEST_RELATIVES_SCIENTIFIC = {
    "Balaena mysticetus": "Eubalaena japonica",
    "Balaenoptera musculus": "Balaenoptera physalus",
    "Balaenoptera physalus": "Balaenoptera musculus",
    "Delphinapterus leucas": "Peponocephala electra",  # As per table
    "Delphinus delphis": "Tursiops aduncus",
    "Eubalaena australis": "Eubalaena glacialis",
    "Eubalaena glacialis": "Eubalaena australis",
    "Eubalaena japonica": "Eubalaena australis / E. glacialis",  # Will take the first or handle as needed
    "Globicephala macrorhynchus": "Pseudorca crassidens",
    "Hyperoodon ampullatus": "Ziphius cavirostris",
    "Mesoplodon densirostris": "Mesoplodon grayi",
    "Mesoplodon europaeus": "Mesoplodon mirus",
    "Mesoplodon grayi": "Mesoplodon densirostris",
    "Mesoplodon mirus": "Mesoplodon europaeus",
    "Monodon monoceros": "Delphinapterus leucas",
    "Neophocaena asiaeorientalis": "Neophocaena phocaenoides",
    "Neophocaena phocaenoides": "Neophocaena asiaeorientalis",
    "Orcaella brevirostris": "Orcinus orca",
    "Orcinus orca": "Orcaella brevirostris",
    "Peponocephala electra": "Delphinapterus leucas",
    "Phocoena phocoena": "Phocoenoides dalli",
    "Phocoena sinus": "Phocoena spinipinnis",
    "Phocoena spinipinnis": "Phocoena sinus",
    "Phocoenoides dalli": "Phocoena phocoena",
    "Physeter macrocephalus": "Platanista gangetica",
    "Platanista gangetica": "Physeter macrocephalus",
    "Pseudorca crassidens": "Globicephala macrorhynchus",
    "Stenella attenuata": "Stenella longirostris",
    "Stenella longirostris": "Stenella attenuata",
    "Tursiops aduncus": "Delphinus delphis",  # Note: Image 2 says Tursiops truncatus -> Tursiops australis
    "Tursiops australis": "Tursiops truncatus",
    "Tursiops truncatus": "Tursiops australis",
    "Ziphius cavirostris": "Hyperoodon ampullatus"
}


# --- Helper Functions ---

def parse_species_common_name_from_header(header_description):
    """
    Parses the common species name from a FASTA header.
    Assumes format like "Common Species Name ---Sequence 1" or just "Common Species Name"
    This function is CRITICAL and might need adjustment based on your exact FASTA header format.
    """
    if "---" in header_description:
        # E.g., "Blue whale ---Sequence 1"
        name_part = header_description.split("---")[0].strip()
        return name_part
    else:
        # Try to find a multi-word common name if no "---"
        # This is heuristic: assumes common names are 1-3 words long at the start
        parts = header_description.split()
        for i in range(min(3, len(parts)), 0, -1):  # Check 3-word, then 2-word, then 1-word
            potential_name = " ".join(parts[:i])
            if potential_name in COMMON_TO_SCIENTIFIC_MAP:
                return potential_name
        # Fallback to the first part of the header if it looks like a single word name
        # or if the full name wasn't in our map.
        # This is risky if not all common names are in COMMON_TO_SCIENTIFIC_MAP
        # print(f"Warning: Could not robustly parse common name from '{header_description}'. Using first part: '{parts[0]}'")
        return parts[0].strip()  # Fallback, might be incorrect


def calculate_score(seq1_str, seq2_str):
    """Calculates score by counting matching characters."""
    score = 0
    min_len = min(len(seq1_str), len(seq2_str))
    for i in range(min_len):
        if seq1_str[i] == seq2_str[i]:
            score += 1
    return score


# --- Main Script ---
def main():
    # 1. Read the reference sequence
    try:
        reference_record = SeqIO.read(REFERENCE_FASTA_PATH, "fasta")
        reference_seq_str = str(reference_record.seq).upper()  # Ensure uppercase for comparison
        print(f"Successfully read reference sequence: {reference_record.id} (Length: {len(reference_seq_str)})")
    except Exception as e:
        print(f"Error reading reference FASTA file '{REFERENCE_FASTA_PATH}': {e}")
        return

    # 2. Process comparison sequences
    species_scores = defaultdict(list)  # Stores list of scores for each species

    # Find all FASTA files in the comparison directory
    fasta_files_to_compare = glob.glob(os.path.join(COMPARISON_SEQUENCES_DIR, "*.fasta"))
    if not fasta_files_to_compare:
        print(f"No FASTA files found in directory: {COMPARISON_SEQUENCES_DIR}")
        return

    print(f"Found {len(fasta_files_to_compare)} FASTA files for comparison.")

    for fasta_file_path in fasta_files_to_compare:
        print(f"Processing file: {fasta_file_path}...")
        try:
            for record in SeqIO.parse(fasta_file_path, "fasta"):
                # Extract common species name from FASTA header
                # This is a critical step and assumes FASTA headers contain recognizable common names
                common_name_parsed = parse_species_common_name_from_header(record.description)

                if not common_name_parsed:
                    print(f"  Skipping record '{record.id}' due to parsing failure for species name.")
                    continue

                # Compare with reference and store score
                current_seq_str = str(record.seq).upper()  # Ensure uppercase
                score = calculate_score(reference_seq_str, current_seq_str)
                species_scores[common_name_parsed].append(score)
                # print(f"  Processed: {common_name_parsed} (ID: {record.id}), Score: {score}")

        except Exception as e:
            print(f"  Error processing file '{fasta_file_path}': {e}")
            continue

    if not species_scores:
        print("No scores were calculated. Check FASTA files and parsing logic.")
        return

    # 3. Calculate average scores and prepare results
    results_data = []
    for common_name, scores_list in species_scores.items():
        if not scores_list:
            continue
        average_score = sum(scores_list) / len(scores_list)

        # Get the scientific name for lookup in the closest relatives table
        scientific_name = COMMON_TO_SCIENTIFIC_MAP.get(common_name)

        closest_relative_common_name = "N/A"  # Default
        if scientific_name:
            closest_relative_scientific = CLOSEST_RELATIVES_SCIENTIFIC.get(scientific_name)
            if closest_relative_scientific:
                # Handle cases like "Eubalaena australis / E. glacialis"
                if "/" in closest_relative_scientific:
                    # Taking the first one, or you can define custom logic
                    closest_relative_scientific = closest_relative_scientific.split('/')[0].strip()

                # Convert closest relative back to common name for output
                closest_relative_common_name = SCIENTIFIC_TO_COMMON_MAP.get(closest_relative_scientific,
                                                                            "Unknown (Scientific: " + closest_relative_scientific + ")")
        else:
            print(
                f"Warning: Common name '{common_name}' not found in COMMON_TO_SCIENTIFIC_MAP. Cannot find closest relative.")

        results_data.append({
            "Species": common_name,
            "Average Score": average_score,
            "Closest Relative (from table)": closest_relative_common_name
        })

    if not results_data:
        print("No results to output.")
        return

    # 4. Sort results by average score (descending)
    results_data.sort(key=lambda x: x["Average Score"], reverse=True)

    # 5. Output to Excel
    df_results = pd.DataFrame(results_data)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(OUTPUT_EXCEL_PATH), exist_ok=True)

    try:
        df_results.to_excel(OUTPUT_EXCEL_PATH, index=False)
        print(f"\nResults successfully written to: {OUTPUT_EXCEL_PATH}")
        print("\nPreview of results:")
        print(df_results.head())
    except Exception as e:
        print(f"Error writing Excel file: {e}")


if __name__ == "__main__":
    main()
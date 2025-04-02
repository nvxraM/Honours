import os
import csv
from collections import defaultdict
from Bio import SeqIO

# Define directories
base_dir = "sequences/gb"  # Base directory containing species folders
output_dir = "sequences/zpapers"  # Directory to store the output CSV
os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist
output_file = os.path.join(output_dir, "titles.csv")  # Path for the output file

# Dictionary to store title occurrences across all species
title_count = defaultdict(int)  # Counts occurrences of each unique title
species_mapping = defaultdict(set)  # Maps titles to the species they belong to

total_files = 0  # Counter for total GenBank files processed

# Iterate over each species folder in the base directory
for species in os.listdir(base_dir):
    species_dir = os.path.join(base_dir, species)
    if os.path.isdir(species_dir):  # Ensure the path is a directory
        for file in os.listdir(species_dir):
            # Process only GenBank files (.gb or .gbk extensions)
            if file.lower().endswith((".gb", ".gbk")):
                total_files += 1  # Increment file counter
                file_path = os.path.join(species_dir, file)
                with open(file_path, "r") as handle:
                    try:
                        # Parse the GenBank file using BioPython
                        for record in SeqIO.parse(handle, "genbank"):
                            if "references" in record.annotations:  # Check if references exist
                                for ref in record.annotations["references"]:
                                    # Extract and store title if available
                                    if hasattr(ref, "title") and ref.title:
                                        title_count[ref.title] += 1  # Increment title count
                                        species_mapping[ref.title].add(species)  # Add species to mapping
                    except Exception as e:
                        print(f"Error processing {file_path}: {e}")  # Print error if file cannot be processed

# Write collected titles and occurrences to a CSV file
with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Title", "Occurrences", "Species"])  # Write header row
    for title, count in title_count.items():
        writer.writerow([title, count, "; ".join(sorted(species_mapping[title]))])  # Write title, occurrence count, and associated species

print(f"Titles extracted and saved to {output_file}")
print(f"Total GenBank files processed: {total_files}")

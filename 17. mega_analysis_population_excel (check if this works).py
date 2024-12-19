import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Base directory now points to the new structure
base_dir = 'sequences/Interactive_group_editor'

# Regex pattern to parse lines.
# Assuming lines are in the format:
# SpeciesName(GenusName)  value  SE
pattern = re.compile(r'^(.+?)\((.*?)\)\s+(\S+)\s+(\S+)$')

# Data structure to hold results:
# results[Genus][Species] = { 'd':..., 'd_SE':..., 'dN':..., 'dN_SE':..., 'dS':..., 'dS_SE':... }
results = {}

def parse_file(filepath, genus, value_key, se_key):
    """
    Parse a given MEGA result file containing lines like:
    SpeciesName(GenusName)  value  SE
    Store the values in the results dictionary.
    """
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                m = pattern.match(line)
                if m:
                    species_name, genus_name, value_str, se_str = m.groups()
                    if genus not in results:
                        results[genus] = {}
                    if species_name not in results[genus]:
                        results[genus][species_name] = {}
                    results[genus][species_name][value_key] = float(value_str)
                    results[genus][species_name][se_key] = float(se_str)

# Scan the directory structure:
# sequences/Interactive_group_editor/[Genus]/[Species]/MEGA_Populations/
for genus in os.listdir(base_dir):
    genus_path = os.path.join(base_dir, genus)
    if not os.path.isdir(genus_path):
        continue

    for species in os.listdir(genus_path):
        species_path = os.path.join(genus_path, species)
        if not os.path.isdir(species_path):
            continue

        mega_pop_path = os.path.join(species_path, "MEGA_Populations")
        if not os.path.isdir(mega_pop_path):
            continue

        # Files we expect:
        # [Species]_d_value.meg
        # [Species]_nonsynonymous_only.meg
        # [Species]_Synonymous_only.meg
        d_file = os.path.join(mega_pop_path, f'{species}_d_value.meg')
        nonsyn_file = os.path.join(mega_pop_path, f'{species}_nonsynonymous_only.meg')
        syn_file = os.path.join(mega_pop_path, f'{species}_Synonymous_only.meg')

        # Parse the files if they exist
        parse_file(d_file, genus, 'd', 'd_SE')
        parse_file(nonsyn_file, genus, 'dN', 'dN_SE')
        parse_file(syn_file, genus, 'dS', 'dS_SE')

# Compute the average GC content and sequence count per species.
# FASTA located at: sequences/Interactive_group_editor/[Genus]/[Species]/[Species]_concatenated_gapped_sequences.fasta
avg_gc_per_species = {}
sequence_count_per_species = {}

for genus in results.keys():
    avg_gc_per_species[genus] = {}
    sequence_count_per_species[genus] = {}

    for species in results[genus].keys():
        # Convert species to underscore form in case needed
        # The code currently expects species_name exactly as in results
        # If the FASTA headers contain underscores, we must ensure consistency.
        # The parsing pattern from the MEGA files should produce a species_name that matches the FASTA headers.
        # If they differ, adjust accordingly. Here we assume the species_name matches headers that end in Genus_species
        species_underscore = species.replace(' ', '_')

        fasta_path = os.path.join(base_dir, genus, species_underscore, f"{species_underscore}_concatenated_gapped_sequences.fasta")
        if os.path.exists(fasta_path):
            species_sequences = []
            for record in SeqIO.parse(fasta_path, "fasta"):
                species_sequences.append(str(record.seq))

            if species_sequences:
                full_seq = ''.join(species_sequences)
                gc_val = gc_fraction(full_seq) * 100
                seq_count = len(species_sequences)
                avg_gc_per_species[genus][species] = gc_val
                sequence_count_per_species[genus][species] = seq_count
            else:
                avg_gc_per_species[genus][species] = None
                sequence_count_per_species[genus][species] = None
        else:
            avg_gc_per_species[genus][species] = None
            sequence_count_per_species[genus][species] = None

# Build a list of rows for the DataFrame
rows = []
for genus, species_data in results.items():
    for species, metrics in species_data.items():
        d = metrics.get('d')
        d_SE = metrics.get('d_SE')
        dN = metrics.get('dN')
        dN_SE = metrics.get('dN_SE')
        dS = metrics.get('dS')
        dS_SE = metrics.get('dS_SE')

        # dN/dS Ratio as blank (or compute if you have both dN and dS)
        dn_ds_ratio = None

        gc_value = avg_gc_per_species.get(genus, {}).get(species, None)
        seq_count = sequence_count_per_species.get(genus, {}).get(species, None)

        # Blank columns for additional data
        avg_gen_time = None
        avg_body_size = None
        avg_life_span = None

        row = {
            'Genus': genus,
            'Species': species,
            'd': d,
            'd S.E.': d_SE,
            'dN': dN,
            'dN S.E.': dN_SE,
            'dS': dS,
            'dS S.E.': dS_SE,
            'dN/dS Ratio': dn_ds_ratio,
            'Avg genome GC (%)': gc_value,
            'Sequence No.': seq_count,
            'Avg. generation time (Years)': avg_gen_time,
            'Avg. body Size (Kg)': avg_body_size,
            'Avg. Life Span (Years)': avg_life_span
        }
        rows.append(row)

# Create DataFrame with the desired columns
df = pd.DataFrame(rows, columns=[
    'Genus', 'Species', 'd', 'd S.E.', 'dN', 'dN S.E.',
    'dS', 'dS S.E.', 'dN/dS Ratio', 'Avg genome GC (%)',
    'Sequence No.', 'Avg. generation time (Years)',
    'Avg. body Size (Kg)', 'Avg. Life Span (Years)'
])

# Ensure output directory exists and write to Excel
output_path = os.path.join(base_dir, 'All_Populations_Analysis.xlsx')
os.makedirs(os.path.dirname(output_path), exist_ok=True)

df.to_excel(output_path, index=False)
print(f"Data successfully written to {output_path}")

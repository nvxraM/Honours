import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Root directory containing sequences
root_dir = 'sequences/CDS_Genus'

# Regex pattern to parse lines.
# This pattern:
# - Captures species name (including spaces) until '('
# - Captures genus name inside parentheses
# - Captures two numerical values (value and SE)
pattern = re.compile(r'^(.+?)\((.*?)\)\s+(\S+)\s+(\S+)$')

# Data structure: results[Genus][Species] = { 'd':..., 'd_SE':..., 'dN':..., 'dN_SE':..., 'dS':..., 'dS_SE':... }
results = {}

def parse_file(filepath, genus, value_key, se_key):
    """
    Parse a given file containing lines like:
    SpeciesName(GenusName)  value  SE
    and store the values in results.
    """
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                m = pattern.match(line)
                if m:
                    species_name, genus_name, value_str, se_str = m.groups()
                    if species_name not in results[genus]:
                        results[genus][species_name] = {}
                    results[genus][species_name][value_key] = float(value_str)
                    results[genus][species_name][se_key] = float(se_str)

# Parse all text files for each genus
for genus in os.listdir(root_dir):
    genus_path = os.path.join(root_dir, genus, 'MEGA_Grouped_by_Species')
    if not os.path.isdir(genus_path):
        continue

    d_file = os.path.join(genus_path, f'{genus}_d_value.txt')
    nonsyn_file = os.path.join(genus_path, f'{genus}_nonsynonymous_only.txt')
    syn_file = os.path.join(genus_path, f'{genus}_Synonymous_only.txt')

    if genus not in results:
        results[genus] = {}

    parse_file(d_file, genus, 'd', 'd_SE')
    parse_file(nonsyn_file, genus, 'dN', 'dN_SE')
    parse_file(syn_file, genus, 'dS', 'dS_SE')

# Compute the average GC content and sequence count per species
# We have a single concatenated FASTA for each genus:
# sequences/CDS_Genus/<Genus>/CDS_nucleotide_gapped/<Genus>_concatenated_gapped_sequences.fasta
avg_gc_per_species = {}
sequence_count_per_species = {}

for genus in results.keys():
    # Dictionary to hold sequences for each species in this genus
    species_sequences = {}
    fasta_path = os.path.join(root_dir, genus, 'CDS_nucleotide_gapped', f'{genus}_concatenated_gapped_sequences.fasta')

    if os.path.exists(fasta_path):
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Header format assumed: >Accession_SpeciesName
            # Example: >KC572788.1_Balaenoptera_physalus
            header = record.id
            parts = header.split('_', 1)
            if len(parts) == 2:
                # The second part should be the species name with underscores
                species_in_fasta = parts[1]
                if species_in_fasta not in species_sequences:
                    species_sequences[species_in_fasta] = []
                species_sequences[species_in_fasta].append(str(record.seq))

        # Now compute GC and sequence count for each species
        genus_gc = {}
        genus_count = {}
        for sp_name, seq_list in species_sequences.items():
            full_seq = ''.join(seq_list)
            genus_gc[sp_name] = gc_fraction(full_seq) * 100
            genus_count[sp_name] = len(seq_list)

        avg_gc_per_species[genus] = genus_gc
        sequence_count_per_species[genus] = genus_count
    else:
        # No FASTA file found for this genus
        avg_gc_per_species[genus] = {}
        sequence_count_per_species[genus] = {}

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

        # dN/dS Ratio as blank
        dn_ds_ratio = None

        # Convert species to underscore form to look up in avg_gc_per_species and sequence counts
        species_underscore = species.replace(' ', '_')
        gc_value = avg_gc_per_species.get(genus, {}).get(species_underscore, None)
        seq_count = sequence_count_per_species.get(genus, {}).get(species_underscore, None)

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

# Ensure output directory exists
output_path = 'sequences/MEGA/All_Grouped_by_Species.xlsx'
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Write to Excel
df.to_excel(output_path, index=False)
print(f"Data successfully written to {output_path}")

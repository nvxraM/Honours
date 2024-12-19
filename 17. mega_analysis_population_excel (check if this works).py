import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

base_dir = 'sequences/Interactive_group_editor'

# Matches lines like: "Balaena mysticetus Population 1    0.00171694  0.00029752"
pattern = re.compile(r'^(\S+)\s+(\S+)\s+Population\s+(\d+)\s+(\S+)\s+(\S+)$')

results = {}

def parse_file(filepath, genus, species, value_key, se_key):
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                m = pattern.match(line)
                if m:
                    line_genus, line_species, population_str, value_str, se_str = m.groups()

                    pop_num = int(population_str)

                    if genus not in results:
                        results[genus] = {}
                    if species not in results[genus]:
                        results[genus][species] = {}

                    if pop_num not in results[genus][species]:
                        results[genus][species][pop_num] = {}

                    results[genus][species][pop_num][value_key] = float(value_str)
                    results[genus][species][pop_num][se_key] = float(se_str)

def gc_content_third_position(seq_list):
    third_position_bases = []
    for seq_str in seq_list:
        codon_third_bases = seq_str[2::3]
        third_position_bases.append(codon_third_bases)

    all_third_positions = "".join(third_position_bases)
    if all_third_positions:
        return gc_fraction(all_third_positions) * 100
    return None

# Scan directories
for genus in os.listdir(base_dir):
    genus_path = os.path.join(base_dir, genus)
    if not os.path.isdir(genus_path):
        continue

    for species_dir_name in os.listdir(genus_path):
        species_path = os.path.join(genus_path, species_dir_name)
        if not os.path.isdir(species_path):
            continue

        mega_pop_path = os.path.join(species_path, "MEGA_Populations")
        if not os.path.isdir(mega_pop_path):
            continue

        d_file = os.path.join(mega_pop_path, f'{species_dir_name}_d_value.txt')
        nonsyn_file = os.path.join(mega_pop_path, f'{species_dir_name}_nonsynonymous_only.txt')
        syn_file = os.path.join(mega_pop_path, f'{species_dir_name}_Synonymous_only.txt')

        if os.path.exists(d_file):
            parse_file(d_file, genus, species_dir_name, 'd', 'd_SE')
        if os.path.exists(nonsyn_file):
            parse_file(nonsyn_file, genus, species_dir_name, 'dN', 'dN_SE')
        if os.path.exists(syn_file):
            parse_file(syn_file, genus, species_dir_name, 'dS', 'dS_SE')

avg_gc_per_species = {}
sequence_count_per_species = {}

for genus in results.keys():
    avg_gc_per_species[genus] = {}
    sequence_count_per_species[genus] = {}

    for species in results[genus].keys():
        fasta_path = os.path.join(base_dir, genus, species, f"{species}_concatenated_gapped_sequences.fasta")
        if os.path.exists(fasta_path):
            species_sequences = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]
            if species_sequences:
                gc_val = gc_content_third_position(species_sequences)
                seq_count = len(species_sequences)
                avg_gc_per_species[genus][species] = gc_val
                sequence_count_per_species[genus][species] = seq_count
            else:
                avg_gc_per_species[genus][species] = None
                sequence_count_per_species[genus][species] = None
        else:
            avg_gc_per_species[genus][species] = None
            sequence_count_per_species[genus][species] = None

rows = []
for genus, species_data in results.items():
    for species, pop_data in species_data.items():
        for population, metrics in pop_data.items():
            d = metrics.get('d')
            d_SE = metrics.get('d_SE')
            dN = metrics.get('dN')
            dN_SE = metrics.get('dN_SE')
            dS = metrics.get('dS')
            dS_SE = metrics.get('dS_SE')

            dn_ds_ratio = None
            gc_value = avg_gc_per_species.get(genus, {}).get(species, None)
            seq_count = sequence_count_per_species.get(genus, {}).get(species, None)

            avg_gen_time = None
            avg_body_size = None
            avg_life_span = None

            row = {
                'Genus': genus,
                'Species': species,
                'Population': population,
                'd': d,
                'd S.E.': d_SE,
                'dN': dN,
                'dN S.E.': dN_SE,
                'dS': dS,
                'dS S.E.': dS_SE,
                'dN/dS Ratio': dn_ds_ratio,
                'GC3 (%)': gc_value,  # Renamed column here
                'Sequence No.': seq_count,
                'Avg. generation time (Years)': avg_gen_time,
                'Avg. body Size (Kg)': avg_body_size,
                'Avg. Life Span (Years)': avg_life_span
            }
            rows.append(row)

df = pd.DataFrame(rows, columns=[
    'Genus', 'Species', 'Population', 'd', 'd S.E.', 'dN', 'dN S.E.',
    'dS', 'dS S.E.', 'dN/dS Ratio', 'GC3 (%)',
    'Sequence No.', 'Avg. generation time (Years)',
    'Avg. body Size (Kg)', 'Avg. Life Span (Years)'
])

output_path = os.path.join(base_dir, 'All_Populations_Analysis.xlsx')
os.makedirs(os.path.dirname(output_path), exist_ok=True)

df.to_excel(output_path, index=False)
print(f"Data successfully written to {output_path}")

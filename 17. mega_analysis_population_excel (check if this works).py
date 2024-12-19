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
        # Extract every 3rd base starting from the third base (index 2)
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

        # Parse metric files (d, dN, dS)
        if os.path.exists(d_file):
            parse_file(d_file, genus, species_dir_name, 'd', 'd_SE')
        if os.path.exists(nonsyn_file):
            parse_file(nonsyn_file, genus, species_dir_name, 'dN', 'dN_SE')
        if os.path.exists(syn_file):
            parse_file(syn_file, genus, species_dir_name, 'dS', 'dS_SE')

        # Parse the grp file to map each accession to its population
        grp_file = os.path.join(mega_pop_path, f'{species_dir_name}_population.grp')
        population_to_accessions = {}
        if os.path.exists(grp_file):
            with open(grp_file, 'r', encoding='utf-8') as gf:
                for line in gf:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    # Example line:
                    # PP761291.1_Delphinus_delphis=Delphinus_delphis_Population_2
                    # Left side is accession, right side contains population info
                    if '=' in line:
                        left, right = line.split('=', 1)
                        # right might look like "Delphinus_delphis_Population_2"
                        # Extract population number from the right part
                        pop_match = re.search(r'_Population_(\d+)$', right)
                        if pop_match:
                            pop_num = int(pop_match.group(1))
                            if pop_num not in population_to_accessions:
                                population_to_accessions[pop_num] = []
                            population_to_accessions[pop_num].append(left)

        # Now read the fasta file and assign sequences to populations
        fasta_path = os.path.join(base_dir, genus, species_dir_name,
                                  f"{species_dir_name}_concatenated_gapped_sequences.fasta")
        if os.path.exists(fasta_path) and os.path.getsize(fasta_path) > 0:
            # Create a dictionary to hold sequences per population
            population_sequences = {pop_num: [] for pop_num in population_to_accessions.keys()}

            for record in SeqIO.parse(fasta_path, "fasta"):
                # record.id should match something like 'PP761291.1_Delphinus_delphis'
                # Check which population this record belongs to
                for pop_num, acc_list in population_to_accessions.items():
                    if record.id in acc_list:
                        population_sequences[pop_num].append(str(record.seq))

            # Compute GC content at 3rd codon position for each population
            for pop_num, seqs in population_sequences.items():
                if seqs:
                    gc_val = gc_content_third_position(seqs)
                else:
                    gc_val = None
                # Store the GC3 value in results if that population already exists
                # (It should, since we parsed d/dN/dS files, but if not, create it)
                if genus not in results:
                    results[genus] = {}
                if species_dir_name not in results[genus]:
                    results[genus][species_dir_name] = {}
                if pop_num not in results[genus][species_dir_name]:
                    results[genus][species_dir_name][pop_num] = {}
                results[genus][species_dir_name][pop_num]['GC3'] = gc_val
                results[genus][species_dir_name][pop_num]['Sequence_Count'] = len(seqs)
        else:
            # If no fasta found or it's empty, still ensure GC3 entries exist
            # This is only relevant if there were populations defined
            for pop_num in population_to_accessions.keys():
                if genus not in results:
                    results[genus] = {}
                if species_dir_name not in results[genus]:
                    results[genus][species_dir_name] = {}
                if pop_num not in results[genus][species_dir_name]:
                    results[genus][species_dir_name][pop_num] = {}
                # No sequences, so no GC3
                results[genus][species_dir_name][pop_num]['GC3'] = None
                results[genus][species_dir_name][pop_num]['Sequence_Count'] = None

# Now generate the DataFrame
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
            gc_value = metrics.get('GC3')
            seq_count = metrics.get('Sequence_Count')

            dn_ds_ratio = None
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
                'Avg GC3 (%)': gc_value,
                'Sequence No.': seq_count,
                'Avg. generation time (Years)': avg_gen_time,
                'Avg. body Size (Kg)': avg_body_size,
                'Avg. Life Span (Years)': avg_life_span
            }
            rows.append(row)

df = pd.DataFrame(rows, columns=[
    'Genus', 'Species', 'Population', 'd', 'd S.E.', 'dN', 'dN S.E.',
    'dS', 'dS S.E.', 'dN/dS Ratio', 'Avg GC3 (%)',
    'Sequence No.', 'Avg. generation time (Years)',
    'Avg. body Size (Kg)', 'Avg. Life Span (Years)'
])

output_path = os.path.join(base_dir, 'All_Populations_Analysis.xlsx')
os.makedirs(os.path.dirname(output_path), exist_ok=True)

df.to_excel(output_path, index=False)
print(f"Data successfully written to {output_path}")

import os
from pathlib import Path

STANDARD_BASES = set("ATCG-")
DIGIT_TO_WORD = {
    '1': 'one', '2': 'two', '3': 'three', '4': 'four', '5': 'five',
    '6': 'six', '7': 'seven', '8': 'eight', '9': 'nine', '10': 'ten'
}

def replace_ambiguous_bases(seq):
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def convert_pop_format(pop, species):
    # Only convert pop names, don't touch sample names
    if pop.startswith(species + "_Population_"):
        num = pop.replace(species + "_Population_", "")
        word = DIGIT_TO_WORD.get(num, num)
        return f"{species}_Population_{word}"
    else:
        return f"{species}_{pop}"

def read_grp(grp_path, species):
    mapping = {}
    with open(grp_path) as f:
        for line in f:
            if '=' in line:
                sample, pop = line.strip().split('=', 1)
                sample_key = sample.strip()
                pop_value = convert_pop_format(pop.strip(), species)
                mapping[sample_key] = pop_value
    return mapping

def read_fasta(fasta_path):
    seqs = []
    with open(fasta_path) as f:
        lines = [l.strip() for l in f]
    name, seq = None, ""
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, replace_ambiguous_bases(seq)))
            name = line[1:].strip()
            seq = ""
        else:
            seq += line
    if name is not None:
        seqs.append((name, replace_ambiguous_bases(seq)))
    return seqs

def write_mega(out_path, seqs, title="MEGA Alignment"):
    with open(out_path, 'w') as f:
        f.write("#MEGA\n")
        f.write(f"!Title {title};\n")
        f.write("!Format DataType=DNA Haploid GeneticCode=2;\n\n")
        for name, seq in seqs:
            f.write(f"#{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_sg_file(out_path, seqs, mapping):
    with open(out_path, 'w') as f:
        for name, _ in seqs:
            pop = mapping.get(name)
            if pop:
                f.write(f"{name}\t{pop}\n")
            else:
                print(f"Warning: {name} not found in .grp mapping!")

def write_alleles_file(out_path, seqs):
    with open(out_path, 'w') as f:
        for name, seq in seqs:
            f.write(f">{name} {seq}\n")
        f.write("//\n")  # Append // at the end

def process_all_species_haploid(base_path, output_base):
    base = Path(base_path)
    output_base = Path(output_base)
    for genus_folder in base.iterdir():
        if not genus_folder.is_dir():
            continue
        for species_folder in genus_folder.iterdir():
            if not species_folder.is_dir():
                continue
            species_name = species_folder.name
            fasta_path = species_folder / f"{species_name}_concatenated_gapped_sequences.fasta"
            grp_path = species_folder / "MEGA_Populations" / f"{species_name}_population.grp"
            if not fasta_path.exists() or not grp_path.exists():
                print(f"Skipping {species_folder} (missing FASTA or GRP file)")
                continue
            mapping = read_grp(grp_path, species_name)
            seqs = read_fasta(fasta_path)
            pops = set(mapping.get(name) for name, _ in seqs if mapping.get(name) is not None)
            species_out = output_base / species_folder.name
            species_out.mkdir(parents=True, exist_ok=True)
            # Create the Results directory (empty)
            (species_out / "Results").mkdir(parents=True, exist_ok=True)
            sg_file = species_out / f"{species_folder.name}.SG.txt"
            write_sg_file(sg_file, seqs, mapping)
            for pop in pops:
                pop_seqs = [(name, seq) for name, seq in seqs if mapping.get(name, None) == pop]
                if pop_seqs:
                    outname = f"{species_folder.name}_{pop}.meg"
                    write_mega(species_out / outname, pop_seqs, title=f"{species_folder.name} {pop} Haploid Alignment")
            alleles_path = species_out / f"{species_folder.name}.alleles"
            write_alleles_file(alleles_path, seqs)
    print("Wrote MEGA, SG, and alleles files for all populations/species.")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Mega_Haploid"

if __name__ == "__main__":
    Path(output_base).mkdir(parents=True, exist_ok=True)
    process_all_species_haploid(input_base, output_base)

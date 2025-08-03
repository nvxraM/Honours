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

def write_nexus_with_traits(out_path, seqs, mapping, species_name):
    ntax = len(seqs)
    nchar = len(seqs[0][1]) if seqs else 0
    with open(out_path, 'w') as f:
        f.write("#NEXUS\n\n")
        f.write("Begin data;\n")
        f.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        f.write("    Format datatype=dna symbols=\"ACTGN-\" missing=? gap=-;\n")
        f.write("    Matrix\n")
        for name, seq in seqs:
            f.write(f"{name.replace(' ', '_'):30s} {seq}\n")
        f.write("    ;\nEnd;\n\n")
        # Add population trait block
        f.write("Begin traits;\n")
        f.write("    Dimensions NTraits=1;\n")
        f.write("    Format labels=yes missing=? separator=Comma;\n")
        f.write("    TraitLabels population;\n")
        f.write("    Matrix\n")
        for name, _ in seqs:
            pop = mapping.get(name, "?")
            f.write(f"    {name.replace(' ', '_'):30s} {pop}\n")
        f.write("    ;\nEnd;\n")

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
            species_out = output_base / species_folder.name
            species_out.mkdir(parents=True, exist_ok=True)
            # Write NEXUS file with traits block for populations
            nexus_out = species_out / f"{species_folder.name}.nexus"
            write_nexus_with_traits(nexus_out, seqs, mapping, species_name)
    print("Wrote NEXUS files with population info for all species.")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Nexus"

if __name__ == "__main__":
    Path(output_base).mkdir(parents=True, exist_ok=True)
    process_all_species_haploid(input_base, output_base)

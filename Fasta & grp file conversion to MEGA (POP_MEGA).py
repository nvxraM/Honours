import os
from pathlib import Path



#the fix for the file not getting accepted when running a multi on dnaSP is to change the SG file so its like the following
#KY026766.1_Balaena_mysticetus_0.	Balaena_mysticetus_0._Population_one
#and the other alleles file should looke like
#>MN145937.1_Balaena_mysticetus_0.
#weird but it worked, only reason I didnt fix this because the multi didnt get the standard deviations....FML

STANDARD_BASES = set("ATCG-")

DIGIT_TO_WORD = {
    '1': '_one', '2': '_two', '3': '_three', '4': '_four', '5': '_five',
    '6': '_six', '7': '_seven', '8': '_eight', '9': '_nine', '10': '_ten'
    # Extend as needed!
}

def replace_ambiguous_bases(seq):
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def popnum_to_word(popname):
    for num, word in DIGIT_TO_WORD.items():
        if popname.endswith('_' + num):
            return popname[:-(len(num)+1)] + word
    return popname

def read_grp(grp_path):
    mapping = {}
    with open(grp_path) as f:
        for line in f:
            if '=' in line:
                sample, pop = line.strip().split('=', 1)
                mapping[sample.strip()] = pop.strip()
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

def process_all_species_diploid(base_path, diploid_base):
    base = Path(base_path)
    diploid_base = Path(diploid_base)
    for genus_folder in base.iterdir():
        if not genus_folder.is_dir():
            continue
        for species_folder in genus_folder.iterdir():
            if not species_folder.is_dir():
                continue
            fasta_path = species_folder / f"{species_folder.name}_concatenated_gapped_sequences.fasta"
            grp_path = species_folder / "MEGA_Populations" / f"{species_folder.name}_population.grp"
            if not fasta_path.exists() or not grp_path.exists():
                print(f"Skipping {species_folder} (missing FASTA or GRP file)")
                continue
            mapping = read_grp(grp_path)
            seqs = read_fasta(fasta_path)
            # Get all unique populations
            pops = set(mapping[sample] for sample, _ in seqs if sample in mapping)
            pop_folder = diploid_base / species_folder.name / "populations"
            pop_folder.mkdir(parents=True, exist_ok=True)
            for pop in pops:
                pop_word = popnum_to_word(pop)
                pop_seqs = [(name, seq) for name, seq in seqs if mapping.get(name, None) == pop]
                if pop_seqs:
                    outname = f"{species_folder.name}_{pop_word}.meg"
                    write_mega(pop_folder / outname, pop_seqs, title=f"{species_folder.name} {pop_word} Diploid Alignment")
    print("Wrote diploid MEGA files with each population in its own file.")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
diploid_base = "sequences/Species_POP_Mega_Diploid"

if __name__ == "__main__":
    Path(diploid_base).mkdir(parents=True, exist_ok=True)
    process_all_species_diploid(input_base, diploid_base)

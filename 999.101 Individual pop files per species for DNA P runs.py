import os
from pathlib import Path

# Set of standard nucleotides and gap
STANDARD_BASES = set("ATCG-")

def replace_ambiguous_bases(seq):
    """Replace ambiguous nucleotides with N."""
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def read_grp(grp_path):
    """Read group (population) mapping file."""
    mapping = {}
    with open(grp_path) as f:
        for line in f:
            if '=' in line:
                sample, pop = line.strip().split('=', 1)
                mapping[sample.strip()] = pop.strip()
    return mapping

def read_fasta(fasta_path):
    """Read a FASTA file and return a list of (name, cleaned sequence) tuples."""
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
    """Write MEGA format file, wrapping sequences at 60 characters."""
    with open(out_path, 'w') as f:
        f.write("#MEGA\n")
        f.write(f"!Title {title};\n")
        f.write("!Description ChromosomalLocation=Mitochondrial;\n")  # <-- Mitochondrial annotation
        f.write("!Format DataType=DNA Haploid GeneticCode=2;\n\n")
        for name, seq in seqs:
            f.write(f"#{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def process_all_species_haploid(base_path, haploid_base):
    """
    For each species, writes a separate MEGA (.meg) file for each population,
    directly inside the species folder (no subfolder).
    """
    base = Path(base_path)
    haploid_base = Path(haploid_base)
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
            pops = set(mapping[sample] for sample, _ in seqs if sample in mapping)
            species_out = haploid_base / species_folder.name  # Output folder = species folder
            species_out.mkdir(parents=True, exist_ok=True)
            for pop in pops:
                pop_seqs = [(name, seq) for name, seq in seqs if mapping.get(name, None) == pop]
                if pop_seqs:
                    outname = f"{species_folder.name}_{pop}.meg"
                    write_mega(species_out / outname, pop_seqs, title=f"{species_folder.name} {pop} Haploid Mitochondrial Alignment")
    print("Wrote haploid MEGA files (mitochondrial) with each population in its own file.")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
haploid_base = "sequences/Species_POP_Mega_Haploid_split_split"

if __name__ == "__main__":
    Path(haploid_base).mkdir(parents=True, exist_ok=True)
    process_all_species_haploid(input_base, haploid_base)

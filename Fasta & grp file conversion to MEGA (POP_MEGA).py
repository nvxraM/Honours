import os
from pathlib import Path

STANDARD_BASES = set("ATCG-")

def get_pop_suffix(pop):
    # Unused now, but kept if you want to reference population info elsewhere
    for s in pop.split('_'):
        if s.isdigit():
            return s
    return pop.split('_')[-1]

def read_grp(grp_path):
    mapping = {}
    with open(grp_path) as f:
        for line in f:
            if '=' in line:
                sample, pop = line.strip().split('=', 1)
                mapping[sample.strip()] = get_pop_suffix(pop.strip())
    return mapping

def replace_ambiguous_bases(seq):
    # Replace any non-ATCG- with 'N'
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def read_fasta(fasta_path, mapping):
    seqs = []
    with open(fasta_path) as f:
        lines = [l.strip() for l in f]
    name, seq = None, ""
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, replace_ambiguous_bases(seq)))
            orig_name = line[1:].strip()
            # --- Only use orig_name, do not add suffix ---
            name = orig_name
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
        f.write("!Format DataType=DNA Haploid GeneticCode=2;\n\n")  # mtDNA code, haploid
        for name, seq in seqs:
            f.write(f"#{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_fasta(out_path, seqs):
    """Write output in standard FASTA format"""
    with open(out_path, 'w') as f:
        for name, seq in seqs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def write_sg(out_path, fasta_path, mapping):
    """Writes SG file with original sample names and population numbers (space separated!)"""
    with open(fasta_path) as f:
        lines = [l.strip() for l in f if l.startswith(">")]
    with open(out_path, "w") as sg:
        for line in lines:
            sample = line[1:].strip()
            pop = mapping.get(sample, "1")  # Default to 1 if missing
            sg.write(f"{sample} {pop}\n")   # single space!

def process_all_species(base_path, output_path):
    base = Path(base_path)
    output = Path(output_path)
    output.mkdir(parents=True, exist_ok=True)
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
            seqs = read_fasta(fasta_path, mapping)
            # MEGA
            out_mega = output / f"{species_folder.name}.meg"
            write_mega(out_mega, seqs, title=species_folder.name)
            # FASTA (.fa and .fas)
            out_fa = output / f"{species_folder.name}.fa"
            out_fas = output / f"{species_folder.name}.fas"
            write_fasta(out_fa, seqs)
            write_fasta(out_fas, seqs)
            # SG
            out_sg = output / f"{species_folder.name}.SG.txt"
            write_sg(out_sg, fasta_path, mapping)
            print(f"Wrote {out_mega}, {out_fa}, {out_fas}, and {out_sg}")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Mega"

# ------------- RUN THE SCRIPT ----------
if __name__ == "__main__":
    process_all_species(input_base, output_base)

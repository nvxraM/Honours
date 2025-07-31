import os
from pathlib import Path

def get_pop_suffix(pop):
    """Extract the population number from the population string."""
    # E.g., Balaena_mysticetus_Population_1 -> 1
    for s in pop.split('_'):
        if s.isdigit():
            return s
    # If not found, just use the last part after underscore
    return pop.split('_')[-1]

def read_grp(grp_path):
    # Returns dict: {sample: popnum}
    mapping = {}
    with open(grp_path) as f:
        for line in f:
            if '=' in line:
                sample, pop = line.strip().split('=', 1)
                mapping[sample.strip()] = get_pop_suffix(pop.strip())
    return mapping

def read_fasta(fasta_path, mapping):
    # Returns list of tuples: (header, sequence)
    seqs = []
    with open(fasta_path) as f:
        lines = [l.strip() for l in f]
    name, seq = None, ""
    for line in lines:
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, seq))
            orig_name = line[1:].strip()
            suffix = mapping.get(orig_name, "")
            name = f"{orig_name}_{suffix}" if suffix else orig_name
            seq = ""
        else:
            seq += line
    if name is not None:
        seqs.append((name, seq))
    return seqs

def write_nexus(out_path, seqs):
    ntax = len(seqs)
    nchar = len(seqs[0][1]) if seqs else 0
    with open(out_path, 'w') as f:
        f.write("#NEXUS\n\n")
        f.write("BEGIN DATA;\n")
        f.write(f"    DIMENSIONS NTAX={ntax} NCHAR={nchar};\n")
        f.write("    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        f.write("    MATRIX\n")
        for name, seq in seqs:
            f.write(f"{name}\t{seq}\n")
        f.write("    ;\nEND;\n")

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
            out_nexus = output / f"{species_folder.name}.nex"
            write_nexus(out_nexus, seqs)
            print(f"Wrote {out_nexus}")

# -------- CONFIGURE THESE PATHS --------
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Nexus"

# ------------- RUN THE SCRIPT ----------
if __name__ == "__main__":
    process_all_species(input_base, output_base)

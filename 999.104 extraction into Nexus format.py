import os
import re
from pathlib import Path

STANDARD_BASES = set("ATCG-")

WORD_TO_NUM = {
    'one': '1', 'two': '2', 'three': '3', 'four': '4', 'five': '5',
    'six': '6', 'seven': '7', 'eight': '8', 'nine': '9', 'ten': '10'
}

def replace_ambiguous_bases(seq: str) -> str:
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def gp_from(pop: str) -> str:
    """
    Extracts the population number (supports 'Population_1' or 'Population_one'),
    regardless of any species prefix/suffix. Returns 'Gp#'. If no match, returns
    a safe token from the input (spaces -> underscores).
    """
    m = re.search(r'Population_(one|two|three|four|five|six|seven|eight|nine|ten|\d+)$', pop, re.IGNORECASE)
    if not m:
        return pop.replace(" ", "_")
    tok = m.group(1)
    if tok.isdigit():
        num = tok
    else:
        num = WORD_TO_NUM.get(tok.lower(), tok)
    return f"Gp{num}"

def read_grp(grp_path: Path) -> dict:
    """
    Read MEGA .grp file: lines like 'SampleID = Population_1'
    Returns dict: sample -> raw population token (e.g., 'Population_1')
    """
    mapping = {}
    with open(grp_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            sample, pop = line.split('=', 1)
            mapping[sample.strip()] = pop.strip()
    return mapping

def read_fasta(fasta_path: Path):
    """
    Read FASTA; returns list of (name, seq) with ambiguous bases replaced by N.
    """
    seqs = []
    with open(fasta_path, encoding='utf-8') as f:
        name, seq = None, []
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, replace_ambiguous_bases(''.join(seq))))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if name is not None:
            seqs.append((name, replace_ambiguous_bases(''.join(seq))))
    return seqs

def write_nexus(out_path: Path, seqs, label: str = ""):
    """
    Write a NEXUS with MEGA/BEAST-friendly DATA block.
    Optional label is written as a comment (e.g., 'Gp1').
    """
    ntax = len(seqs)
    nchar = len(seqs[0][1]) if seqs else 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("#NEXUS\n\n")
        f.write("Begin DATA;\n")
        f.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        # MEGA expects symbols without '-' here; gap is set separately
        f.write('    Format datatype=dna symbols="ACTGN" missing=? gap=-;\n')
        f.write("Matrix\n")
        for name, seq in seqs:
            f.write(f"{name.replace(' ', '_'):30s} {seq}\n")
        f.write(";\nEnd;\n")
        if label:
            f.write(f"\n[Group: {label}]\n")  # harmless comment for MEGA/BEAST

def process_populations_to_nexus(base_path: str, output_base: str):
    """
    Walk species directories under base_path, read FASTA + .grp, and write
    per-population NEXUS files to output_base as:
        sequences/Species_POP_Nexus/<Species>/GpX/<Species>_GpX.nexus
    """
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
            grp_path   = species_folder / "MEGA_Populations" / f"{species_name}_population.grp"

            if not fasta_path.exists() or not grp_path.exists():
                print(f"Skipping {species_folder} (missing FASTA or GRP file)")
                continue

            mapping = read_grp(grp_path)
            if not mapping:
                print(f"Skipping {species_folder} (empty or invalid GRP mapping)")
                continue

            seqs = read_fasta(fasta_path)
            if not seqs:
                print(f"Skipping {species_folder} (no sequences in FASTA)")
                continue

            species_out = output_base / species_folder.name
            species_out.mkdir(parents=True, exist_ok=True)

            # Population-level results folder only
            (species_out / "BSP_Population_Level_Results").mkdir(exist_ok=True)

            # Group sequences by population (as per mapping)
            pop_to_seqs = {}
            for name, seq in seqs:
                pop = mapping.get(name)
                if not pop:
                    # Sample not present in .grp mapping; skip it.
                    continue
                pop_to_seqs.setdefault(pop, []).append((name, seq))

            # Write per-population NEXUS with GpX naming everywhere
            for pop, pop_seqs in pop_to_seqs.items():
                if not pop_seqs:
                    continue
                gp = gp_from(pop)  # e.g., "Gp1"
                full_out_folder = species_out / gp
                full_out_folder.mkdir(parents=True, exist_ok=True)

                # filename uses species + GpX (species only once)
                nexus_out = full_out_folder / f"{species_name}_{gp}.nexus"

                # optional label: keep it short (GpX)
                write_nexus(nexus_out, pop_seqs, label=gp)

            print(f"Wrote population NEXUS files for {species_name}")

    print("Done.")

# --- paths ---
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Nexus"

if __name__ == "__main__":
    Path(output_base).mkdir(parents=True, exist_ok=True)
    process_populations_to_nexus(input_base, output_base)

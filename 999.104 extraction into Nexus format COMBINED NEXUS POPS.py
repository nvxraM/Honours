import os
import re
from pathlib import Path

STANDARD_BASES = set("ATCG-")

def replace_ambiguous_bases(seq: str) -> str:
    return ''.join(base if base.upper() in STANDARD_BASES else 'N' for base in seq)

def sanitize_name(name: str) -> str:
    # Nexus-friendly: spaces -> underscores; strip punctuation some parsers dislike
    n = re.sub(r"\s+", "_", name.strip())
    n = re.sub(r"[,\[\]\(\)';:]", "_", n)
    return n

def read_fasta(fasta_path: Path):
    """
    Read FASTA; returns list of (orig_name, safe_name, cleaned_seq).
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
                    seqs.append((name, sanitize_name(name), replace_ambiguous_bases(''.join(seq))))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if name is not None:
            seqs.append((name, sanitize_name(name), replace_ambiguous_bases(''.join(seq))))
    return seqs

def write_nexus_species_only(out_path: Path, named_seqs):
    """
    Write ONE NEXUS per species with ONLY the DATA block (no TRAITS/SETS/populations).
    """
    if not named_seqs:
        raise ValueError("No sequences to write.")

    # Ensure aligned
    lengths = {len(seq) for _, _, seq in named_seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not the same length: found lengths {sorted(lengths)}")
    nchar = lengths.pop()
    ntax = len(named_seqs)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("#NEXUS\n\n")
        f.write("Begin DATA;\n")
        f.write(f"    Dimensions ntax={ntax} nchar={nchar};\n")
        f.write('    Format datatype=dna symbols="ACTGN" missing=? gap=-;\n')
        f.write("Matrix\n")
        for _, safe_name, seq in named_seqs:
            f.write(f"{safe_name:30s} {seq}\n")
        f.write(";\nEnd;\n")

def process_species_to_combined_nexus(base_path: str, output_base: str):
    """
    Walk species directories under base_path, read FASTA, and write ONE
    species-level NEXUS (no population info) to:
        sequences/Species_POP_Nexus_Combined/<Species>/<Species>_combined.nexus

    Also creates:
        sequences/Species_POP_Nexus_Combined/<Species>/BSP_Species_Level_Results/
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

            if not fasta_path.exists():
                print(f"Skipping {species_folder} (missing FASTA: {fasta_path.name})")
                continue

            seqs = read_fasta(fasta_path)
            if not seqs:
                print(f"Skipping {species_folder} (no sequences in FASTA)")
                continue

            # Keep all sequences; no population filtering/labels
            named_seqs = seqs

            # Output paths
            species_out = output_base / species_name
            species_out.mkdir(parents=True, exist_ok=True)

            # Requested extra directory inside each species folder
            bsp_species_results = species_out / "BSP_Species_Level_Results"
            bsp_species_results.mkdir(parents=True, exist_ok=True)

            nexus_out = species_out / f"{species_name}_combined.nexus"
            write_nexus_species_only(nexus_out, named_seqs)
            print(f"Wrote species-level NEXUS (no populations) for {species_name}: {nexus_out}")
            print(f"Ensured results directory exists: {bsp_species_results}")

    print("Done.")

# --- paths ---
input_base = "sequences/Interactive_group_editor"
output_base = "sequences/Species_POP_Nexus_Combined"

if __name__ == "__main__":
    Path(output_base).mkdir(parents=True, exist_ok=True)
    process_species_to_combined_nexus(input_base, output_base)

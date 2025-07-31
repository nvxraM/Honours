import os
from pathlib import Path
from Bio import AlignIO
from collections import Counter, defaultdict
import numpy as np
import csv
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# ---- Haplotype Diversity (Hd) ----
# Hd = (n/(n-1)) * (1 - sum(pi^2))
# n = number of individuals in the population
# pi = frequency of each haplotype (proportion of individuals with that haplotype)
# This is Nei's gene diversity: the probability that two randomly chosen haplotypes are different.
def haplotype_diversity(hap_counts):
    n = sum(hap_counts)
    if n <= 1:
        return 0.0
    freq = np.array(hap_counts) / n
    Hd = (n / (n-1)) * (1 - np.sum(freq ** 2))
    return Hd

# ---- Nucleotide Diversity (Pi) ----
# Pi = average number of nucleotide differences per site between any two DNA sequences in the sample
# For each pair of sequences:
#   - Count the number of nucleotide differences
#   - Normalize by the number of sites compared (ignoring ambiguous/gap sites)
#   - Average over all possible pairs
def nucleotide_diversity(seqs):
    n = len(seqs)
    if n <= 1:
        return 0.0
    L = len(seqs[0])
    total_pi = 0
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            # Only count positions where both bases are A/T/G/C (ignore gaps/ambiguous)
            diffs = sum(a != b for a, b in zip(seqs[i], seqs[j]) if a in "ATGC" and b in "ATGC")
            length = sum((a in "ATGC" and b in "ATGC") for a, b in zip(seqs[i], seqs[j]))
            if length > 0:
                total_pi += diffs / length
                count += 1
    if count == 0:
        return 0.0
    # Average over all pairs
    return total_pi / count

# ---- Number of Variable (Segregating) Sites (VS) ----
# This is simply the count of alignment columns (positions) where there is more than one unique nucleotide (A/T/G/C)
def variable_sites(seqs):
    nsites = 0
    alignment = np.array([list(s) for s in seqs])
    for col in alignment.T:
        bases = set([x for x in col if x in "ATGC"])
        if len(bases) > 1:
            nsites += 1
    return nsites

# ---- Main per-file processing ----
def process_nexus(nexus_path):
    try:
        aln = AlignIO.read(nexus_path, "nexus")
    except Exception as e:
        print(f"Error reading {nexus_path}: {e}")
        return []
    pops = defaultdict(list)
    for rec in aln:
        # Assumes the population is specified as the last underscore-separated element (e.g., ..._1)
        if "_" in rec.id:
            pop = rec.id.split("_")[-1]
            pops[pop].append(str(rec.seq).upper())
    species = Path(nexus_path).stem
    pop_results = []
    for pop, seqs in pops.items():
        n = len(seqs)
        vs = variable_sites(seqs)
        # ---- Haplotype Number (HN) ----
        # HN is the number of unique haplotypes (unique sequences) in the population
        haps = Counter(seqs)
        hn = len(haps)
        hd = haplotype_diversity(list(haps.values()))
        pi = nucleotide_diversity(seqs)
        pop_results.append({
            "Species": species,
            "Population": pop,
            "N": n,  # Number of individuals/sequences
            "Variable_Sites": vs,  # VS: Number of segregating sites
            "Haplotype_Number": hn,  # HN: Number of unique haplotypes
            "Haplotype_Diversity": f"{hd:.8f}",  # Hd: Probability two sequences have different haplotypes
            "Nucleotide_Diversity": f"{pi:.8f}"  # Pi: Mean proportion of differences per site
        })
    return pop_results

def main():
    folder = Path("sequences/Species_POP_Nexus")
    all_results = []
    nexus_files = list(folder.glob("*.nex"))
    print(f"Processing {len(nexus_files)} .nex files (multi-threaded)...")
    max_threads = 25  # ~80% of your 32 logical threads
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_file = {executor.submit(process_nexus, f): f for f in nexus_files}
        for future in tqdm(as_completed(future_to_file), total=len(nexus_files), desc="Analyzing species", unit="species"):
            res = future.result()
            all_results.extend(res)
    # Write to CSV
    out_csv = folder / "population_statistics.csv"
    with open(out_csv, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Species", "Population", "N", "Variable_Sites",
            "Haplotype_Number", "Haplotype_Diversity", "Nucleotide_Diversity"
        ])
        writer.writeheader()
        for row in all_results:
            writer.writerow(row)
    print(f"Results written to {out_csv}")

if __name__ == "__main__":
    main()

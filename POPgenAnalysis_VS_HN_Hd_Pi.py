import os
from pathlib import Path
from Bio import AlignIO
from collections import Counter, defaultdict
import numpy as np
import csv
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import sqrt

# ----- Haplotype Diversity (Hd) -----
# Hd (Nei's 1987 gene diversity) is the probability that two randomly selected haplotypes are different.
# It measures the uniqueness of haplotypes in your sample, i.e., the genetic diversity at the haplotype level.
# Formula: Hd = (n/(n-1)) * (1 - sum(pi^2))
#    - n: number of individuals/sequences
#    - pi: frequency (proportion) of each haplotype
def haplotype_diversity(hap_counts):
    n = sum(hap_counts)
    if n <= 1:
        return 0.0
    freq = np.array(hap_counts) / n  # Calculate frequencies of each haplotype
    Hd = (n / (n-1)) * (1 - np.sum(freq ** 2))
    return Hd

# ----- Standard Deviation of Haplotype Diversity (Hd) -----
# This describes how much Hd would vary between different samples of the same size from this population.
# Formula (common approximation, matches DnaSP): Var(Hd) = Hd*(1-Hd)/(n-1)
# SD(Hd) = sqrt(Var(Hd))
def hd_std(hd, n):
    if n <= 1:
        return 0.0
    var = (hd * (1 - hd)) / (n - 1)
    return sqrt(var) if var >= 0 else 0.0

# ----- Nucleotide Diversity (Pi) -----
# Pi (Nei & Li, 1979) measures the average proportion of nucleotide differences per site between all possible pairs of sequences.
# 1. For each pair of sequences, count the number of differences and divide by the number of valid (A/T/G/C) sites.
# 2. Average across all pairs.
def nucleotide_diversity(seqs):
    n = len(seqs)
    if n <= 1:
        return 0.0, []
    pairwise_pis = []
    for i in range(n):
        for j in range(i+1, n):
            # Only compare positions where both have a clear base (A/T/G/C)
            diffs = sum(a != b for a, b in zip(seqs[i], seqs[j]) if a in "ATGC" and b in "ATGC")
            length = sum((a in "ATGC" and b in "ATGC") for a, b in zip(seqs[i], seqs[j]))
            if length > 0:
                pairwise_pis.append(diffs / length)
    if not pairwise_pis:
        return 0.0, []
    pi = np.mean(pairwise_pis)
    return pi, pairwise_pis

# ----- Standard Deviation of Nucleotide Diversity (Pi) -----
# This quantifies the variability in Pi among all possible pairs of sequences.
# It's simply the standard deviation of the list of all pairwise Pi values.
def pi_std(pairwise_pis):
    if len(pairwise_pis) < 2:
        return 0.0
    return float(np.std(pairwise_pis, ddof=1))

# ----- Number of Variable Sites (VS) -----
# The number of alignment columns (sites) with more than one type of nucleotide (A/T/G/C) in your sample.
def variable_sites(seqs):
    nsites = 0
    alignment = np.array([list(s) for s in seqs])
    for col in alignment.T:
        bases = set([x for x in col if x in "ATGC"])
        if len(bases) > 1:
            nsites += 1
    return nsites

# ----- Per-file Processing -----
# Reads a .nex file, assigns sequences to populations, and computes all stats for each population.
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
        # ----- Haplotype Number (HN) -----
        # The count of unique haplotype sequences in the population.
        haps = Counter(seqs)
        hn = len(haps)
        hd = haplotype_diversity(list(haps.values()))
        hd_sd = hd_std(hd, n)
        pi, pairwise_pis = nucleotide_diversity(seqs)
        pi_sd = pi_std(pairwise_pis)
        pop_results.append({
            "Species": species,
            "Population": pop,
            "N": n,  # Number of sequences in this population
            "Variable_Sites": vs,  # Number of variable (segregating) sites
            "Haplotype_Number": hn,  # Number of unique haplotypes (HN)
            "Haplotype_Diversity": f"{hd:.8f}",  # Hd (probability two are different)
            "Hd_StdDev": f"{hd_sd:.8f}",        # Standard deviation of Hd
            "Nucleotide_Diversity": f"{pi:.8f}",  # Pi (average pairwise difference)
            "Pi_StdDev": f"{pi_sd:.8f}"           # Standard deviation of Pi
        })
    return pop_results

def main():
    # Directory containing your .nex files
    folder = Path("sequences/Species_POP_Nexus")
    all_results = []
    nexus_files = list(folder.glob("*.nex"))
    print(f"Processing {len(nexus_files)} .nex files (multi-threaded)...")
    max_threads = 25  # Use 80% of 32 logical cores
    # Run processing in parallel threads for speed
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_file = {executor.submit(process_nexus, f): f for f in nexus_files}
        for future in tqdm(as_completed(future_to_file), total=len(nexus_files), desc="Analyzing species", unit="species"):
            res = future.result()
            all_results.extend(res)
    # Write results to CSV file
    out_csv = folder / "zpopulation_statistics.csv"
    with open(out_csv, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Species", "Population", "N", "Variable_Sites",
            "Haplotype_Number", "Haplotype_Diversity", "Hd_StdDev",
            "Nucleotide_Diversity", "Pi_StdDev"
        ])
        writer.writeheader()
        for row in all_results:
            writer.writerow(row)
    print(f"Results written to {out_csv}")

if __name__ == "__main__":
    main()

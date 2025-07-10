#!/usr/bin/env python3
"""
score_csr_full_report_with_reference.py

1) For each CSR .afa alignment under sequences/CDS_CSR/<Genus>/aligned/:
     • Use the first record as the “closest‐relative” reference.
     • For each other record (samples), compute:
         – Synonymous mutation count
         – Nonsynonymous mutation count
         – Analyzed codon count
         – Positions of synonymous mutations (1-based codon index[SAMPLE_CODON])
         – Positions of nonsynonymous mutations (1-based codon index[SAMPLE_CODON])
2) Per species: record the reference ID, sum all sample synonymous-counts,
   nonsynonymous-counts, and analyzed-codon-counts for the summary report.
3) Grand sum across all individuals of all species for these S/NS metrics
   for the summary report.
4) Output:
   - CSR_species_summary.xlsx   (one row per species + an “All_species” row)
   - CSR_all_species_details.xlsx (one consolidated file, with each species' details on a separate sheet,
                                   sheet named after the species, including codon details in position strings)

Usage:
   python score_csr_full_report_with_reference.py
"""

import sys
import re  # For sanitizing sheet names
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# ─── CONFIG ────────────────────────────────────────────────────────────────

LOCAL_BASE = Path("sequences/CDS_CSR")
SUMMARY_XLSX = LOCAL_BASE / "CSR_species_summary.xlsx"
ALL_DETAILS_XLSX = LOCAL_BASE / "CSR_all_species_details.xlsx"

# --- Codon Definitions for S/NS analysis ---
SYN_TARGET_CODONS = frozenset([
    "TCT", "TCC", "TCA", "TCG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
    "CGT", "CGC", "CGA", "CGG", "ACT", "ACC", "ACA", "ACG", "GTT", "GTC", "GTA", "GTG",
    "GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG"
])

NONSYN_SPECIFIC_CODONS = frozenset([
    "CTA", "CTG", "CGA", "CGG", "TTA", "TTG"
])

# ─── HELPERS ───────────────────────────────────────────────────────────────

def sanitize_sheet_name(name: str) -> str:
    """Sanitizes a string to be a valid Excel sheet name."""
    name = re.sub(r'[\\/*?:\[\]]', '', name)
    return name[:31]


def calculate_syn_nonsyn_counts(ref_seq_str: str, samp_seq_str: str):
    """
    Calculates synonymous and nonsynonymous mutation counts and their positions.
    Returns (syn_count, nonsyn_count, analyzed_codons_count,
             syn_positions_list_with_codons, nonsyn_positions_list_with_codons).
    Positions are 1-based codon indices. Codons are the sample codons.
    Each element in the position lists is a tuple (position, sample_codon).
    """
    syn_count = 0
    nonsyn_count = 0
    analyzed_codons_count = 0
    syn_positions_codons = []
    nonsyn_positions_codons = []

    ref_s = ref_seq_str.upper()
    samp_s = samp_seq_str.upper()
    min_len = min(len(ref_s), len(samp_s))
    num_codons = min_len // 3

    for i_codon in range(num_codons):
        idx = i_codon * 3
        ref_c = ref_s[idx: idx + 3]
        samp_c = samp_s[idx: idx + 3]

        # Validate codon (no gaps or Ns)
        valid = len(ref_c) == 3 and len(samp_c) == 3 and all(b not in '-N' for b in ref_c + samp_c)
        if not valid:
            continue

        analyzed_codons_count += 1
        pos = i_codon + 1

        if ref_c == samp_c:
            continue

        # Synonymous
        if ref_c in SYN_TARGET_CODONS and ref_c[:2] == samp_c[:2] and ref_c[2] != samp_c[2]:
            syn_count += 1
            syn_positions_codons.append((pos, samp_c))
            continue

        # Nonsynonymous specific
        if ref_c in NONSYN_SPECIFIC_CODONS and ref_c[0] == samp_c[0] and ref_c[2] == samp_c[2] and ref_c[1] != samp_c[1]:
            nonsyn_count += 1
            nonsyn_positions_codons.append((pos, samp_c))
            continue

        # General nonsynonymous
        if ref_c[1] != samp_c[1] or ref_c[2] != samp_c[2]:
            nonsyn_count += 1
            nonsyn_positions_codons.append((pos, samp_c))

    return syn_count, nonsyn_count, analyzed_codons_count, syn_positions_codons, nonsyn_positions_codons


def process_alignment(afa_path: Path, summary_rows: list):
    genus = afa_path.parent.parent.name
    species = afa_path.stem.removeprefix("aligned").lstrip("_")
    records = list(SeqIO.parse(str(afa_path), "fasta"))

    if len(records) < 2:
        return None, None

    ref_rec = records[0]
    ref_id = ref_rec.id
    ref_seq = str(ref_rec.seq)

    raw_syn = []
    raw_nonsyn = []
    raw_codon_counts = []
    detail_rows = []

    # Reference row
    _, _, ref_codons, _, _ = calculate_syn_nonsyn_counts(ref_seq, ref_seq)
    detail_rows.append({
        "Role": "Reference", "Sequence_ID": ref_id,
        "Synonymous_Count": 0, "Nonsynonymous_Count": 0,
        "Analyzed_Codons": ref_codons,
        "Synonymous_Mutation_Positions": "", "Nonsynonymous_Mutation_Positions": ""
    })

    # Sample rows
    for rec in records[1:]:
        syn, nonsyn, codons, syn_pos, nonsyn_pos = calculate_syn_nonsyn_counts(ref_seq, str(rec.seq))
        raw_syn.append(syn)
        raw_nonsyn.append(nonsyn)
        raw_codon_counts.append(codons)

        detail_rows.append({
            "Role": "Sample", "Sequence_ID": rec.id,
            "Synonymous_Count": syn, "Nonsynonymous_Count": nonsyn,
            "Analyzed_Codons": codons,
            "Synonymous_Mutation_Positions": ",".join(f"{p}[{c}]" for p, c in syn_pos),
            "Nonsynonymous_Mutation_Positions": ",".join(f"{p}[{c}]" for p, c in nonsyn_pos)
        })

    n_samples = len(raw_syn)
    if n_samples > 0:
        sum_syn = sum(raw_syn)
        sum_nonsyn = sum(raw_nonsyn)
        sum_codons = sum(raw_codon_counts)
        summary_rows.append({
            "Genus": genus, "Species": species, "Reference_ID": ref_id,
            "Num_Samples": n_samples,
            "Sum_Syn_Count": sum_syn,
            "Sum_Nonsyn_Count": sum_nonsyn,
            "Sum_Analyzed_Codons": sum_codons
        })

    detail_df = pd.DataFrame(detail_rows)
    cols = ["Role", "Sequence_ID", "Synonymous_Count", "Nonsynonymous_Count", "Analyzed_Codons",
            "Synonymous_Mutation_Positions", "Nonsynonymous_Mutation_Positions"]
    return species, detail_df[cols]


def main():
    summary = []
    details = {}
    LOCAL_BASE.mkdir(parents=True, exist_ok=True)

    for genus_dir in LOCAL_BASE.iterdir():
        if not genus_dir.is_dir(): continue
        aligned = genus_dir / "aligned"
        if not aligned.is_dir(): continue
        for afa in aligned.glob("*.afa"):
            print(f"Processing {afa}...")
            name, df = process_alignment(afa, summary)
            if name and df is not None:
                details[name] = df

    # Write summary
    if summary:
        df_sum = pd.DataFrame(summary)
        df_sum.sort_values(["Genus", "Species"], inplace=True)
        total_samples = df_sum["Num_Samples"].sum()
        total_syn = df_sum["Sum_Syn_Count"].sum()
        total_nonsyn = df_sum["Sum_Nonsyn_Count"].sum()
        total_codons = df_sum["Sum_Analyzed_Codons"].sum()
        all_row = {
            "Genus": "All", "Species": "All_species", "Reference_ID": "",
            "Num_Samples": total_samples,
            "Sum_Syn_Count": total_syn,
            "Sum_Nonsyn_Count": total_nonsyn,
            "Sum_Analyzed_Codons": total_codons
        }
        df_sum = pd.concat([df_sum, pd.DataFrame([all_row])], ignore_index=True)
        cols = ["Genus", "Species", "Reference_ID", "Num_Samples",
                "Sum_Syn_Count", "Sum_Nonsyn_Count", "Sum_Analyzed_Codons"]
        df_sum[cols].to_excel(SUMMARY_XLSX, index=False)
        print(f"✔ Wrote summary to {SUMMARY_XLSX}")
    else:
        print("No summary data generated.", file=sys.stderr)

    # Write details
    if details:
        with pd.ExcelWriter(ALL_DETAILS_XLSX, engine='openpyxl') as writer:
            for sp in sorted(details):
                details[sp].to_excel(writer, sheet_name=sanitize_sheet_name(sp), index=False)
        print(f"✔ Wrote details to {ALL_DETAILS_XLSX}")
    else:
        print("No detail data generated.", file=sys.stderr)

if __name__ == "__main__":
    main()

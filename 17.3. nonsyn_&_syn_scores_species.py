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
2) Per species: record the reference ID, average all sample synonymous-counts,
   nonsynonymous-counts, and analyzed-codon-counts for the summary report.
3) Grand (weighted) average across all individuals of all species for these S/NS metrics
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
    syn_positions_codons = []  # Will store (position, sample_codon)
    nonsyn_positions_codons = []  # Will store (position, sample_codon)

    ref_s = ref_seq_str.upper()
    samp_s = samp_seq_str.upper()
    min_len = min(len(ref_s), len(samp_s))
    num_codons = min_len // 3

    for i_codon in range(num_codons):
        idx = i_codon * 3
        ref_c = ref_s[idx: idx + 3]
        samp_c = samp_s[idx: idx + 3]

        is_codon_valid_for_sns = True
        if len(ref_c) < 3 or len(samp_c) < 3:
            is_codon_valid_for_sns = False
        else:
            for k_base in range(3):
                if ref_c[k_base] in '-N' or samp_c[k_base] in '-N':
                    is_codon_valid_for_sns = False
                    break

        if not is_codon_valid_for_sns:
            continue

        analyzed_codons_count += 1
        codon_pos_1_based = i_codon + 1

        if ref_c == samp_c:
            continue

        mutation_counted_this_codon = False

        if ref_c in SYN_TARGET_CODONS:
            if ref_c[:2] == samp_c[:2] and ref_c[2] != samp_c[2]:
                syn_count += 1
                syn_positions_codons.append((codon_pos_1_based, samp_c))  # Store pos and sample codon
                mutation_counted_this_codon = True

        if mutation_counted_this_codon:
            continue

        if ref_c in NONSYN_SPECIFIC_CODONS:
            if ref_c[0] == samp_c[0] and ref_c[2] == samp_c[2] and ref_c[1] != samp_c[1]:
                nonsyn_count += 1
                nonsyn_positions_codons.append((codon_pos_1_based, samp_c))  # Store pos and sample codon
                mutation_counted_this_codon = True

        if mutation_counted_this_codon:
            continue

        condition_not_both_special = not (ref_c in NONSYN_SPECIFIC_CODONS and samp_c in NONSYN_SPECIFIC_CODONS)

        if condition_not_both_special:
            if ref_c[1] != samp_c[1] or ref_c[2] != samp_c[2]:
                nonsyn_count += 1
                nonsyn_positions_codons.append((codon_pos_1_based, samp_c))  # Store pos and sample codon

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

    raw_syn_counts_summary = []
    raw_nonsyn_counts_summary = []
    raw_analyzed_codons_summary = []
    current_species_detail_rows = []

    _, _, ref_analyzed_codons, _, _ = calculate_syn_nonsyn_counts(ref_seq, ref_seq)
    current_species_detail_rows.append({
        "Role": "Reference", "Sequence_ID": ref_id,
        "Synonymous_Count": 0, "Nonsynonymous_Count": 0,
        "Analyzed_Codons": ref_analyzed_codons,
        "Synonymous_Mutation_Positions": "", "Nonsynonymous_Mutation_Positions": ""
    })

    for rec in records[1:]:
        samp_seq_str = str(rec.seq)
        # syn_p_codons and nonsyn_p_codons will be lists of (pos, codon_str) tuples
        syn_c, nonsyn_c, analyzed_c, syn_p_codons, nonsyn_p_codons = calculate_syn_nonsyn_counts(ref_seq, samp_seq_str)

        raw_syn_counts_summary.append(syn_c)
        raw_nonsyn_counts_summary.append(nonsyn_c)
        raw_analyzed_codons_summary.append(analyzed_c)

        # Format position strings with codons
        syn_pos_str = ",".join([f"{pos}[{codon}]" for pos, codon in syn_p_codons])
        nonsyn_pos_str = ",".join([f"{pos}[{codon}]" for pos, codon in nonsyn_p_codons])

        current_species_detail_rows.append({
            "Role": "Sample", "Sequence_ID": rec.id,
            "Synonymous_Count": syn_c, "Nonsynonymous_Count": nonsyn_c,
            "Analyzed_Codons": analyzed_c,
            "Synonymous_Mutation_Positions": syn_pos_str,
            "Nonsynonymous_Mutation_Positions": nonsyn_pos_str
        })

    n_samples = len(raw_syn_counts_summary)
    if n_samples > 0:
        avg_syn_count = sum(raw_syn_counts_summary) / n_samples
        avg_nonsyn_count = sum(raw_nonsyn_counts_summary) / n_samples
        avg_analyzed_codons = sum(raw_analyzed_codons_summary) / n_samples
        summary_rows.append({
            "Genus": genus, "Species": species, "Reference_ID": ref_id,
            "Num_Samples": n_samples,
            "Avg_Syn_Count": round(avg_syn_count, 2),
            "Avg_Nonsyn_Count": round(avg_nonsyn_count, 2),
            "Avg_Analyzed_Codons": round(avg_analyzed_codons, 2),
        })

    species_detail_df = None
    if current_species_detail_rows:
        df = pd.DataFrame(current_species_detail_rows)
        detail_cols_order = [
            "Role", "Sequence_ID",
            "Synonymous_Count", "Nonsynonymous_Count", "Analyzed_Codons",
            "Synonymous_Mutation_Positions", "Nonsynonymous_Mutation_Positions"
        ]
        existing_cols = [col for col in detail_cols_order if col in df.columns]
        species_detail_df = df[existing_cols]
    return species, species_detail_df


# ─── MAIN ─────────────────────────────────────────────────────────────────

def main():
    summary_rows = []
    all_species_detail_data = {}

    LOCAL_BASE.mkdir(parents=True, exist_ok=True)

    for genus_dir in LOCAL_BASE.iterdir():
        if not genus_dir.is_dir(): continue
        aligned_dir = genus_dir / "aligned"
        if not aligned_dir.is_dir(): continue
        for afa_path in aligned_dir.glob("*.afa"):
            print(f"Processing {afa_path}...")
            species_name, species_df = process_alignment(afa_path, summary_rows)
            if species_name and species_df is not None and not species_df.empty:
                if species_name in all_species_detail_data:
                    print(
                        f"Warning: Overwriting detail data for species {species_name}. Multiple .afa files might map to it.",
                        file=sys.stderr)
                all_species_detail_data[species_name] = species_df

    if not summary_rows:
        print("No summary data generated from .afa files.", file=sys.stderr)
    else:
        df_sum = pd.DataFrame(summary_rows)
        df_sum.sort_values(["Genus", "Species"], inplace=True)
        total_samples = df_sum["Num_Samples"].sum()
        if total_samples > 0:
            weighted_avg_syn = (df_sum["Avg_Syn_Count"] * df_sum["Num_Samples"]).sum() / total_samples
            weighted_avg_nonsyn = (df_sum["Avg_Nonsyn_Count"] * df_sum["Num_Samples"]).sum() / total_samples
            weighted_avg_analyzed_codons = (df_sum["Avg_Analyzed_Codons"] * df_sum["Num_Samples"]).sum() / total_samples
        else:
            weighted_avg_syn, weighted_avg_nonsyn, weighted_avg_analyzed_codons = 0, 0, 0
        all_species_summary_row = {
            "Genus": "All", "Species": "All_species", "Reference_ID": "",
            "Num_Samples": total_samples,
            "Avg_Syn_Count": round(weighted_avg_syn, 2),
            "Avg_Nonsyn_Count": round(weighted_avg_nonsyn, 2),
            "Avg_Analyzed_Codons": round(weighted_avg_analyzed_codons, 2),
        }
        df_sum = pd.concat([df_sum, pd.DataFrame([all_species_summary_row])], ignore_index=True)
        summary_cols_order = [
            "Genus", "Species", "Reference_ID", "Num_Samples",
            "Avg_Syn_Count", "Avg_Nonsyn_Count", "Avg_Analyzed_Codons"
        ]
        df_sum = df_sum.reindex(columns=summary_cols_order)
        df_sum.to_excel(SUMMARY_XLSX, index=False)
        print(f"✔ Wrote summary to {SUMMARY_XLSX}")

    if not all_species_detail_data:
        print("No detail data generated for any species to write to Excel.", file=sys.stderr)
    else:
        with pd.ExcelWriter(ALL_DETAILS_XLSX, engine='openpyxl') as writer:
            sorted_species_names = sorted(all_species_detail_data.keys())
            for species_name in sorted_species_names:
                df_detail_species = all_species_detail_data[species_name]
                safe_sheet_name = sanitize_sheet_name(species_name)
                df_detail_species.to_excel(writer, sheet_name=safe_sheet_name, index=False)
        print(f"✔ Wrote all species details to {ALL_DETAILS_XLSX}")


if __name__ == "__main__":
    main()
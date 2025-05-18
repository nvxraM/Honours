#!/usr/bin/env python3
"""
score_csr_full_report_with_reference.py

1) For each CSR .afa alignment under sequences/CDS_CSR/<Genus>/aligned/:
     • Use the first record as the “closest‐relative” reference.
     • For each other record (samples), compute:
         – Valid positions (skip '-' or 'N')
         – Mismatch positions
         – Percent identity = matches/valid * 100
         – Mismatch count
2) Per species: record the reference ID, average all sample percent‐IDs and mismatch‐counts,
   collect raw per‐sample percent‐IDs into a comma‐separated string.
3) Grand (weighted) average across all individuals of all species for both metrics.
4) Output:
   - CSR_species_summary.xlsx   (one row per species + an “All_species” row), including Reference_ID
   - CSR_species_detail_<SPECIES>.xlsx (full per‐sequence breakdown)

Usage:
   python score_csr_full_report_with_reference.py
"""

import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# ─── CONFIG ────────────────────────────────────────────────────────────────

LOCAL_BASE     = Path("sequences/CDS_CSR")
SUMMARY_XLSX   = LOCAL_BASE / "CSR_species_summary.xlsx"

DETAIL_SPECIES = "Balaena_Balaena_mysticetus"
DETAIL_XLSX    = LOCAL_BASE / f"CSR_species_detail_{DETAIL_SPECIES}.xlsx"

# ─── HELPERS ───────────────────────────────────────────────────────────────

def percent_identity_and_mismatches(ref_seq: str, samp_seq: str):
    """
    Walk two aligned sequences column-by-column (same length).
    Skip any column where either char is '-' or 'N'.
    Return (percent_id, mismatch_count, mismatch_positions_list).
    """
    r = ref_seq.upper()
    s = samp_seq.upper()
    valid = matches = 0
    mismatches = []
    for i, (a, b) in enumerate(zip(r, s), start=1):
        if a in "-N" or b in "-N":
            continue
        valid += 1
        if a == b:
            matches += 1
        else:
            mismatches.append(i)
    pid = (matches / valid * 100) if valid else 0.0
    return round(pid, 2), len(mismatches), mismatches

def process_alignment(afa_path, detail_rows, summary_rows):
    genus     = afa_path.parent.parent.name
    species   = afa_path.stem.removeprefix("aligned").lstrip("_")
    records   = list(SeqIO.parse(str(afa_path), "fasta"))
    if len(records) < 2:
        return

    # reference record is the first one
    ref_rec  = records[0]
    ref_id   = ref_rec.id
    ref_seq  = str(ref_rec.seq)

    raw_pids    = []
    raw_mcounts = []

    # Detail: include reference row
    if species == DETAIL_SPECIES:
        valid_ref = sum(1 for c in ref_seq.upper() if c not in "-N")
        detail_rows.append({
            "Role":               "Reference",
            "Sequence_ID":        ref_id,
            "Sequence":           ref_seq,
            "Valid_Positions":    valid_ref,
            "Mismatch_Count":     0,
            "Percent_ID":         100.0,
            "Mismatch_Positions": ""
        })

    # Score each sample against the reference
    for rec in records[1:]:
        pid, mcount, mpos = percent_identity_and_mismatches(ref_seq, str(rec.seq))
        raw_pids.append(pid)
        raw_mcounts.append(mcount)
        if species == DETAIL_SPECIES:
            detail_rows.append({
                "Role":               "Sample",
                "Sequence_ID":        rec.id,
                "Sequence":           str(rec.seq),
                "Valid_Positions":    None,
                "Mismatch_Count":     mcount,
                "Percent_ID":         pid,
                "Mismatch_Positions": ",".join(map(str, mpos))
            })

    n = len(raw_pids)
    if n == 0:
        return

    avg_pid     = sum(raw_pids)    / n
    avg_mcount  = sum(raw_mcounts) / n
    avg_mpercent= round(100 - avg_pid, 2)

    summary_rows.append({
        "Genus":                genus,
        "Species":              species,
        "Reference_ID":         ref_id,
        "Num_Samples":          n,
        "Avg_Percent_ID":       round(avg_pid, 2),
        "Avg_Mismatch_Count":   round(avg_mcount, 2),
        "Avg_Mismatch_Percent": avg_mpercent,
        "Raw_Percent_IDs":      ", ".join(f"{x:.2f}" for x in raw_pids)
    })

# ─── MAIN ─────────────────────────────────────────────────────────────────

def main():
    summary_rows = []
    detail_rows  = []

    for genus_dir in LOCAL_BASE.iterdir():
        if not genus_dir.is_dir():
            continue
        aligned_dir = genus_dir / "aligned"
        if not aligned_dir.is_dir():
            continue
        for afa in aligned_dir.glob("*.afa"):
            process_alignment(afa, detail_rows, summary_rows)

    if not summary_rows:
        print("No data found.", file=sys.stderr)
        sys.exit(1)

    # Build summary DataFrame
    df_sum = pd.DataFrame(summary_rows)
    df_sum.sort_values("Avg_Percent_ID", ascending=False, inplace=True)

    # Compute weighted grand averages
    total_samples      = df_sum["Num_Samples"].sum()
    weighted_avg_pid   = (df_sum["Avg_Percent_ID"] * df_sum["Num_Samples"]).sum() / total_samples
    weighted_avg_mcount= (df_sum["Avg_Mismatch_Count"] * df_sum["Num_Samples"]).sum() / total_samples
    weighted_avg_mpct  = round(100 - weighted_avg_pid, 2)

    # Append “All_species” row (no reference for all)
    df_sum.loc[len(df_sum)] = {
        "Genus":                "All",
        "Species":              "All_species",
        "Reference_ID":         "",
        "Num_Samples":          total_samples,
        "Avg_Percent_ID":       round(weighted_avg_pid, 2),
        "Avg_Mismatch_Count":   round(weighted_avg_mcount, 2),
        "Avg_Mismatch_Percent": weighted_avg_mpct,
        "Raw_Percent_IDs":      ""
    }

    # Write summary workbook
    df_sum.to_excel(SUMMARY_XLSX, index=False)
    print(f"✔ Wrote summary to {SUMMARY_XLSX}")

    # Write detail workbook for DETAIL_SPECIES
    if detail_rows:
        df_det = pd.DataFrame(detail_rows)
        df_det.to_excel(DETAIL_XLSX, index=False)
        print(f"✔ Wrote detail for {DETAIL_SPECIES} to {DETAIL_XLSX}")
    else:
        print(f"No detail rows for {DETAIL_SPECIES}", file=sys.stderr)

if __name__ == "__main__":
    main()

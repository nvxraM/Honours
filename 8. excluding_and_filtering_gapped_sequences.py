import os
import shutil
from Bio import SeqIO
from collections import defaultdict

def exclude_gapped_and_enforce_minimum(base_dir, min_sequences=5):
    """
    For each species directory under the base directory (sequences/CDS/[Species]):
    1. Identify the concatenated gapped sequences FASTA file in:
       sequences/CDS/[Species]/CDS_nucleotide_gapped/[Species]_concatenated_gapped_sequences.fasta.
    2. Check each sequence for gaps ('-'):
       - Write sequences containing gaps to an excluded FASTA file under `sequences/excluded/[Genus]/`.
       - Keep sequences without gaps in the original concatenated FASTA file.
    3. Remove sequences with gaps from `sequences/fasta/[Species].fasta`.
    4. Remove sequences with gaps from individual GenBank files in `sequences/gb/[Species]`.
    5. After removing gapped sequences, enforce a minimum number of sequences per species.
       If a species has fewer than `min_sequences` sequences total (based on `sequences/fasta/[Species].fasta`),
       move that species to the `excluded` directory and remove from main directories.
    6. Print a summary report of how many sequences each species has, grouped by genus.

    Directory structure assumed:
    sequences/
      ├─ CDS/
      │  ├─ [Species]/
      │  │  └─ CDS_nucleotide_gapped/
      │  │     └─ [Species]_concatenated_gapped_sequences.fasta
      ├─ fasta/
      │  ├─ [Species].fasta
      ├─ gb/
      │  ├─ [Species]/ (contains individual .gb files)
      ├─ excluded/
      │  └─ [Genus]/
    """

    cds_path = os.path.join(base_dir, "CDS")
    fasta_path = os.path.join(base_dir, "fasta")
    gb_path = os.path.join(base_dir, "gb")
    excluded_base = os.path.join(base_dir, "excluded")

    # List of species directories under CDS
    species_list = [d for d in os.listdir(cds_path) if os.path.isdir(os.path.join(cds_path, d))]

    # ------------------------------------------------------------
    # STEP 1 & 2: PROCESS EACH SPECIES TO EXCLUDE GAPPED SEQUENCES
    # ------------------------------------------------------------
    for species in species_list:
        species_dir = os.path.join(cds_path, species)
        gapped_dir = os.path.join(species_dir, "CDS_nucleotide_gapped")
        gapped_file = os.path.join(gapped_dir, f"{species}_concatenated_gapped_sequences.fasta")

        if not os.path.exists(gapped_file):
            # No concatenated file for this species
            continue

        records = list(SeqIO.parse(gapped_file, "fasta"))

        # Separate sequences into those with and without gaps
        no_gap_records = []
        gap_records = []

        for rec in records:
            if "-" in str(rec.seq):
                gap_records.append(rec)
            else:
                no_gap_records.append(rec)

        # If there are gapped sequences, exclude them
        if gap_records:
            # Attempt to extract genus from the first gapped record
            first_gap_id_parts = gap_records[0].id.split("_")
            if len(first_gap_id_parts) < 2:
                # If format isn't as expected, fallback to species name
                # Assuming species name is Genus_species
                genus_name = species.split("_")[0]
            else:
                genus_name = first_gap_id_parts[0]

            # Create excluded genus directory if it doesn't exist
            excluded_genus_dir = os.path.join(excluded_base, genus_name)
            os.makedirs(excluded_genus_dir, exist_ok=True)

            # Write gapped sequences to an excluded FASTA file
            excluded_gapped_fasta = os.path.join(excluded_genus_dir, f"excluded_{genus_name}_gapped_sequences.fasta")
            with open(excluded_gapped_fasta, "a") as excl_out:
                SeqIO.write(gap_records, excl_out, "fasta")

            # -----------------------------------------------------
            # STEP 3: REMOVE GAPPED SEQUENCES FROM [Species].fasta
            # -----------------------------------------------------
            species_fasta = os.path.join(fasta_path, f"{species}.fasta")
            if os.path.exists(species_fasta):
                all_species_recs = list(SeqIO.parse(species_fasta, "fasta"))
                gap_ids = set(r.id for r in gap_records)
                # Keep only those not in gap_ids
                retained_species_recs = [r for r in all_species_recs if r.id not in gap_ids]
                with open(species_fasta, "w") as out:
                    SeqIO.write(retained_species_recs, out, "fasta")

            # ----------------------------------------------------
            # STEP 4: REMOVE GAPPED SEQUENCES FROM GB FILES
            # ----------------------------------------------------
            gb_species_dir = os.path.join(gb_path, species)
            if os.path.exists(gb_species_dir):
                gap_ids = set(r.id for r in gap_records)
                # We assume each GB file corresponds to a single record ID
                for gb_file in os.listdir(gb_species_dir):
                    gb_id = os.path.splitext(gb_file)[0]
                    if gb_id in gap_ids:
                        os.remove(os.path.join(gb_species_dir, gb_file))

        # Rewrite the original gapped file with only no-gap sequences
        with open(gapped_file, "w") as out:
            SeqIO.write(no_gap_records, out, "fasta")

    # --------------------------------------------------------------------
    # STEP 5: ENFORCE MINIMUM NUMBER OF SEQUENCES PER SPECIES
    # After removing gapped sequences, we now ensure each species has at
    # least `min_sequences` sequences.
    # --------------------------------------------------------------------

    # Count sequences per species from their FASTA files
    species_counts = {}
    if os.path.exists(fasta_path):
        fasta_files = [f for f in os.listdir(fasta_path) if f.endswith(".fasta")]
    else:
        fasta_files = []

    for fasta_file in fasta_files:
        species_name = os.path.splitext(fasta_file)[0]
        fasta_file_path = os.path.join(fasta_path, fasta_file)
        records = list(SeqIO.parse(fasta_file_path, "fasta"))
        count = len(records)
        species_counts[species_name] = count

    # Now apply the minimum sequence rule
    for sp, count in species_counts.items():
        if count < min_sequences:
            # Need to exclude this species
            species_fasta = os.path.join(fasta_path, f"{sp}.fasta")
            if not os.path.exists(species_fasta):
                continue
            records = list(SeqIO.parse(species_fasta, "fasta"))

            # Extract genus name
            if records:
                # Try to get genus from the first record
                first_id_parts = records[0].id.split("_")
                if len(first_id_parts) < 2:
                    # Fallback: parse genus from species name
                    genus_name = sp.split("_")[0]
                else:
                    genus_name = first_id_parts[0]
            else:
                # No records left, fallback to species name
                genus_name = sp.split("_")[0]

            print(f"Species {sp} under genus {genus_name} has less than {min_sequences} samples. Moving to excluded directory.")

            # Create excluded genus directory
            excluded_genus_dir = os.path.join(excluded_base, genus_name)
            os.makedirs(excluded_genus_dir, exist_ok=True)

            # Move this species' records to an excluded FASTA file if there are any records
            if records:
                excluded_fasta = os.path.join(excluded_genus_dir, f"excluded_{genus_name}.fasta")
                with open(excluded_fasta, "a") as out_excl:
                    SeqIO.write(records, out_excl, "fasta")

            # Clear this species from the main species FASTA
            with open(species_fasta, "w") as out:
                pass

            # Move GB files
            gb_species_dir = os.path.join(gb_path, sp)
            if os.path.exists(gb_species_dir):
                excluded_gb_dir = os.path.join(excluded_genus_dir, "gb", sp)
                os.makedirs(excluded_gb_dir, exist_ok=True)
                for gb_file in os.listdir(gb_species_dir):
                    os.rename(
                        os.path.join(gb_species_dir, gb_file),
                        os.path.join(excluded_gb_dir, gb_file)
                    )
                # Remove the now-empty gb species directory
                os.rmdir(gb_species_dir)

            # Move CDS files
            cds_species_dir = os.path.join(cds_path, sp)
            if os.path.exists(cds_species_dir):
                excluded_cds_dir = os.path.join(excluded_genus_dir, "CDS", sp)
                os.makedirs(excluded_cds_dir, exist_ok=True)
                for file in os.listdir(cds_species_dir):
                    os.rename(
                        os.path.join(cds_species_dir, file),
                        os.path.join(excluded_cds_dir, file)
                    )
                # Remove the now-empty cds species directory
                os.rmdir(cds_species_dir)

    # --------------------------------------------------------------------
    # STEP 6: PRINT SUMMARY REPORT
    # Recount the sequences now that some may have been excluded.
    # --------------------------------------------------------------------
    final_counts = {}
    if os.path.exists(fasta_path):
        fasta_files = [f for f in os.listdir(fasta_path) if f.endswith(".fasta")]
        for fasta_file in fasta_files:
            species_name = os.path.splitext(fasta_file)[0]
            fasta_file_path = os.path.join(fasta_path, fasta_file)
            records = list(SeqIO.parse(fasta_file_path, "fasta"))
            count = len(records)

            # Derive genus from species name
            # Assuming species_name is in form Genus_species
            genus = species_name.split("_")[0]

            if genus not in final_counts:
                final_counts[genus] = {}
            final_counts[genus][species_name] = count

    # Print the final summary
    for genus, species_dict in final_counts.items():
        print("............................................................................................................")
        print(genus)
        for sp, cnt in species_dict.items():
            print(f"\tSpecies {sp} has {cnt} sequences")


# Base directory containing the sequences
base_directory = "sequences"

# Run the combined exclusion and enforcement with a minimum of 5 sequences
exclude_gapped_and_enforce_minimum(base_directory, min_sequences=5)

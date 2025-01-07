#!/usr/bin/env python3
"""

#####################################################################
*******************************WARNING*******************************
I have added the Rhino in manually, as it was not included in the original script.
sequences/gb/Hippopotamus_amphibius/Hippopotamus_amphibius.gb
sequences/fasta/Hippopotamus_amphibius/Hippopotamus_amphibius.fasta

These will need to be added to complete the intended operation/analysis
and then removed if trying to run the script as is previous to this,
unless you intend to include the Rhino in the analysis.

#####################################################################


CDS_Reference Pipeline + Post-Combination Processing

Steps:
  1) Parse GenBank files in 'gb' to extract CDS for NC_ references, embedding species in FASTA headers.
  2) Upload & align protein sequences on a remote server (MUSCLE).
  3) Add codon-based gaps to each gene’s nucleotide sequences.
  4) Combine those gapped sequences into a single file,
      'CDS_Reference_concatenated_gapped_sequences.fasta'
  5) Finally, re-process that file: for every codon (3-nt chunk),
     if there's ANY '-' or 'N' in that codon (across any record),
     we remove that codon from **all** sequences.
     This ensures no internal partial gaps remain in the final data.

The final file after Step 5 has no codons containing '-' or 'N'.
"""

import os
import shutil
import subprocess
import concurrent.futures

import paramiko
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ------------------------------------------------------------------
#  GLOBAL CONFIG
# ------------------------------------------------------------------

BASE_DIRECTORY = "sequences"                     # local base directory
REMOTE_REF_DIR = "/home/mfreeman/USCServer/CDS_Reference/"  # remote directory on server
LOCAL_REF_DIR  = os.path.join(BASE_DIRECTORY, "CDS_Reference")

SSH_SERVER = "203.101.229.234"
SSH_PORT   = 22
SSH_USER   = "mfreeman"
SSH_KEY    = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"

# Exclude ND6 from final concatenation?
EXCLUDE_ND6 = True

GENE_ORDER = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"
]


# ------------------------------------------------------------------
#  PART 1: EXTRACT REFERENCE CDS FROM GENBANK (WITH SPECIES NAME)
# ------------------------------------------------------------------

def log_message(log_file_path, message):
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + "\n")

def remove_stop_codon(sequence, accession_id, log_file_path):
    """
    Removes full/partial stop codons at the end so the final length is multiple of 3.
    Logs details.
    """
    stop_codons = ["TAA", "TAG", "AGA", "AGG"]

    # Check for full stop codon
    if len(sequence) >= 3 and sequence[-3:] in stop_codons:
        log_message(log_file_path,
                    f"{accession_id} [{sequence[-3:]}] Full stop codon removed.")
        sequence = sequence[:-3]
        if len(sequence) % 3 == 0:
            return sequence

    # If still not divisible by 3, remove partial codons
    while len(sequence) > 0 and len(sequence) % 3 != 0:
        if len(sequence) >= 2 and sequence[-2:] in ["TA", "TG", "AG"] and (len(sequence) - 2) % 3 == 0:
            log_message(log_file_path,
                        f"{accession_id} [{sequence[-2:]}] Partial stop codon (2 nt) removed.")
            sequence = sequence[:-2]
        elif len(sequence) >= 1 and sequence[-1:] in ["T", "A", "G"] and (len(sequence) - 1) % 3 == 0:
            log_message(log_file_path,
                        f"{accession_id} [{sequence[-1:]}] Partial stop codon (1 nt) removed.")
            sequence = sequence[:-1]
        else:
            break
    return sequence

def normalize_gene_name(gene_name):
    """
    Ensure standard naming, e.g. NADH1->ND1, COB->CYTB, etc.
    """
    g = gene_name.upper()
    # NADH synonyms
    if g == "NADH1":
        return "ND1"
    elif g in ["NADH2", "NADH DEHYDROGENASE SUBUNIT 2"]:
        return "ND2"
    elif g == "NADH3":
        return "ND3"
    elif g == "NADH4":
        return "ND4"
    elif g == "NADH4L":
        return "ND4L"
    elif g == "NADH5":
        return "ND5"
    elif g == "NADH6":
        return "ND6"
    # ATP
    elif g in ["ATPASE 8", "ATPASE8"]:
        return "ATP8"
    elif g in ["ATPASE 6", "ATPASE6", "ATP"]:
        return "ATP6"
    # COX
    elif g in ["COI", "COXI", "CO1"]:
        return "COX1"
    elif g in ["COII", "COXII", "CO2"]:
        return "COX2"
    elif g in ["COIII", "COXIII", "CO3"]:
        return "COX3"
    # Cytb
    elif g in ["COB", "CYTB"]:
        return "CYTB"
    return g

def normalize_combined_gene_name(gene_name):
    """
    Final standard mapping: COB->CYTB, etc.
    """
    mapping = {
        "NAD1": "ND1", "ND1": "ND1",
        "NAD2": "ND2", "ND2": "ND2",
        "NAD3": "ND3", "ND3": "ND3",
        "NAD4": "ND4", "ND4": "ND4",
        "NAD4L": "ND4L", "ND4L": "ND4L",
        "NAD5": "ND5", "ND5": "ND5",
        "NAD6": "ND6", "ND6": "ND6",
        "COI": "COX1", "CO1": "COX1", "COX1": "COX1",
        "COII": "COX2", "CO2": "COX2", "COX2": "COX2",
        "COIII": "COX3", "CO3": "COX3", "COX3": "COX3",
        "ATP6": "ATP6", "ATPASE6": "ATP6",
        "ATP8": "ATP8", "ATPASE8": "ATP8",
        "COB": "CYTB", "CYTB": "CYTB"
    }
    return mapping.get(gene_name, gene_name)

def process_reference_genbank_files(base_directory):
    """
    1) Scan base_directory/gb for .gb files.
    2) For each record with accession "NC_", extract CDS with translation.
    3) Embed species name in FASTA header.
    4) Write into:
         base_directory/CDS_Reference/CDS_nucleotide/[GENE].fasta
         base_directory/CDS_Reference/CDS_protein/[GENE].fasta
    """
    gb_directory = os.path.join(base_directory, "gb")
    ref_directory = os.path.join(base_directory, "CDS_Reference")
    nucleotide_dir = os.path.join(ref_directory, "CDS_nucleotide")
    protein_dir = os.path.join(ref_directory, "CDS_protein")
    log_file_path = os.path.join(ref_directory, "CDS_Reference_Log.txt")

    os.makedirs(nucleotide_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)

    if not os.path.exists(gb_directory):
        log_message(log_file_path, f"No directory {gb_directory} found.")
        return

    for root, dirs, files in os.walk(gb_directory):
        for file in files:
            if not file.endswith(".gb"):
                continue
            file_path = os.path.join(root, file)
            with open(file_path, 'r') as gb_handle:
                for record in SeqIO.parse(gb_handle, "gb"):
                    accession_id = record.id.split()[0]
                    if not accession_id.startswith("NC_"):
                        continue

                    organism = record.annotations.get('organism', 'Unknown').replace(' ', '_')

                    cds_features = [
                        f for f in record.features
                        if f.type == "CDS" and "translation" in f.qualifiers
                    ]
                    for feature in cds_features:
                        if "gene" in feature.qualifiers:
                            raw_gene_name = feature.qualifiers["gene"][0]
                        elif "product" in feature.qualifiers:
                            raw_gene_name = feature.qualifiers["product"][0]
                        else:
                            continue

                        gene_name = normalize_gene_name(raw_gene_name)
                        gene_name = normalize_combined_gene_name(gene_name)

                        if len(record.seq) == 0:
                            log_message(log_file_path,
                                        f"No seq data for {accession_id}, skipping.")
                            continue

                        # Nucleotide
                        nuc_seq = str(feature.extract(record.seq))
                        nuc_seq = remove_stop_codon(nuc_seq, accession_id, log_file_path)

                        # Protein
                        prot_seq = feature.qualifiers["translation"][0]

                        # Insert species name
                        header = f">{accession_id}_{organism}_{gene_name}"

                        # Write to gene-level FASTA
                        with open(os.path.join(nucleotide_dir, f"{gene_name}.fasta"), "a") as nf:
                            nf.write(f"{header}\n{nuc_seq}\n")

                        with open(os.path.join(protein_dir, f"{gene_name}.fasta"), "a") as pf:
                            pf.write(f"{header}\n{prot_seq}\n")

def step1_extract_references():
    process_reference_genbank_files(BASE_DIRECTORY)
    print("[Step 1] Reference CDS (with species name) extracted.")


# ------------------------------------------------------------------
#  PART 2: REMOTE MUSCLE ALIGNMENT
# ------------------------------------------------------------------

def run_command(command):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode

def create_ssh_client(server, port, user, key_file):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, username=user, key_filename=key_file)
    return client

def execute_command_via_ssh(client, command):
    stdin, stdout, stderr = client.exec_command(command)
    exit_code = stdout.channel.recv_exit_status()
    output = stdout.read().decode()
    error = stderr.read().decode()
    return output, error, exit_code

def ensure_reference_directories_on_server(ssh_client, remote_base_directory):
    cmd = f"mkdir -p {remote_base_directory}CDS_protein/aligned"
    _, error, exit_code = execute_command_via_ssh(ssh_client, cmd)
    if exit_code != 0:
        print(f"Error ensuring directories: {error}")
    else:
        print("Remote reference directories ensured.")

def run_muscle_alignment(client, fasta_file):
    aligned_dir = os.path.join(os.path.dirname(fasta_file), "aligned")
    base_filename = os.path.basename(fasta_file)
    afa_filename = base_filename.replace('.fasta', '.afa')
    output_path = os.path.join(aligned_dir, afa_filename).replace("\\", "/")

    # Skip if it already exists
    out, err, status = execute_command_via_ssh(client, f"test -f {output_path} && echo exists")
    if out.strip() == "exists":
        print(f"Skipping {fasta_file} (already aligned).")
        return

    # MUSCLE command (adapt if needed)
    muscle_cmd = (
        f"singularity exec /RDS/Q1233/singularity/muscle.sif "
        f"muscle -in {fasta_file} -out {output_path}"
    )

    max_attempts = 3
    for attempt in range(1, max_attempts+1):
        print(f"[MUSCLE] Aligning {fasta_file} (Attempt {attempt}/{max_attempts})...")
        out, err, status = execute_command_via_ssh(client, muscle_cmd)
        if status == 0:
            print(f"Successfully aligned -> {output_path}")
            return
        else:
            if "resource temporarily unavailable" in err.lower():
                if attempt < max_attempts:
                    print(f"Retrying {fasta_file} due to resource error...")
                    continue
                else:
                    print(f"Gave up after {max_attempts} attempts.")
            else:
                print(f"MUSCLE error: {err.strip()}")
            return

def execute_muscle_on_reference(client, remote_ref_dir):
    find_cmd = f"find {remote_ref_dir}CDS_protein -type f -name '*.fasta'"
    output, error, exit_code = execute_command_via_ssh(client, find_cmd)
    if exit_code != 0:
        print(f"Error finding .fasta: {error}")
        return

    files_list = [f.strip() for f in output.splitlines() if f.strip()]

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(run_muscle_alignment, client, f) for f in files_list]
        concurrent.futures.wait(futures)

    print("[Step 2] All reference .fasta alignments completed.")

def upload_reference_files(ssh_client, local_ref_dir, remote_ref_dir):
    sftp = ssh_client.open_sftp()
    try:
        protein_dir = os.path.join(local_ref_dir, "CDS_protein")
        if not os.path.exists(protein_dir):
            print("No local CDS_protein folder to upload.")
            return

        for file in os.listdir(protein_dir):
            if file.endswith('.fasta'):
                local_fp = os.path.join(protein_dir, file)
                remote_fp = f"{remote_ref_dir}CDS_protein/{file}"
                print(f"Uploading {local_fp} -> {remote_fp}")
                sftp.put(local_fp, remote_fp)

    except Exception as e:
        print(f"Error uploading: {str(e)}")
    finally:
        sftp.close()

def download_reference_files(ssh_client, remote_ref_dir, local_ref_dir):
    sftp = ssh_client.open_sftp()
    try:
        local_aligned_dir = os.path.join(local_ref_dir, "CDS_protein", "aligned")
        os.makedirs(local_aligned_dir, exist_ok=True)

        remote_aligned_dir = f"{remote_ref_dir}CDS_protein/aligned"
        try:
            files = sftp.listdir(remote_aligned_dir)
            afa_files = [f for f in files if f.endswith('.afa')]
            if not afa_files:
                print("No .afa files to download.")
                return
        except IOError as e:
            print(f"Failed listing remote aligned dir: {str(e)}")
            return

        for file in afa_files:
            remote_fp = f"{remote_aligned_dir}/{file}"
            local_fp = os.path.join(local_aligned_dir, file)
            tmp_fp = os.path.join(local_aligned_dir, "temp_" + file)

            if os.path.exists(local_fp):
                print(f"{file} already exists locally; skipping.")
                continue

            print(f"Downloading {remote_fp} -> {local_fp}")
            try:
                sftp.get(remote_fp, tmp_fp)
                if os.path.getsize(tmp_fp) > 0:
                    shutil.move(tmp_fp, local_fp)
                else:
                    print("Downloaded file is empty, removing.")
                    os.remove(tmp_fp)
            except Exception as e:
                print(f"Error downloading {file}: {str(e)}")
                if os.path.exists(tmp_fp):
                    os.remove(tmp_fp)

    finally:
        sftp.close()

def step2_remote_alignment():
    client = create_ssh_client(SSH_SERVER, SSH_PORT, SSH_USER, SSH_KEY)
    ensure_reference_directories_on_server(client, REMOTE_REF_DIR)
    upload_reference_files(client, LOCAL_REF_DIR, REMOTE_REF_DIR)
    execute_muscle_on_reference(client, REMOTE_REF_DIR)
    download_reference_files(client, REMOTE_REF_DIR, LOCAL_REF_DIR)
    client.close()
    print("[Step 2] Remote alignment & download complete.")


# ------------------------------------------------------------------
#  PART 3: ADDING CODON-BASED GAPS
# ------------------------------------------------------------------

def add_codon_gaps_to_references(base_dir):
    """
    For each .afa in CDS_protein/aligned, find the matching .fasta in CDS_nucleotide,
    insert '---' for each gap in the protein alignment, then write to CDS_nucleotide_gapped.
    """
    aligned_dir = os.path.join(base_dir, "CDS_protein", "aligned")
    nuc_dir     = os.path.join(base_dir, "CDS_nucleotide")
    gapped_dir  = os.path.join(base_dir, "CDS_nucleotide_gapped")

    if not os.path.exists(aligned_dir):
        print("No aligned dir found.")
        return

    os.makedirs(gapped_dir, exist_ok=True)

    for filename in os.listdir(aligned_dir):
        if not filename.endswith(".afa"):
            continue

        gene_name = os.path.splitext(filename)[0]  # e.g. ND1
        protein_path = os.path.join(aligned_dir, filename)
        nucleotide_path = os.path.join(nuc_dir, f"{gene_name}.fasta")

        if not os.path.exists(nucleotide_path):
            print(f"No matching nucleotide file for {gene_name}")
            continue

        prot_recs = list(SeqIO.parse(protein_path, "fasta"))
        nuc_recs  = list(SeqIO.parse(nucleotide_path, "fasta"))

        # ID-based matching
        nuc_dict = {r.id: r for r in nuc_recs}
        mod_recs = []

        for prot_rec in prot_recs:
            nuc_rec = nuc_dict.get(prot_rec.id)
            if not nuc_rec:
                print(f"No matching nt record for {prot_rec.id}")
                continue

            gene_seq = str(nuc_rec.seq)
            new_seq = []
            idx = 0

            for aa in prot_rec.seq:
                if aa == '-':
                    # Insert codon gap
                    new_seq.append('---')
                else:
                    codon = gene_seq[idx:idx+3]
                    new_seq.append(codon)
                    idx += 3

            mod_seq = "".join(new_seq)
            new_rec = SeqRecord(
                Seq(mod_seq),
                id=nuc_rec.id,
                description=nuc_rec.description
            )
            mod_recs.append(new_rec)

        outpath = os.path.join(gapped_dir, f"{gene_name}.fasta")
        SeqIO.write(mod_recs, outpath, 'fasta')
        print(f"Gapped {gene_name} -> {outpath}")

def step3_add_codon_gaps():
    add_codon_gaps_to_references(LOCAL_REF_DIR)
    print("[Step 3] Codon-based gaps applied to reference nucleotides.")


# ------------------------------------------------------------------
#  PART 4: COMBINE GAPPED SEQUENCES INTO ONE FILE
# ------------------------------------------------------------------

def get_filtered_gene_order():
    if EXCLUDE_ND6:
        return [g for g in GENE_ORDER if g != "ND6"]
    else:
        return GENE_ORDER

def step4_combine():
    """
    1. For each gene in get_filtered_gene_order, read the gapped FASTA in
       CDS_nucleotide_gapped,
    2. Remove any trailing _GENE from the record ID,
    3. Concatenate everything for each unique accession,
    4. Write single combined file: 'CDS_Reference_concatenated_gapped_sequences.fasta'
    """
    gapped_dir   = os.path.join(LOCAL_REF_DIR, "CDS_nucleotide_gapped")
    output_file  = os.path.join(gapped_dir, "CDS_Reference_concatenated_gapped_sequences.fasta")

    genes = get_filtered_gene_order()
    combined_dict = {}  # {accession_without_gene: concatenated_seq}

    for gene in genes:
        gene_file = os.path.join(gapped_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} not found; skipping.")
            continue

        with open(gene_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                accession = record.id
                # Strip _GENE from the end if present
                for gnm in GENE_ORDER:
                    suff = f"_{gnm}"
                    if accession.endswith(suff):
                        accession = accession[: -len(suff)]
                        break

                if accession not in combined_dict:
                    combined_dict[accession] = ""
                combined_dict[accession] += str(record.seq)

    # Write the combined file (unprocessed)
    with open(output_file, "w") as out:
        for acc, seq in combined_dict.items():
            out.write(f">{acc}\n{seq}\n")

    print("[Step 4] Created combined file ->", output_file)


# ------------------------------------------------------------------
#  PART 5: REMOVE ANY CODONS WITH '-' OR 'N'
# ------------------------------------------------------------------

def step5_remove_codons_with_gaps_or_N():
    """
    Reads the newly combined file
      'CDS_Reference_concatenated_gapped_sequences.fasta'
    and removes any codon (3-nt chunk) that has at least one '-' or 'N'.
    Then rewrites the same file with the final cleaned sequences.
    """
    gapped_dir  = os.path.join(LOCAL_REF_DIR, "CDS_nucleotide_gapped")
    combined_fp = os.path.join(gapped_dir, "CDS_Reference_concatenated_gapped_sequences.fasta")

    if not os.path.exists(combined_fp):
        print(f"No combined file found at {combined_fp}")
        return

    # Read all sequences into a dict
    seq_dict = {}
    for record in SeqIO.parse(combined_fp, "fasta"):
        seq_dict[record.id] = str(record.seq)

    if not seq_dict:
        print("No sequences found in combined file. Nothing to process.")
        return

    # 1. Determine the shortest length among all sequences
    min_len = min(len(s) for s in seq_dict.values())
    # 2. Trim all sequences to that length
    for acc in seq_dict:
        seq_dict[acc] = seq_dict[acc][:min_len]

    # 3. Ensure multiple of 3 (trim off extras if needed)
    remainder = min_len % 3
    if remainder != 0:
        min_len -= remainder
        for acc in seq_dict:
            seq_dict[acc] = seq_dict[acc][:min_len]

    if min_len == 0:
        # If we lost everything, just overwrite with empty
        with open(combined_fp, "w") as out:
            for acc in seq_dict:
                out.write(f">{acc}\n\n")
        print("[Step 5] All sequences ended up 0-length. File overwritten with empties.")
        return

    # 4. Now remove any codon that contains '-' or 'N'
    num_codons = min_len // 3
    codon_mask = [True] * num_codons

    # Build a mask across all sequences
    for i in range(num_codons):
        start = i * 3
        end   = start + 3
        # If ANY sequence has '-' or 'N' in that codon, mask = False
        for seq in seq_dict.values():
            codon = seq[start:end]
            if '-' in codon or 'N' in codon:
                codon_mask[i] = False
                break

    # 5. Rebuild sequences, skipping masked codons
    for acc, fullseq in seq_dict.items():
        parts = []
        for i in range(num_codons):
            if codon_mask[i]:
                start = i*3
                end   = start+3
                parts.append(fullseq[start:end])
        seq_dict[acc] = "".join(parts)

    # 6. Optionally, ensure final length is consistent
    final_min_len = min(len(s) for s in seq_dict.values()) if seq_dict else 0
    for acc in seq_dict:
        seq_dict[acc] = seq_dict[acc][:final_min_len]

    # Overwrite the combined file with final cleaned sequences
    with open(combined_fp, "w") as out:
        for acc, sequence in seq_dict.items():
            out.write(f">{acc}\n{sequence}\n")

    print("[Step 5] Removed any codons containing '-' or 'N'. Final file overwritten:")

    # Show final stats
    for acc in seq_dict:
        print(f" {acc} => length {len(seq_dict[acc])}")


# ------------------------------------------------------------------
#  MAIN
# ------------------------------------------------------------------

def main():
    # Step 1: Extract from GenBank
    step1_extract_references()

    # Step 2: Remote MUSCLE alignment
    step2_remote_alignment()

    # Step 3: Add codon-based gaps
    step3_add_codon_gaps()

    # Step 4: Combine all genes into one file
    step4_combine()

    # Step 5: Remove any codons containing '-' or 'N'
    step5_remove_codons_with_gaps_or_N()

if __name__ == "__main__":
    main()

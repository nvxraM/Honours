#!/usr/bin/env python3
"""
CDS_Reference Pipeline (Unified Ingroup + Outgroup, ONLY NC_ references)

This script processes *all* GenBank files from:
  - sequences/gb/...
  - sequences/outgroups/gb/...      <--- NOTE the "outgroups" folder name
BUT it only extracts sequences where the record.id starts with "NC_".

Any record that does not start with NC_ is skipped.

The pipeline:
  1) Extracts CDS from these NC_ references -> local (CDS_nucleotide, CDS_protein)
  2) Uploads & aligns protein sequences on remote server (MUSCLE)
  3) Inserts codon-based gaps in the nucleotide sequences
  4) Combines them into one big FASTA
  5) Removes any codons that contain '-' or 'N' in any record
"""

import os
import shutil
import subprocess
import concurrent.futures

import paramiko
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# -----------------------------
#  GLOBAL SETTINGS
# -----------------------------

# Now referencing "outgroups" instead of "outgroup"
GB_DIRECTORIES = [
    os.path.join("sequences", "gb"),
    os.path.join("sequences", "outgroups", "gb")
]

BASE_DIRECTORY = "sequences"
REMOTE_REF_DIR = "/home/mfreeman/USCServer/CDS_Reference/"
LOCAL_REF_DIR  = os.path.join(BASE_DIRECTORY, "CDS_Reference")

SSH_SERVER = "203.101.229.234"
SSH_PORT   = 22
SSH_USER   = "mfreeman"
SSH_KEY    = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"

# Whether to exclude ND6 from the final combined alignment
EXCLUDE_ND6 = True

GENE_ORDER = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"
]


# -----------------------------
#  UTILITY FUNCTIONS
# -----------------------------

def log_message(log_file_path, message):
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + "\n")


def remove_stop_codon(sequence, accession_id, log_file_path):
    """
    Removes full/partial stop codons at the end so the final length is multiple of 3.
    Logs details to a file.
    """
    stop_codons = ["TAA", "TAG", "AGA", "AGG"]

    # Check for a full (3-nt) stop codon
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


def run_command(command):
    """
    Simple local shell command runner.
    """
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode


# ------------------------------------------------------------
# STEP 1: EXTRACT ALL REFERENCE CDS (ONLY NC_)
# ------------------------------------------------------------

def process_all_genbank_files():
    """
    1) Walk both GB_DIRECTORIES for .gb files.
    2) For each record, if record.id starts with NC_, then extract CDS.
       Otherwise skip it (not a reference).
    3) Write to sequences/CDS_Reference/CDS_nucleotide/ and /CDS_protein/
    """
    ref_dir = os.path.join(BASE_DIRECTORY, "CDS_Reference")
    nuc_dir = os.path.join(ref_dir, "CDS_nucleotide")
    prot_dir = os.path.join(ref_dir, "CDS_protein")
    log_file = os.path.join(ref_dir, "CDS_Reference_Log.txt")

    os.makedirs(nuc_dir, exist_ok=True)
    os.makedirs(prot_dir, exist_ok=True)

    for gb_dir in GB_DIRECTORIES:
        if not os.path.exists(gb_dir):
            continue  # skip missing directories

        for root, dirs, files in os.walk(gb_dir):
            for file in files:
                if not file.endswith(".gb"):
                    continue
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as gb_handle:
                    for record in SeqIO.parse(gb_handle, "gb"):
                        accession_id = record.id.split()[0]
                        # Only use references that start with NC_
                        if not accession_id.startswith("NC_"):
                            continue

                        organism = record.annotations.get('organism', 'Unknown')
                        organism = organism.replace(' ', '_')

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
                                log_message(log_file, f"No seq data for {accession_id}, skipping.")
                                continue

                            # Extract nucleotide from feature
                            nuc_seq = str(feature.extract(record.seq))
                            nuc_seq = remove_stop_codon(nuc_seq, accession_id, log_file)

                            # Protein
                            prot_seq = feature.qualifiers["translation"][0]

                            # Build header with species + gene
                            header = f">{accession_id}_{organism}_{gene_name}"

                            # Append to gene FASTA
                            with open(os.path.join(nuc_dir, f"{gene_name}.fasta"), "a") as nf:
                                nf.write(f"{header}\n{nuc_seq}\n")

                            with open(os.path.join(prot_dir, f"{gene_name}.fasta"), "a") as pf:
                                pf.write(f"{header}\n{prot_seq}\n")


def step1_extract_references():
    process_all_genbank_files()
    print("[Step 1] Extracted CDS from NC_ references (both ingroup + outgroup).")


# ------------------------------------------------------------
# STEP 2: REMOTE MUSCLE ALIGNMENT
# ------------------------------------------------------------

def ensure_remote_directories(ssh_client, remote_base_directory):
    cmd = f"mkdir -p {remote_base_directory}CDS_protein/aligned"
    _, error, exit_code = execute_command_via_ssh(ssh_client, cmd)
    if exit_code != 0:
        print(f"Error ensuring remote directories: {error}")
    else:
        print("[Remote] reference directories ensured.")


def upload_protein_fastas(ssh_client, local_ref_dir, remote_ref_dir):
    """
    Upload all protein .fasta from local to remote.
    """
    sftp = ssh_client.open_sftp()
    try:
        protein_dir = os.path.join(local_ref_dir, "CDS_protein")
        if not os.path.exists(protein_dir):
            print("[Upload] No local CDS_protein folder found.")
            return

        for file in os.listdir(protein_dir):
            if file.endswith('.fasta'):
                local_fp = os.path.join(protein_dir, file)
                remote_fp = f"{remote_ref_dir}CDS_protein/{file}"
                print(f"[Upload] {local_fp} -> {remote_fp}")
                sftp.put(local_fp, remote_fp)
    except Exception as e:
        print(f"Error uploading: {str(e)}")
    finally:
        sftp.close()


def run_muscle_alignment(ssh_client, fasta_file):
    """
    For a given .fasta on the remote server, produce an aligned .afa using MUSCLE.
    """
    aligned_dir = os.path.join(os.path.dirname(fasta_file), "aligned")
    base_filename = os.path.basename(fasta_file)
    afa_filename = base_filename.replace('.fasta', '.afa')
    output_path = os.path.join(aligned_dir, afa_filename).replace("\\", "/")

    # Skip if .afa already exists
    out, err, status = execute_command_via_ssh(ssh_client, f"test -f {output_path} && echo exists")
    if out.strip() == "exists":
        print(f"[MUSCLE] Skipping {fasta_file}, already aligned.")
        return

    muscle_cmd = (
        f"singularity exec /RDS/Q1233/singularity/muscle.sif "
        f"muscle -in {fasta_file} -out {output_path}"
    )

    max_attempts = 3
    for attempt in range(1, max_attempts+1):
        print(f"[MUSCLE] Aligning {fasta_file} (Attempt {attempt}/{max_attempts})...")
        out, err, status = execute_command_via_ssh(ssh_client, muscle_cmd)
        if status == 0:
            print(f"[MUSCLE] Success -> {output_path}")
            return
        else:
            if "resource temporarily unavailable" in err.lower():
                if attempt < max_attempts:
                    print("[MUSCLE] Resource error, retrying...")
                    continue
                else:
                    print("[MUSCLE] Gave up after 3 attempts.")
            else:
                print(f"[MUSCLE] Error: {err.strip()}")
            return


def execute_muscle_on_remote(ssh_client, remote_ref_dir):
    """
    1) Find all .fasta in remote CDS_protein
    2) Align each with MUSCLE in parallel
    """
    find_cmd = f"find {remote_ref_dir}CDS_protein -type f -name '*.fasta'"
    output, error, exit_code = execute_command_via_ssh(ssh_client, find_cmd)
    if exit_code != 0:
        print(f"[Remote] Error finding .fasta: {error}")
        return

    files_list = [line.strip() for line in output.splitlines() if line.strip()]

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        for f in files_list:
            futures.append(executor.submit(run_muscle_alignment, ssh_client, f))
        concurrent.futures.wait(futures)

    print("[Step 2] All protein FASTAs aligned on remote server.")


def download_aligned_fastas(ssh_client, remote_ref_dir, local_ref_dir):
    """
    Download the .afa alignments back to local.
    """
    sftp = ssh_client.open_sftp()
    try:
        local_aligned_dir = os.path.join(local_ref_dir, "CDS_protein", "aligned")
        os.makedirs(local_aligned_dir, exist_ok=True)

        remote_aligned_dir = f"{remote_ref_dir}CDS_protein/aligned"
        try:
            files = sftp.listdir(remote_aligned_dir)
            afa_files = [f for f in files if f.endswith('.afa')]
            if not afa_files:
                print("[Download] No .afa files to download.")
                return
        except IOError as e:
            print(f"[Download] Failed listing remote aligned dir: {str(e)}")
            return

        for file in afa_files:
            remote_fp = f"{remote_aligned_dir}/{file}"
            local_fp = os.path.join(local_aligned_dir, file)
            tmp_fp = os.path.join(local_aligned_dir, "temp_" + file)

            if os.path.exists(local_fp):
                print(f"[Download] {file} already exists locally; skipping.")
                continue

            print(f"[Download] {remote_fp} -> {local_fp}")
            try:
                sftp.get(remote_fp, tmp_fp)
                if os.path.getsize(tmp_fp) > 0:
                    shutil.move(tmp_fp, local_fp)
                else:
                    print("[Download] File is empty, removing.")
                    os.remove(tmp_fp)
            except Exception as e:
                print(f"[Download] Error downloading {file}: {str(e)}")
                if os.path.exists(tmp_fp):
                    os.remove(tmp_fp)

    finally:
        sftp.close()


def step2_remote_alignment():
    client = create_ssh_client(SSH_SERVER, SSH_PORT, SSH_USER, SSH_KEY)
    ensure_remote_directories(client, REMOTE_REF_DIR)
    upload_protein_fastas(client, LOCAL_REF_DIR, REMOTE_REF_DIR)
    execute_muscle_on_remote(client, REMOTE_REF_DIR)
    download_aligned_fastas(client, REMOTE_REF_DIR, LOCAL_REF_DIR)
    client.close()
    print("[Step 2] Remote MUSCLE alignment + download complete.")


# ------------------------------------------------------------
# STEP 3: ADD CODON-BASED GAPS TO NUCLEOTIDES
# ------------------------------------------------------------

def add_codon_gaps_to_nucleotides(base_dir):
    """
    For each .afa in CDS_protein/aligned, find matching .fasta in CDS_nucleotide,
    insert '---' for each gap in the protein alignment, and write to CDS_nucleotide_gapped.
    """
    aligned_dir = os.path.join(base_dir, "CDS_protein", "aligned")
    nuc_dir     = os.path.join(base_dir, "CDS_nucleotide")
    gapped_dir  = os.path.join(base_dir, "CDS_nucleotide_gapped")

    if not os.path.exists(aligned_dir):
        print("[Step 3] No aligned directory found. Skipping codon-gap step.")
        return

    os.makedirs(gapped_dir, exist_ok=True)

    for filename in os.listdir(aligned_dir):
        if not filename.endswith(".afa"):
            continue

        gene_name = os.path.splitext(filename)[0]
        protein_path = os.path.join(aligned_dir, filename)
        nucleotide_path = os.path.join(nuc_dir, f"{gene_name}.fasta")

        if not os.path.exists(nucleotide_path):
            print(f"[Step 3] No matching nucleotide file for {gene_name}, skipping.")
            continue

        prot_recs = list(SeqIO.parse(protein_path, "fasta"))
        nuc_recs  = list(SeqIO.parse(nucleotide_path, "fasta"))

        # Build dict for quick ID lookup
        nuc_dict = {r.id: r for r in nuc_recs}
        mod_recs = []

        for prot_rec in prot_recs:
            if prot_rec.id not in nuc_dict:
                print(f"[Step 3] No matching NT record for {prot_rec.id}, skipping.")
                continue

            nuc_rec = nuc_dict[prot_rec.id]
            nt_seq_str = str(nuc_rec.seq)
            new_seq = []
            idx = 0

            for aa in prot_rec.seq:
                if aa == '-':
                    # Insert a 3-nt gap
                    new_seq.append('---')
                else:
                    codon = nt_seq_str[idx: idx+3]
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
        SeqIO.write(mod_recs, outpath, "fasta")
        print(f"[Step 3] Gapped {gene_name} -> {outpath}")


def step3_add_codon_gaps():
    add_codon_gaps_to_nucleotides(LOCAL_REF_DIR)
    print("[Step 3] Codon-based gaps applied to nucleotides.")


# ------------------------------------------------------------
# STEP 4: COMBINE GAPPED SEQUENCES INTO A SINGLE FILE
# ------------------------------------------------------------

def get_filtered_gene_order():
    """
    Return GENE_ORDER with ND6 optionally excluded.
    """
    if EXCLUDE_ND6:
        return [g for g in GENE_ORDER if g != "ND6"]
    else:
        return GENE_ORDER


def step4_combine():
    """
    Read each gene's gapped FASTA (CDS_nucleotide_gapped/[GENE].fasta),
    concatenate all genes for each sample, produce one big FASTA:
    CDS_Reference_concatenated_gapped_sequences.fasta
    """
    gapped_dir  = os.path.join(LOCAL_REF_DIR, "CDS_nucleotide_gapped")
    output_file = os.path.join(gapped_dir, "CDS_Reference_concatenated_gapped_sequences.fasta")

    genes = get_filtered_gene_order()
    combined_dict = {}  # { sample_id: combined_seq }

    for gene in genes:
        gene_fasta = os.path.join(gapped_dir, f"{gene}.fasta")
        if not os.path.exists(gene_fasta):
            print(f"[Step 4] Missing gapped file for {gene}, skipping.")
            continue

        for record in SeqIO.parse(gene_fasta, "fasta"):
            sample_id = record.id
            # Also strip out trailing _GENE if present
            for gnm in GENE_ORDER:
                suffix = f"_{gnm}"
                if sample_id.endswith(suffix):
                    sample_id = sample_id[:-len(suffix)]
                    break

            seq_str = str(record.seq)
            if sample_id not in combined_dict:
                combined_dict[sample_id] = ""
            combined_dict[sample_id] += seq_str

    with open(output_file, "w") as out:
        for sample_id, seq in combined_dict.items():
            out.write(f">{sample_id}\n{seq}\n")

    print("[Step 4] Combined all gapped genes into:", output_file)


# ------------------------------------------------------------
# STEP 5: REMOVE ANY CODONS WITH '-' OR 'N'
# ------------------------------------------------------------

def step5_remove_codons_with_gaps_or_N():
    """
    We open the single combined FASTA, and for each codon (3-nt chunk):
     - If ANY sample has a '-' or 'N' in that codon, remove that codon from ALL samples.
    """
    gapped_dir   = os.path.join(LOCAL_REF_DIR, "CDS_nucleotide_gapped")
    combined_fp  = os.path.join(gapped_dir, "CDS_Reference_concatenated_gapped_sequences.fasta")

    if not os.path.exists(combined_fp):
        print("[Step 5] No combined file found, skipping.")
        return

    seq_dict = {}
    for record in SeqIO.parse(combined_fp, "fasta"):
        seq_dict[record.id] = str(record.seq)

    if not seq_dict:
        print("[Step 5] Combined file is empty. Nothing to clean.")
        return

    # 1. Find the minimum sequence length
    min_len = min(len(s) for s in seq_dict.values())

    # 2. Truncate all sequences to min_len
    for k in seq_dict:
        seq_dict[k] = seq_dict[k][:min_len]

    # 3. Make sure length is multiple of 3
    remainder = min_len % 3
    if remainder != 0:
        min_len -= remainder
        for k in seq_dict:
            seq_dict[k] = seq_dict[k][:min_len]

    if min_len == 0:
        # everything is length zero
        with open(combined_fp, "w") as out:
            for k in seq_dict:
                out.write(f">{k}\n\n")
        print("[Step 5] All sequences reduced to length 0.")
        return

    num_codons = min_len // 3
    # 4. Build a mask to indicate which codons are "clean" (no '-' or 'N' in any sample)
    codon_mask = [True] * num_codons

    for i in range(num_codons):
        start = i * 3
        end   = start + 3
        for seq in seq_dict.values():
            codon = seq[start:end]
            if '-' in codon or 'N' in codon:
                codon_mask[i] = False
                break

    # 5. Rebuild each sequence with only the codons that pass
    for k, fullseq in seq_dict.items():
        new_codons = []
        for i in range(num_codons):
            if codon_mask[i]:
                start = i * 3
                end   = start + 3
                new_codons.append(fullseq[start:end])
        seq_dict[k] = "".join(new_codons)

    # 6. Optionally ensure final min length is consistent
    final_len = min(len(s) for s in seq_dict.values()) if seq_dict else 0
    for k in seq_dict:
        seq_dict[k] = seq_dict[k][:final_len]

    # Overwrite the combined file
    with open(combined_fp, "w") as out:
        for k, s in seq_dict.items():
            out.write(f">{k}\n{s}\n")

    print("[Step 5] Removed codons with '-' or 'N'. Final alignment saved to:")
    print(" ", combined_fp)

    # Quick stats
    for k in seq_dict:
        print(f"  {k} => length {len(seq_dict[k])}")


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------

def main():
    # Step 1: Extract from .gb files (ONLY if accession_id.startswith("NC_"))
    step1_extract_references()

    # Step 2: Align protein sequences on remote server
    step2_remote_alignment()

    # Step 3: Insert codon-based gaps into nucleotide sequences
    step3_add_codon_gaps()

    # Step 4: Combine all gapped sequences into one multi-gene FASTA
    step4_combine()

    # Step 5: Remove codons that contain '-' or 'N' in any sample
    step5_remove_codons_with_gaps_or_N()

if __name__ == "__main__":
    main()

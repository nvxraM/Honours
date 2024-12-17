#!/usr/bin/env python3
"""
Revised Script:

This script:
1. Processes GenBank files to extract reference CDS_Genus (from "NC_" records) into one common directory:
   sequences/CDS_Reference/All_References/...

2. From these references, it combines all protein sequences by gene into:
   sequences/CDS_Reference/All_References/CDS_Reference_protein/

3. Aligns them using MUSCLE on a remote server and downloads the aligned files.

4. Inserts codon-based gaps into the nucleotide sequences based on protein alignment, storing gapped
   sequences in:
   sequences/CDS_Reference/All_References/CDS_Reference_nucleotide_gapped/

5. Concatenates all gapped gene sequences into one large FASTA:
   sequences/CDS_Reference/All_MEGA_Grouped_Individually/All_concatenated_gapped_sequences.fasta

6. Creates a single .grp file for these sequences, grouped by Genus.

7. Runs MEGA analyses once on the combined dataset.

8. Aligns all combined sequences together locally with MUSCLE and removes columns with gaps or Ns at the codon level,
   leaving a perfectly aligned, codon-based, gap-free, and N-free final dataset.
"""

import os
import re
import shutil
import subprocess
import concurrent.futures
import paramiko
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

#########################################
# Constants and Configuration Parameters
#########################################

GENE_ORDER = [
    "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3",
    "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6"
]

EXCLUDE_ND6 = True


#########################################
# Utility and Helper Functions
#########################################

def log_message(log_file_path, message):
    with open(log_file_path, 'a') as log_file:
        log_file.write(message + "\n")


def normalize_gene_name(gene_name):
    gene_name = gene_name.upper()
    conversion_map = {
        "NADH1": "ND1", "NADH2": "ND2", "NADH3": "ND3", "NADH4": "ND4",
        "NADH4L": "ND4L", "NADH5": "ND5", "NADH6": "ND6",
        "ATPASE 8": "ATP8", "ATPASE8": "ATP8",
        "ATPASE 6": "ATP6", "ATPASE6": "ATP6", "ATP": "ATP6",
        "COI": "COX1", "COXI": "COX1", "CO1": "COX1",
        "COII": "COX2", "COXII": "COX2", "CO2": "COX2",
        "COIII": "COX3", "COXIII": "COX3", "CO3": "COX3"
    }
    return conversion_map.get(gene_name, gene_name)


def normalize_combined_gene_name(gene_name):
    conversion_map = {
        "NAD1": "ND1", "NAD2": "ND2", "NAD3": "ND3", "NAD4": "ND4",
        "NAD4L": "ND4L", "NAD5": "ND5", "NAD6": "ND6",
        "CO1": "COX1", "CO2": "COX2", "CO3": "COX3",
        "ATPASE6": "ATP6", "ATPASE8": "ATP8", "COB": "CYTB"
    }
    return conversion_map.get(gene_name, gene_name)


def remove_stop_codon(sequence, accession_id, species_name, log_file_path):
    stop_codons = ["TAA", "TAG", "AGA", "AGG"]
    if len(sequence) >= 3 and sequence[-3:] in stop_codons:
        log_message(log_file_path, f"{species_name} {accession_id} [{sequence[-3:]}] Full stop codon removed.")
        sequence = sequence[:-3]
        if len(sequence) % 3 == 0:
            return sequence

    while len(sequence) > 0 and (len(sequence) % 3 != 0):
        last_base = sequence[-1:]
        last_two = sequence[-2:]
        if len(sequence) >= 2 and last_two in ["TA", "TG", "AG"] and (len(sequence) - 2) % 3 == 0:
            log_message(log_file_path, f"{species_name} {accession_id} [{last_two}] Partial stop codon removed.")
            sequence = sequence[:-2]
        elif len(sequence) >= 1 and last_base in ["T", "A", "G"] and (len(sequence) - 1) % 3 == 0:
            log_message(log_file_path, f"{species_name} {accession_id} [{last_base}] Partial stop codon removed.")
            sequence = sequence[:-1]
        else:
            break

    return sequence


def get_filtered_gene_order():
    if EXCLUDE_ND6:
        return [g for g in GENE_ORDER if g != "ND6"]
    return GENE_ORDER


#########################################
# Step 1: Process GenBank Files into All_References
#########################################

def find_fasta_file(base_directory, species_folder):
    from Bio import SeqIO
    species_fasta_dir = os.path.join(base_directory, "fasta", species_folder)
    if os.path.exists(species_fasta_dir):
        for fn in os.listdir(species_fasta_dir):
            if fn.endswith('.fasta'):
                return os.path.join(species_fasta_dir, fn)

    fasta_base_dir = os.path.join(base_directory, "fasta")
    for root, _, files in os.walk(fasta_base_dir):
        for fn in files:
            if fn.endswith('.fasta'):
                fasta_path = os.path.join(root, fn)
                with open(fasta_path, 'r') as fh:
                    for record in SeqIO.parse(fh, "fasta"):
                        if species_folder in record.id:
                            return fasta_path
    return None


def extract_and_save(records, gene_directory, protein_directory, fasta_records, species_name, log_file_path):
    """
    Extract CDS_Genus from records starting with 'NC_', save nucleotide and protein sequences.
    Store them all in the All_References directories rather than genus-level directories.
    """
    for record in records:
        if not record.id.startswith("NC_"):
            continue
        cds_features = [f for f in record.features if f.type == "CDS_Genus" and "translation" in f.qualifiers]

        for feature in cds_features:
            if "gene" in feature.qualifiers:
                gene_name = normalize_gene_name(feature.qualifiers['gene'][0])
            elif "product" in feature.qualifiers:
                product_name = feature.qualifiers['product'][0].lower()
                if product_name.startswith("nadh dehydrogenase subunit"):
                    subunit_number = product_name.split()[-1]
                    if subunit_number.isdigit():
                        gene_name = f"ND{subunit_number}"
                    else:
                        continue
                else:
                    continue
            else:
                continue

            combined_gene_name = normalize_combined_gene_name(gene_name)
            accession_id = record.id.split()[0]

            if len(record.seq) == 0:
                if accession_id in fasta_records:
                    seq_to_use = fasta_records[accession_id].seq
                    source = 'FASTA record'
                else:
                    log_message(log_file_path, f"Sequence not found for {accession_id}, gene {gene_name}.")
                    continue
            else:
                seq_to_use = record.seq
                source = 'GenBank record'

            nucleotide_seq = str(feature.extract(seq_to_use))
            log_message(log_file_path, f"Extracted {accession_id}, gene {gene_name} from {source}.")
            nucleotide_seq = remove_stop_codon(nucleotide_seq, accession_id,
                                               record.annotations.get('organism', 'Unknown'), log_file_path)

            protein_seq = feature.qualifiers['translation'][0]
            # Record ID: include species_name to keep them unique
            species_gene_header = f">{record.id}_{species_name}"

            # Append nucleotide to gene file
            gene_file = os.path.join(gene_directory, f"{combined_gene_name}.fasta")
            with open(gene_file, 'a') as gf:
                gf.write(f"{species_gene_header}\n{nucleotide_seq}\n")

            # Append protein to protein file
            protein_file = os.path.join(protein_directory, f"{combined_gene_name}.fasta")
            with open(protein_file, 'a') as pf:
                pf.write(f"{species_gene_header}\n{protein_seq}\n")


def recheck_sequence_lengths(directory, log_file_path):
    """
    Recheck sequences in All_References/CDS_Reference_nucleotide to ensure no trailing stops.
    """
    for gene_file in os.listdir(directory):
        if gene_file.endswith('.fasta'):
            file_path = os.path.join(directory, gene_file)
            adjusted = []
            with open(file_path, 'r') as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    species_name = rec.id.split('_', 1)[1] if '_' in rec.id else rec.id
                    cleaned_seq = remove_stop_codon(str(rec.seq), rec.id, species_name, log_file_path)
                    adjusted.append((rec.id, cleaned_seq))
            with open(file_path, 'w') as fh:
                for rid, seq in adjusted:
                    fh.write(f">{rid}\n{seq}\n")


def process_genbank_files(base_directory):
    from Bio import SeqIO
    input_directory = os.path.join(base_directory, "gb")

    # Output directories for all references
    all_ref_dir = os.path.join(base_directory, "CDS_Reference", "All_References")
    gene_dir = os.path.join(all_ref_dir, 'CDS_Reference_nucleotide')
    protein_dir = os.path.join(all_ref_dir, 'CDS_Reference_protein')
    os.makedirs(gene_dir, exist_ok=True)
    os.makedirs(protein_dir, exist_ok=True)
    log_file_path = os.path.join(base_directory, "CDS_Reference_Acquisition_Log.txt")

    if not os.path.exists(input_directory):
        log_message(log_file_path, f"Input directory {input_directory} does not exist.")
        return

    # Process each species
    for species_folder in os.listdir(input_directory):
        species_path = os.path.join(input_directory, species_folder)
        if os.path.isdir(species_path):
            fasta_file_path = find_fasta_file(base_directory, species_folder)
            if not fasta_file_path:
                log_message(log_file_path, f"FASTA file for {species_folder} not found.")
                continue

            with open(fasta_file_path, 'r') as f:
                fasta_records = {r.id: r for r in SeqIO.parse(f, "fasta")}

            for genbank_file in os.listdir(species_path):
                if genbank_file.endswith('.gb'):
                    gb_file_path = os.path.join(species_path, genbank_file)
                    with open(gb_file_path, 'r') as gb_handle:
                        for record in SeqIO.parse(gb_handle, "gb"):
                            extract_and_save([record], gene_dir, protein_dir, fasta_records, species_folder,
                                             log_file_path)
            log_message(log_file_path, f"Processed CDS_Reference for species: {species_folder}")

    # Recheck stop codons removal
    recheck_sequence_lengths(gene_dir, log_file_path)


#########################################
# Step 2: Combine All Genes (Now Already in All_References)
#########################################

def create_combined_protein_files(base_directory):
    # Already done by extraction. No further action needed.
    pass


def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
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


def ensure_remote_combined_dir(ssh_client, remote_directory):
    cmd = f"mkdir -p {remote_directory}CDS_Reference_protein/aligned"
    _, error, exit_code = execute_command_via_ssh(ssh_client, cmd)
    if exit_code != 0:
        print(f"Error ensuring remote combined directory: {error}")
    else:
        print("Remote combined directory ensured.")


def upload_combined_files(ssh_client, local_directory, remote_directory):
    sftp = ssh_client.open_sftp()
    ref_protein_dir = os.path.join(local_directory, "CDS_Reference", "All_References", "CDS_Reference_protein")
    if not os.path.exists(ref_protein_dir):
        print("No protein directory found.")
        sftp.close()
        return

    remote_combined_dir = f"{remote_directory}CDS_Reference_protein"
    _, dir_err, dir_status = execute_command_via_ssh(ssh_client, f'mkdir -p {remote_combined_dir}/aligned')
    if dir_status != 0:
        print(f"Error ensuring remote directory: {dir_err.strip()}")

    for gene_name in get_filtered_gene_order():
        local_file_path = os.path.join(ref_protein_dir, f"{gene_name}.fasta")
        if os.path.exists(local_file_path):
            remote_file_path = os.path.join(remote_combined_dir, f"{gene_name}.fasta").replace("\\", "/")
            out, err, st = execute_command_via_ssh(ssh_client, f"test -f {remote_file_path} && echo exists")
            if out.strip() == 'exists':
                print(f"{remote_file_path} already exists on server. Skipping upload.")
            else:
                print(f"Uploading {local_file_path} to {remote_file_path}")
                sftp.put(local_file_path, remote_file_path)
        else:
            print(f"No sequences for {gene_name}. Skipping.")
    sftp.close()


def run_muscle_alignment_on_combined(client, remote_directory):
    remote_combined_dir = f"{remote_directory}CDS_Reference_protein"
    out, err, status = execute_command_via_ssh(client, f"ls {remote_combined_dir}/*.fasta")
    if status != 0:
        print(f"Error listing gene files: {err.strip()}")
        return

    files_list = [line.strip() for line in out.split('\n') if line.strip()]
    if not files_list:
        print("No gene FASTA files found on remote server.")
        return

    def align_file(file_path):
        dirname = os.path.dirname(file_path)
        aligned_dir = os.path.join(dirname, "aligned")
        base_filename = os.path.basename(file_path)
        afa_filename = base_filename.replace('.fasta', '.afa')
        output_path = os.path.join(aligned_dir, afa_filename).replace("\\", "/")

        out_check, err_check, status_check = execute_command_via_ssh(client, f"test -f {output_path} && echo exists")
        if out_check.strip() == 'exists':
            print(f"Skipping {file_path}, already aligned.")
            return

        muscle_command = f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {file_path} -out {output_path}"
        max_attempts = 3
        for attempt in range(1, max_attempts + 1):
            print(f"Aligning {file_path} using MUSCLE... (Attempt {attempt}/{max_attempts})")
            out_exec, err_exec, status_exec = execute_command_via_ssh(client, muscle_command)
            if status_exec == 0:
                print(f"Successfully aligned {file_path} -> {output_path}")
                return
            else:
                if "resource temporarily unavailable" in err_exec.lower():
                    if attempt < max_attempts:
                        print(f"Failed: {err_exec.strip()} - Retrying...")
                    else:
                        print(f"Failed after {max_attempts} attempts: {err_exec.strip()}")
                    continue
                else:
                    print(f"Failed to align {file_path}: {err_exec.strip()}")
                    return

    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(align_file, f) for f in files_list]
        concurrent.futures.wait(futures)

    print("All combined gene alignments completed.")


def download_aligned_combined_files(ssh_client, local_directory, remote_directory):
    sftp = ssh_client.open_sftp()
    ref_protein_dir = os.path.join(local_directory, "CDS_Reference", "All_References", "CDS_Reference_protein")
    aligned_dir_local = os.path.join(ref_protein_dir, "aligned")
    os.makedirs(aligned_dir_local, exist_ok=True)

    remote_aligned_dir = f"{remote_directory}CDS_Reference_protein/aligned"
    try:
        files = sftp.listdir(remote_aligned_dir)
        afa_files = [f for f in files if f.endswith('.afa')]
        if not afa_files:
            print(f"No .afa files to download in {remote_aligned_dir}")
        else:
            for file_name in afa_files:
                remote_file_path = os.path.join(remote_aligned_dir, file_name).replace("\\", "/")
                local_file_path = os.path.join(aligned_dir_local, file_name)
                if os.path.exists(local_file_path):
                    print(f"{file_name} already exists locally. Skipping.")
                    continue
                temp_file_path = os.path.join(aligned_dir_local, "temp_" + file_name)
                print(f"Downloading: {remote_file_path} -> {temp_file_path}")
                sftp.get(remote_file_path, temp_file_path)
                if os.path.getsize(temp_file_path) > 0:
                    shutil.move(temp_file_path, local_file_path)
                    print(f"Downloaded and verified {file_name}")
                else:
                    print(f"Downloaded file {file_name} is empty.")
                    os.remove(temp_file_path)
    except IOError as e:
        print(f"Failed to list files in {remote_aligned_dir}: {str(e)}")

    sftp.close()


#########################################
# Step 3: Add Gaps to Nucleotide Sequences (Globally)
#########################################

def add_gaps_to_all_references(base_dir):
    global_aligned_dir = os.path.join(base_dir, "CDS_Reference", "All_References", "CDS_Reference_protein", "aligned")
    if not os.path.exists(global_aligned_dir):
        print(f"Global aligned directory {global_aligned_dir} does not exist.")
        return

    nucleotide_dir = os.path.join(base_dir, "CDS_Reference", "All_References", "CDS_Reference_nucleotide")
    gapped_dir = os.path.join(base_dir, "CDS_Reference", "All_References", "CDS_Reference_nucleotide_gapped")
    os.makedirs(gapped_dir, exist_ok=True)

    # Load all nucleotide sequences by gene into memory
    all_nucleotide_map = {}
    for gene in get_filtered_gene_order():
        gene_file = os.path.join(nucleotide_dir, f"{gene}.fasta")
        if os.path.exists(gene_file):
            gene_map = {rec.id: str(rec.seq) for rec in SeqIO.parse(gene_file, "fasta")}
            all_nucleotide_map[gene] = gene_map
        else:
            print(f"No nucleotide file for {gene}, skipping gap insertion for this gene.")

    # Apply gaps based on protein alignments
    for gene in get_filtered_gene_order():
        gene_afa = os.path.join(global_aligned_dir, f"{gene}.afa")
        if not os.path.exists(gene_afa):
            print(f"No alignment found for {gene}. Skipping gap insertion.")
            continue

        if gene not in all_nucleotide_map:
            continue

        protein_records = list(SeqIO.parse(gene_afa, 'fasta'))
        nucleotide_map = all_nucleotide_map[gene]
        new_records = []

        for protein_record in protein_records:
            if protein_record.id not in nucleotide_map:
                continue

            gene_seq = nucleotide_map[protein_record.id]
            new_gene_seq = []
            gene_index = 0

            for aa in protein_record.seq:
                if aa == '-':
                    new_gene_seq.append('---')
                else:
                    if gene_index + 3 <= len(gene_seq):
                        codon = gene_seq[gene_index:gene_index + 3]
                        gene_index += 3
                        new_gene_seq.append(codon)
                    else:
                        new_gene_seq.append('---')

            new_gene_seq_str = ''.join(new_gene_seq)
            new_gene_record = SeqRecord(Seq(new_gene_seq_str), id=protein_record.id, description="")
            new_records.append(new_gene_record)

        # Save gapped nucleotide alignment for this gene
        gapped_gene_file = os.path.join(gapped_dir, f"{gene}.fasta")
        SeqIO.write(new_records, gapped_gene_file, 'fasta')
        print(f"Gapped gene sequences for {gene} saved to {gapped_gene_file}")


#########################################
# Step 4: Concatenate Gapped Sequences
#########################################

def process_gene_sequences(gene_sequences):
    if not gene_sequences:
        return gene_sequences

    # Remove trailing gaps
    for acc in gene_sequences:
        gene_sequences[acc] = gene_sequences[acc].rstrip('-')

    # Truncate to shortest length
    min_length = min(len(seq) for seq in gene_sequences.values())
    for acc in gene_sequences:
        gene_sequences[acc] = gene_sequences[acc][:min_length]

    # Make length divisible by 3
    remainder = min_length % 3
    if remainder != 0:
        min_length -= remainder
        for acc in gene_sequences:
            gene_sequences[acc] = gene_sequences[acc][:min_length]

    if min_length == 0:
        return {acc: "" for acc in gene_sequences}

    # Remove codons containing 'N'
    num_codons = min_length // 3
    codon_mask = [True] * num_codons
    for i in range(num_codons):
        start = i * 3
        end = start + 3
        for seq in gene_sequences.values():
            codon = seq[start:end]
            if 'N' in codon:
                codon_mask[i] = False
                break

    for acc in gene_sequences:
        seq = gene_sequences[acc]
        filtered_codons = [seq[i * 3:(i * 3) + 3] for i in range(num_codons) if codon_mask[i]]
        gene_sequences[acc] = "".join(filtered_codons)

    if gene_sequences:
        final_min_length = min(len(s) for s in gene_sequences.values())
        for acc in gene_sequences:
            gene_sequences[acc] = gene_sequences[acc][:final_min_length]

    return gene_sequences


def concatenate_all_gapped_sequences(base_dir):
    gapped_dir = os.path.join(base_dir, "CDS_Reference", "All_References", "CDS_Reference_nucleotide_gapped")
    all_dir = os.path.join(base_dir, "CDS_Reference", "All_MEGA_Grouped_Individually")
    os.makedirs(all_dir, exist_ok=True)

    output_file = os.path.join(all_dir, "All_concatenated_gapped_sequences.fasta")

    # Load and concatenate all genes for each sequence
    concatenated_sequences = {}
    filtered_genes = get_filtered_gene_order()

    # First, gather all sequences by gene
    for gene in filtered_genes:
        gene_file = os.path.join(gapped_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            print(f"Warning: {gene_file} does not exist. Skipping this gene.")
            continue

        gene_sequences = {}
        for record in SeqIO.parse(gene_file, "fasta"):
            acc = record.id
            # Remove appended gene name if any
            for g_name in GENE_ORDER:
                if f"_{g_name}" in acc:
                    acc = acc.replace(f"_{g_name}", "")
                    break
            gene_sequences[acc] = str(record.seq)

        processed_sequences = process_gene_sequences(gene_sequences)

        for acc, seq in processed_sequences.items():
            if acc not in concatenated_sequences:
                concatenated_sequences[acc] = ""
            concatenated_sequences[acc] += seq

    # Write concatenated file
    with open(output_file, "w") as output_fh:
        for acc, sequence in concatenated_sequences.items():
            output_fh.write(f">{acc}\n{sequence}\n")

    return output_file


#########################################
# Create single grp file (modified to group by Genus)
#########################################

def create_single_grp_file(combined_fasta):
    if combined_fasta is None:
        return None

    grp_dir = os.path.dirname(combined_fasta)
    grp_file = os.path.join(grp_dir, "All_sequences.grp")

    if os.path.exists(grp_file):
        print(f"{grp_file} already exists. Skipping .grp creation.")
        return grp_file

    with open(combined_fasta, "r") as f, open(grp_file, "w") as grp:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip()
                # Split to extract genus
                # Example seq_id: NC_021435.1_Ziphius_cavirostris
                # Splitting by underscore: ["NC", "021435.1", "Ziphius", "cavirostris"]
                seq_parts = seq_id.split("_")
                if len(seq_parts) >= 4:
                    genus = seq_parts[2]
                else:
                    # Fallback if unexpected format
                    genus = "UnknownGenus"

                # Write: seq_id=Genus
                grp.write(f"{seq_id}={genus}\n")

    print(f"Created single grp file: {grp_file}")
    return grp_file


#########################################
# Run MEGA once on the combined dataset
#########################################

def run_single_mega_analysis(combined_fasta, grp_file):
    megacc_executable = r"C:\Program Files\MEGA11\megacc.exe"
    if not os.path.isfile(megacc_executable):
        print(f"Error: MEGA executable not found at {megacc_executable}")
        return

    settings_files = [
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_nucleotide.mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(nonsynonymous_only).mao",
        r"C:\Users\freem\OneDrive\Documents\USC\Honours\Analysis\MegaCC settings\distance_estimation_within_grp_avg_syn-nonsynonymous(Synonymous_only).mao"
    ]

    output_dir = os.path.dirname(grp_file)
    for mao_file in settings_files:
        if not os.path.isfile(mao_file):
            print(f"Warning: Settings file not found: {mao_file}")
            continue

        analysis_type = os.path.basename(mao_file).replace(".mao", "")
        if analysis_type == "distance_estimation_within_grp_avg_nucleotide":
            analysis_type = "d_value"
        elif "nonsynonymous_only" in analysis_type:
            analysis_type = "nonsynonymous_only"
        elif "Synonymous_only" in analysis_type:
            analysis_type = "Synonymous_only"

        output_file = os.path.join(output_dir, f"All_{analysis_type}.meg")
        if os.path.exists(output_file):
            print(f"MEGA result {output_file} already exists. Skipping this analysis.")
            continue

        command = [
            megacc_executable,
            "-a", mao_file,
            "-d", combined_fasta,
            "-g", grp_file,
            "-o", output_file
        ]

        print(f"Running MEGA analysis {analysis_type}...")
        try:
            subprocess.run(command, check=True)
            print(f"Analysis '{analysis_type}' completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running analysis '{analysis_type}': {e}")


#########################################
# Compile Results to Excel (Optional)
#########################################

def compile_results_to_excel():
    root_dir = 'sequences/CDS_Reference/All_MEGA_Grouped_Individually'
    if not os.path.exists(root_dir):
        print(f"{root_dir} not found, skipping Excel compilation.")
        return

    pattern = re.compile(r'^(.+?)\((.*?)\)\s+(\S+)\s+(\S+)$')
    results = {}

    def parse_file(filepath, value_key, se_key):
        if os.path.exists(filepath):
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    m = pattern.match(line)
                    if m:
                        species_name, _, value_str, se_str = m.groups()
                        if species_name not in results:
                            results[species_name] = {}
                        results[species_name][value_key] = float(value_str)
                        results[species_name][se_key] = float(se_str)

    d_file = os.path.join(root_dir, "All_d_value.meg")
    nonsyn_file = os.path.join(root_dir, "All_nonsynonymous_only.meg")
    syn_file = os.path.join(root_dir, "All_Synonymous_only.meg")

    parse_file(d_file, 'd', 'd_SE')
    parse_file(nonsyn_file, 'dN', 'dN_SE')
    parse_file(syn_file, 'dS', 'dS_SE')

    combined_fasta = os.path.join(root_dir, "All_concatenated_gapped_sequences.fasta")
    if not os.path.exists(combined_fasta):
        print("No combined FASTA found, skipping GC computation.")
        return

    species_sequences = {}
    for record in SeqIO.parse(combined_fasta, "fasta"):
        parts = record.id.split('_', 1)
        species_name = parts[1] if len(parts) == 2 else parts[0]
        species_sequences.setdefault(species_name, []).append(str(record.seq))

    avg_gc_per_species = {}
    sequence_count_per_species = {}
    for sp_name, seq_list in species_sequences.items():
        full_seq = ''.join(seq_list)
        avg_gc_per_species[sp_name] = gc_fraction(full_seq) * 100
        sequence_count_per_species[sp_name] = len(seq_list)

    rows = []
    for species, metrics in results.items():
        d = metrics.get('d')
        d_SE = metrics.get('d_SE')
        dN = metrics.get('dN')
        dN_SE = metrics.get('dN_SE')
        dS = metrics.get('dS')
        dS_SE = metrics.get('dS_SE')
        gc_value = avg_gc_per_species.get(species)
        seq_count = sequence_count_per_species.get(species)

        row = {
            'Species': species,
            'd': d,
            'd S.E.': d_SE,
            'dN': dN,
            'dN S.E.': dN_SE,
            'dS': dS,
            'dS S.E.': dS_SE,
            'Avg genome GC (%)': gc_value,
            'Sequence No.': seq_count
        }
        rows.append(row)

    df = pd.DataFrame(rows, columns=[
        'Species', 'd', 'd S.E.', 'dN', 'dN S.E.',
        'dS', 'dS S.E.', 'Avg genome GC (%)', 'Sequence No.'
    ])

    output_path = os.path.join(root_dir, 'All_Grouped_by_Species.xlsx')
    if os.path.exists(output_path):
        print(f"{output_path} already exists. Overwriting.")
    df.to_excel(output_path, index=False)
    print(f"Data successfully written to {output_path}")


#########################################
# Align All Combined Sequences Together Locally & Remove Gaps/Ns (Codon-based)
#########################################

def remove_gaps_and_ns_from_alignment(aligned_file):
    """
    Remove any codon positions from the alignment that contain '-' or 'N'.
    We consider the alignment length should be divisible by 3 (if not, we truncate).
    We remove entire codon triplets if any sequence has 'N' or '-' in that codon.
    """
    if aligned_file is None or not os.path.exists(aligned_file):
        print("No alignment file found to clean.")
        return None

    records = list(SeqIO.parse(aligned_file, "fasta"))
    if not records:
        print("No records found in alignment.")
        return None

    seqs = [str(r.seq) for r in records]
    length = len(seqs[0])
    # Ensure all sequences are same length
    if any(len(s) != length for s in seqs):
        print("Sequences are not the same length, cannot clean reliably.")
        return None

    # Make length divisible by 3 by truncation if needed
    remainder = length % 3
    if remainder != 0:
        length -= remainder
        seqs = [s[:length] for s in seqs]

    num_codons = length // 3
    codon_indices_to_keep = []
    for i in range(num_codons):
        start = i * 3
        end = start + 3
        codons = [s[start:end] for s in seqs]
        # If any codon contains '-' or 'N', we discard this codon position
        if any('N' in c or '-' in c for c in codons):
            continue
        codon_indices_to_keep.append(i)

    if not codon_indices_to_keep:
        print("All codon columns removed. No clean alignment possible.")
        return None

    # Rebuild sequences
    cleaned_seqs = []
    for s in seqs:
        new_seq = []
        for i in codon_indices_to_keep:
            start = i * 3
            end = start + 3
            new_seq.append(s[start:end])
        cleaned_seqs.append("".join(new_seq))

    cleaned_records = []
    for rec, seq_str in zip(records, cleaned_seqs):
        cleaned_records.append(SeqRecord(Seq(seq_str), id=rec.id, description=""))

    cleaned_file = aligned_file.replace(".afa", "_cleaned.fasta")
    SeqIO.write(cleaned_records, cleaned_file, "fasta")
    print(f"Cleaned codon-based alignment written to {cleaned_file}")

    return cleaned_file


#########################################
# Main Execution
#########################################

if __name__ == "__main__":
    base_directory = "sequences/"
    process_genbank_files(base_directory)
    create_combined_protein_files(base_directory)

    # SSH server configuration
    remote_server = "203.101.229.234"
    remote_port = 22
    remote_user = "mfreeman"
    remote_key_file = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"
    remote_dir = "/home/mfreeman/USCServer/CDS_Reference/"

    ssh_client = create_ssh_client(remote_server, remote_port, remote_user, remote_key_file)
    ensure_remote_combined_dir(ssh_client, remote_dir)
    upload_combined_files(ssh_client, base_directory, remote_dir)
    run_muscle_alignment_on_combined(ssh_client, remote_dir)
    download_aligned_combined_files(ssh_client, base_directory, remote_dir)
    ssh_client.close()

    add_gaps_to_all_references(base_directory)
    combined_fasta = concatenate_all_gapped_sequences(base_directory)
    # Final alignment and cleaning (codon-based gap and N removal)
    cleaned_alignment = remove_gaps_and_ns_from_alignment(combined_fasta)
    grp_file = create_single_grp_file(combined_fasta)
    run_single_mega_analysis(combined_fasta, grp_file)
    compile_results_to_excel()

    print("All steps completed. Final cleaned alignment at:")
    print(cleaned_alignment if cleaned_alignment else "No cleaned alignment produced.")

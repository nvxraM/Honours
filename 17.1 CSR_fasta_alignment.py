#!/usr/bin/env python3
"""
align_csr_fastas.py

Script to upload CSR FASTA files to a remote server, run MUSCLE alignments via Singularity,
and download the resulting .afa alignment files back to the local machine.
"""

import os
import paramiko
import concurrent.futures
from pathlib import Path

# ---------------------- Configuration ----------------------
SERVER         = "203.101.229.234"
PORT           = 22
USER           = "mfreeman"
KEY_FILE       = "C:/Users/freem/Downloads/mfreeman-private-key.txt"

LOCAL_BASE     = Path("sequences/CDS_CSR")
REMOTE_BASE    = "/home/mfreeman/USCServer/CDS_CSR"
SINGULARITY_IMG= "/RDS/Q1233/singularity/muscle.sif"
MAX_WORKERS    = 5

# ---------------------- Helper Functions ----------------------

def create_ssh_client(server, port, user, key_file):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port=port, username=user, key_filename=key_file)
    return client

def exec_cmd(ssh_client, cmd):
    stdin, stdout, stderr = ssh_client.exec_command(cmd)
    code = stdout.channel.recv_exit_status()
    return stdout.read().decode(), stderr.read().decode(), code

# ---------------------- Remote Directory Setup ----------------------

def ensure_remote_dirs(ssh_client):
    """Create remote base + aligned subdirs for each genus."""
    cmds = []
    for genus_dir in LOCAL_BASE.iterdir():
        if not genus_dir.is_dir(): continue
        genus = genus_dir.name
        # ensure both the genus folder and its 'aligned' subfolder exist
        cmds.append(f"mkdir -p {REMOTE_BASE}/{genus}/aligned")
    if cmds:
        out, err, code = exec_cmd(ssh_client, " && ".join(cmds))
        if code != 0:
            print(f"[ERROR] ensure_remote_dirs: {err.strip()}")
        else:
            print("[INFO] Remote directories ensured.")

# ---------------------- Upload CSR FASTAs (idempotent) ----------------------

def upload_csr(ssh_client):
    sftp = ssh_client.open_sftp()
    try:
        for genus_dir in LOCAL_BASE.iterdir():
            if not genus_dir.is_dir(): continue
            genus = genus_dir.name
            remote_genus = f"{REMOTE_BASE}/{genus}"
            for fasta in genus_dir.glob("*.fasta"):
                local_path  = str(fasta)
                remote_path = f"{remote_genus}/{fasta.name}"
                # only upload if it doesn't already exist remotely
                try:
                    sftp.stat(remote_path)
                except IOError:
                    print(f"[UPLOAD] {local_path} → {remote_path}")
                    sftp.put(local_path, remote_path)
                else:
                    print(f"[SKIP] remote already has {fasta.name}")
    finally:
        sftp.close()

# ---------------------- MUSCLE Alignment on Server ----------------------

def run_muscle(ssh_client, fasta_path):
    # ensure the aligned/ folder exists
    aligned_dir = os.path.dirname(fasta_path) + "/aligned"
    exec_cmd(ssh_client, f"mkdir -p {aligned_dir}")

    # determine remote output file
    out_file = os.path.join(aligned_dir,
                            os.path.basename(fasta_path).replace('.fasta', '.afa'))

    # skip if it’s already been aligned
    out, err, code = exec_cmd(ssh_client, f"test -f {out_file} && echo EXISTS")
    if "EXISTS" in out:
        print(f"[SKIP] {os.path.basename(out_file)} already aligned.")
        return

    # run MUSCLE via singularity
    cmd = f"singularity exec {SINGULARITY_IMG} muscle -in {fasta_path} -out {out_file}"
    print(f"[RUN] {cmd}")
    out, err, code = exec_cmd(ssh_client, cmd)
    if code == 0:
        print(f"[OK] Aligned {os.path.basename(fasta_path)}")
    else:
        print(f"[FAIL] {fasta_path}: {err.strip()}")

def align_all(ssh_client):
    # find every .fasta under REMOTE_BASE (excluding anything already under aligned/)
    find_cmd = f"find {REMOTE_BASE} -type f -name '*.fasta' ! -path '*/aligned/*'"
    out, err, code = exec_cmd(ssh_client, find_cmd)
    if code != 0:
        print(f"[ERROR] finding fastas: {err.strip()}")
        return

    files = [line for line in out.split() if line]
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        for fasta in files:
            pool.submit(run_muscle, ssh_client, fasta)

# ---------------------- Download Alignments (idempotent) ----------------------

def download_alignments(ssh_client):
    sftp = ssh_client.open_sftp()
    try:
        for genus_dir in LOCAL_BASE.iterdir():
            if not genus_dir.is_dir():
                continue

            # make sure local aligned/ exists
            local_aligned = genus_dir / "aligned"
            local_aligned.mkdir(exist_ok=True)

            # two places to check on the server
            remote_base_genus = f"{REMOTE_BASE}/{genus_dir.name}"
            remote_aligned     = remote_base_genus + "/aligned"
            remote_roots       = [remote_aligned, remote_base_genus]

            for remote_dir in remote_roots:
                try:
                    for fname in sftp.listdir(remote_dir):
                        if not fname.endswith(".afa"):
                            continue
                        remote_file = f"{remote_dir}/{fname}"
                        local_file  = local_aligned / fname

                        if local_file.exists():
                            print(f"[SKIP] {genus_dir.name}/{fname} already downloaded.")
                        else:
                            print(f"[DL] {remote_file} → {local_file}")
                            sftp.get(remote_file, str(local_file))
                except IOError:
                    # only warn about missing aligned/ subfolder
                    if remote_dir == remote_aligned:
                        print(f"[WARN] No remote aligned folder for {genus_dir.name}")
                    # if remote_base_genus is missing (unlikely), just skip

    finally:
        sftp.close()

# ---------------------- Main Workflow ----------------------

if __name__ == "__main__":
    client = create_ssh_client(SERVER, PORT, USER, KEY_FILE)
    ensure_remote_dirs(client)
    upload_csr(client)
    align_all(client)
    download_alignments(client)
    client.close()

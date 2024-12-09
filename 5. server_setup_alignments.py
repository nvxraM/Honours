import paramiko
import subprocess
import os
import shutil
import concurrent.futures

# ------------------------------------------------------------
# Helper Functions
# ------------------------------------------------------------

def run_command(command):
    """
    Run a local shell command and capture its output, error, and exit code.
    Returns a tuple (output, error, exit_code).
    """
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode


def create_ssh_client(server, port, user, key_file):
    """
    Create and return an SSH client connected to the specified server.
    Uses a given user and private key file for authentication.
    """
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, username=user, key_filename=key_file)
    return client


def execute_command_via_ssh(client, command):
    """
    Execute a command on the remote server via the provided SSH client.
    Returns a tuple (stdout, stderr, exit_code).
    """
    stdin, stdout, stderr = client.exec_command(command)
    exit_code = stdout.channel.recv_exit_status()
    output = stdout.read().decode()
    error = stderr.read().decode()
    return output, error, exit_code


# ------------------------------------------------------------
# Directory Management
# ------------------------------------------------------------

def ensure_directories_on_server(ssh_client, local_base_directory, remote_base_directory):
    """
    Ensure that for each species directory found locally, the necessary remote directories exist:
    specifically the `CDS_protein/aligned` directory under each species.

    Steps:
    - Identify species directories from the local base directory.
    - Construct and run a single combined command on the server to create all these directories.
    """
    # Identify species directories locally (each species is a subdirectory)
    species_dirs = [
        name for name in os.listdir(local_base_directory)
        if os.path.isdir(os.path.join(local_base_directory, name))
    ]

    # Prepare a single combined command to create all aligned directories in one go
    cmd_list = []
    for specie in species_dirs:
        cmd_list.append(f"mkdir -p {remote_base_directory}{specie}/CDS_protein/aligned")

    combined_cmd = " && ".join(cmd_list)
    if combined_cmd.strip():
        _, error, exit_code = execute_command_via_ssh(ssh_client, combined_cmd)
        if exit_code != 0:
            print(f"Error ensuring directories: {error}")
        else:
            print("All required directories on server have been ensured.")
    else:
        print("No species directories found locally. No directories created on server.")


# ------------------------------------------------------------
# MUSCLE Alignment Execution with Retry
# ------------------------------------------------------------

def run_muscle_alignment(client, file):
    """
    Run MUSCLE alignment on a single .fasta file if not already aligned.
    The output is stored in 'CDS_protein/aligned' with a .afa extension.

    If the alignment fails due to a "resource temporarily unavailable" error, it will retry
    up to 3 times before giving up.
    """
    # Determine output alignment file path
    dirname = os.path.dirname(file)  # e.g., .../CDS_protein
    aligned_dir = os.path.join(dirname, "aligned")
    base_filename = os.path.basename(file)  # e.g., ND6.fasta
    afa_filename = base_filename.replace('.fasta', '.afa')
    output_path = os.path.join(aligned_dir, afa_filename).replace("\\", "/")

    # Check if alignment already exists on the server
    out, err, status = execute_command_via_ssh(client, f"test -f {output_path} && echo exists")
    if out.strip() == 'exists':
        print(f"Skipping {file}, already aligned.")
        return

    # MUSCLE command via Singularity container
    muscle_command = (
        f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {file} -out {output_path}"
    )

    max_attempts = 3
    for attempt in range(1, max_attempts + 1):
        print(f"Aligning {file} using MUSCLE... (Attempt {attempt}/{max_attempts})")
        out, err, status = execute_command_via_ssh(client, muscle_command)

        if status == 0:
            # Exit code 0 indicates success
            print(f"Successfully aligned {file} -> {output_path}")
            return
        else:
            # Check if it's the known "resource temporarily unavailable" error
            if "resource temporarily unavailable" in err.lower():
                # If we haven't reached max attempts, retry
                if attempt < max_attempts:
                    print(f"Failed to align {file}: {err.strip()} - Retrying...")
                    continue
                else:
                    # After 3 attempts, give up
                    print(f"Failed to align {file} after {max_attempts} attempts: {err.strip()}")
                    print("Tried 3 times giving up on this one...")
                    return
            else:
                # Some other error occurred, no point in retrying
                print(f"Failed to align {file}: {err.strip()}")
                return


def execute_muscle_on_cds(client, base_directory):
    """
    Executes MUSCLE alignment on all *.fasta files found in CDS_protein directories on the server.

    Steps:
    - Use 'find' to locate all .fasta files in subdirectories of 'CDS_protein'.
    - Run alignments in parallel with ThreadPoolExecutor(max_workers=5).
    """
    find_command = f"find {base_directory} -type f -path '*/CDS_protein/*.fasta'"
    output, error, exit_code = execute_command_via_ssh(client, find_command)

    if exit_code != 0:
        print(f"Error retrieving .fasta files from CDS_protein directories: {error.strip()}")
        return

    # Create a list of all found .fasta files
    files_list = [file.strip() for file in output.strip().split('\n') if file.strip()]

    # Use parallelization to speed up the alignment process
    # max_workers=5 requested
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = [executor.submit(run_muscle_alignment, client, file) for file in files_list]
        concurrent.futures.wait(futures)

    print("Alignment process completed for all .fasta files in CDS_protein directories.")


# ------------------------------------------------------------
# File Transfer (Upload and Download)
# ------------------------------------------------------------

def upload_files(ssh_client, local_directory, remote_directory):
    """
    Upload all .fasta files found in local CDS_protein directories to the corresponding
    remote directories. The directory structure under each species is preserved.

    Steps:
    - Recursively search local_directory for files ending in .fasta within CDS_protein dirs.
    - For each .fasta file, replicate the directory structure on the remote server.
    - Upload the file.
    """
    sftp = ssh_client.open_sftp()
    try:
        for root, _, files in os.walk(local_directory):
            # Only upload files from CDS_protein directories
            if "CDS_protein" not in root:
                continue
            for file in files:
                if file.endswith('.fasta'):
                    local_file_path = os.path.join(root, file)
                    # Compute relative path to mirror directory structure on server
                    relative_path = os.path.relpath(local_file_path, local_directory)
                    remote_file_path = os.path.join(remote_directory, relative_path).replace("\\", "/")

                    # Ensure that the remote directory structure exists
                    remote_dir = os.path.dirname(remote_file_path)
                    _, dir_err, dir_status = execute_command_via_ssh(ssh_client, f'mkdir -p {remote_dir}')
                    if dir_status != 0:
                        print(f"Error ensuring remote directory: {dir_err.strip()}")

                    print(f"Uploading {local_file_path} to {remote_file_path}")
                    sftp.put(local_file_path, remote_file_path)

    except Exception as e:
        print(f"Error while uploading files: {str(e)}")
    finally:
        sftp.close()


def download_files(ssh_client, base_directory, local_directory):
    """
    Download aligned .afa files from the remote server back to the local environment.
    These are placed within the corresponding CDS_protein/aligned directories locally.

    Steps:
    - For each species directory locally, check the corresponding remote aligned directory.
    - Download all .afa files if they don't already exist locally.
    - Verify file integrity after download.
    """
    sftp = ssh_client.open_sftp()
    try:
        for species_dir in os.listdir(local_directory):
            local_aligned_dir = os.path.join(local_directory, species_dir, "CDS_protein", "aligned")
            remote_aligned_dir = f"{base_directory}{species_dir}/CDS_protein/aligned"
            os.makedirs(local_aligned_dir, exist_ok=True)

            # Attempt to list .afa files remotely
            try:
                files = sftp.listdir(remote_aligned_dir)
                afa_files = [file for file in files if file.endswith('.afa')]
                if not afa_files:
                    print(f"No .afa files to download in {remote_aligned_dir}")
                    continue
                print(f"Found .afa files in {remote_aligned_dir}: {afa_files}")
            except IOError as e:
                print(f"Failed to list files in {remote_aligned_dir}: {str(e)}")
                continue

            # Download each .afa file if it doesn't exist locally
            for file in afa_files:
                remote_file_path = os.path.join(remote_aligned_dir, file).replace("\\", "/")
                local_file_path = os.path.join(local_aligned_dir, file)
                temp_file_path = os.path.join(local_aligned_dir, "temp_" + file)

                if os.path.exists(local_file_path):
                    # If the file already exists locally, skip to avoid overwriting
                    print(f"{file} already exists locally and won't be downloaded.")
                else:
                    print(f"Downloading: {remote_file_path} to {temp_file_path}")
                    try:
                        sftp.get(remote_file_path, temp_file_path)
                        # Verify downloaded file size is non-zero
                        if os.path.getsize(temp_file_path) > 0:
                            shutil.move(temp_file_path, local_file_path)
                            print(f"Downloaded and verified {file} to {local_file_path}")
                        else:
                            print(f"Downloaded file {file} is empty. Check server content.")
                            os.remove(temp_file_path)
                    except Exception as e:
                        print(f"Failed to download {file}: {str(e)}")
                        if os.path.exists(temp_file_path):
                            os.remove(temp_file_path)
    finally:
        sftp.close()


# ------------------------------------------------------------
# Main Workflow
# ------------------------------------------------------------

if __name__ == "__main__":
    # Configuration for SSH and directories
    server = "203.101.229.234"
    port = 22
    user = "mfreeman"
    key_file = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"

    remote_dir = "/home/mfreeman/USCServer/CDS/"
    local_base_directory = "sequences/CDS"  # Local base directory containing species subdirs
    local_dir = "sequences/CDS"             # Same as above, can be used interchangeably if needed

    # 1. Create an SSH client connection
    client = create_ssh_client(server, port, user, key_file)

    # 2. Ensure the required remote directories exist
    ensure_directories_on_server(client, local_base_directory, remote_dir)

    # 3. Upload local .fasta files to the server (mirroring directory structure)
    upload_files(client, local_dir, remote_dir)

    # 4. Execute MUSCLE alignment on all .fasta files in CDS_protein directories on the server
    execute_muscle_on_cds(client, remote_dir)

    # 5. Download the aligned .afa files back to local directories
    download_files(client, remote_dir, local_dir)

    # 6. Close the SSH connection
    client.close()

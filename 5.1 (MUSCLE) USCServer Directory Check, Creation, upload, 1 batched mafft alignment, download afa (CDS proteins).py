import paramiko
import subprocess
import os
import shutil
import time
import concurrent.futures

def run_command(command):
    """
    Run a shell command locally and return the output, error, and exit code.

    Parameters
    ----------
    command : str
        The shell command to execute.

    Returns
    -------
    tuple
        (output, error, return_code) from the executed command.
    """
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode

def create_ssh_client(server, port, user, key_file):
    """
    Create an SSH client and connect to the remote server using the provided credentials.

    Parameters
    ----------
    server : str
        Hostname or IP address of the remote server.
    port : int
        SSH port number (usually 22).
    user : str
        Username for SSH login.
    key_file : str
        Path to the private key file for authentication.

    Returns
    -------
    paramiko.SSHClient
        An active SSH client session connected to the server.
    """
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, username=user, key_filename=key_file)
    return client

def ensure_directories_on_server(server, port, user, key_file):
    """
    Ensure that the required directory structure exists on the remote server.
    This function checks if directories exist and creates them if necessary.

    The structure is:
    /home/mfreeman/USCServer/sequences/CDS/{species}/CDS_nucleotide
    /home/mfreeman/USCServer/sequences/CDS/{species}/CDS_protein

    Parameters
    ----------
    server : str
        Remote server hostname or IP.
    port : int
        SSH port.
    user : str
        Username for SSH.
    key_file : str
        Path to the SSH private key file.
    """
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/"
    local_base_directory = "sequences/CDS"

    # Connect via SSH
    ssh = create_ssh_client(server, port, user, key_file)

    # Get all species directories from the local sequences/CDS directory
    species_dirs = [
        name for name in os.listdir(local_base_directory)
        if os.path.isdir(os.path.join(local_base_directory, name))
    ]

    # Ensure directories exist on the server for each species
    for specie in species_dirs:
        commands = [
            f'test -d {base_directory}{specie} || mkdir -p {base_directory}{specie}',
            f'test -d {base_directory}{specie}/CDS_nucleotide || mkdir -p {base_directory}{specie}/CDS_nucleotide',
            f'test -d {base_directory}{specie}/CDS_protein || mkdir -p {base_directory}{specie}/CDS_protein'
        ]
        for command in commands:
            output, error = execute_command_via_ssh(ssh, command)
            if error:
                print(f"Error creating directory {specie} on server: {error}")

    ssh.close()
    print("All directories created successfully.")

def execute_command_via_ssh(client, command):
    """
    Execute a command on the remote server via SSH and return its output and error.

    Parameters
    ----------
    client : paramiko.SSHClient
        The connected SSH client instance.
    command : str
        The command to run on the remote server.

    Returns
    -------
    tuple
        (output_str, error_str) from the remote command execution.
    """
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read().decode()
    error = stderr.read().decode()
    return output, error

def execute_muscle_on_cds(server, port, user, key_file):
    """
    Execute MUSCLE alignments on all .fasta files found in the remote server's
    CDS directories. Output alignments are saved as .afa files in the CDS_protein/aligned directory.

    Parameters
    ----------
    server : str
        Remote server hostname or IP.
    port : int
        SSH port.
    user : str
        Username for SSH.
    key_file : str
        Path to the SSH private key file.
    """
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/"
    client = create_ssh_client(server, port, user, key_file)

    try:
        # Find all .fasta files on the server under the base directory
        find_cmd = f"find {base_directory} -type f -name '*.fasta'"
        output, error = execute_command_via_ssh(client, find_cmd)
        if error:
            print(f"Error retrieving file list from {base_directory}: {error}")
            return

        files_list = [file for file in output.strip().split() if file.endswith('.fasta')]

        def run_muscle_alignment(file):
            """
            Run MUSCLE alignment on a single .fasta file.

            Parameters
            ----------
            file : str
                The remote path to the .fasta file.
            """
            # Output path for the alignment file (change directory from CDS_nucleotide to CDS_protein/aligned)
            output_path = file.replace(".fasta", ".afa").replace("CDS_nucleotide", "CDS_protein/aligned")

            # Check if the alignment output file already exists to avoid re-running
            check_cmd = f"test -f {output_path} && echo exists"
            out, err = execute_command_via_ssh(client, check_cmd)
            if out.strip() == 'exists':
                print(f"Skipping {file}, alignment already exists.")
                return

            # Run MUSCLE alignment using the singularity container
            align_cmd = (
                f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {file} -out {output_path}"
            )
            print(f"Aligning {file} using MUSCLE...")
            out, err = execute_command_via_ssh(client, align_cmd)
            if err:
                print(f"Failed to align {file}: {err}")

        # Use a ThreadPoolExecutor to run alignments concurrently if desired
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            futures = [executor.submit(run_muscle_alignment, file) for file in files_list]
            concurrent.futures.wait(futures)

    finally:
        client.close()

    print("Alignment process completed for all specified files.")

def upload_files(ssh_client, local_directory, remote_directory):
    """
    Upload all .fasta files from a local directory to a corresponding directory structure on the remote server.
    This function ensures the remote directory structure is maintained and mirrored.

    Parameters
    ----------
    ssh_client : paramiko.SSHClient
        Connected SSH client.
    local_directory : str
        Local base directory containing the sequences.
    remote_directory : str
        Remote base directory to upload the files to.
    """
    sftp = ssh_client.open_sftp()
    try:
        for root, _, files in os.walk(local_directory):
            for file in files:
                if file.endswith('.fasta'):
                    local_file_path = os.path.join(root, file)
                    relative_path = os.path.relpath(local_file_path, local_directory)
                    remote_file_path = os.path.join(remote_directory, relative_path).replace("\\", "/")
                    remote_dir = os.path.dirname(remote_file_path)

                    # Ensure remote directory exists
                    ssh_client.exec_command(f'mkdir -p {remote_dir}')

                    print(f"Uploading {local_file_path} to {remote_file_path}")
                    sftp.put(local_file_path, remote_file_path)

    except Exception as e:
        print(f"Error while uploading files: {e}")
    finally:
        sftp.close()

def download_files(ssh_client, base_directory, local_directory):
    """
    Download aligned .afa files from the remote server back to the local machine.
    The local directory structure is assumed to mirror the remote structure.

    Parameters
    ----------
    ssh_client : paramiko.SSHClient
        Connected SSH client.
    base_directory : str
        Remote base directory containing aligned files.
    local_directory : str
        Local base directory to store downloaded aligned files.
    """
    sftp = ssh_client.open_sftp()
    try:
        # Iterate over species directories locally to determine which aligned files to download
        for species_dir in os.listdir(local_directory):
            local_species_aligned_dir = os.path.join(local_directory, species_dir, "CDS_protein", "aligned")
            remote_aligned_dir = f"{base_directory}{species_dir}/CDS_protein/aligned"
            os.makedirs(local_species_aligned_dir, exist_ok=True)

            # List remote .afa files for the species
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

            # Download each .afa file, verifying integrity and avoiding duplicates
            for file in afa_files:
                remote_file_path = os.path.join(remote_aligned_dir, file).replace("\\", "/")
                local_file_path = os.path.join(local_species_aligned_dir, file)
                temp_file_path = os.path.join(local_species_aligned_dir, "temp_" + file)

                if os.path.exists(local_file_path):
                    print(f"{file} already exists locally. Skipping download.")
                    continue

                print(f"Downloading {remote_file_path} to {temp_file_path}")
                try:
                    sftp.get(remote_file_path, temp_file_path)
                    # Verify file is not empty
                    if os.path.getsize(temp_file_path) > 0:
                        shutil.move(temp_file_path, local_file_path)
                        print(f"Downloaded and verified {file} to {local_file_path}")
                    else:
                        print(f"Warning: Downloaded file {file} is empty. Removing it.")
                        os.remove(temp_file_path)
                except Exception as e:
                    print(f"Failed to download {file}: {str(e)}")
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)

    finally:
        sftp.close()


if __name__ == "__main__":
    # Configuration
    server = "203.101.229.234"
    port = 22
    user = "mfreeman"
    key_file = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"
    remote_dir = "/home/mfreeman/USCServer/sequences/CDS/"
    local_dir = "sequences/CDS"

    # Step 1: Ensure required directories exist on the server
    ensure_directories_on_server(server, port, user, key_file)

    # Step 2: Create SSH client
    client = create_ssh_client(server, port, user, key_file)

    # Step 3: Upload CDS files from local machine to the server
    upload_files(client, local_dir, remote_dir)

    # Step 4: Execute MUSCLE alignment on the server
    execute_muscle_on_cds(server, port, user, key_file)

    # Step 5: Download aligned .afa files back to local machine
    client = create_ssh_client(server, port, user, key_file)
    download_files(client, remote_dir, local_dir)

    # Step 6: Close the SSH client
    client.close()

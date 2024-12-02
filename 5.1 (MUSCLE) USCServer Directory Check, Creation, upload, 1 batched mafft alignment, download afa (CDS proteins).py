import paramiko
import subprocess
import os
import shutil
import time
import concurrent.futures

# Function to run a shell command and return the output, error, and exit code
def run_command(command):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode

# Function to create an SSH client and connect to the server
def create_ssh_client(server, port, user, key_file):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, username=user, key_filename=key_file)
    return client

# Function to ensure necessary directories exist on the server
def ensure_directories_on_server(server, port, user, key_file):
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/"

    # Create SSH client to connect to the server
    ssh = create_ssh_client(server, port, user, key_file)

    # Get the list of species directories from the local directory
    local_base_directory = "sequences/CDS"
    species_dirs = [name for name in os.listdir(local_base_directory) if os.path.isdir(os.path.join(local_base_directory, name))]

    # Check if the directories already exist on the server, and create if they don't
    for specie in species_dirs:
        commands = [
            f'test -d {base_directory}{specie} || mkdir -p {base_directory}{specie}',
            f'test -d {base_directory}{specie}/CDS_nucleotide || mkdir -p {base_directory}{specie}/CDS_nucleotide',
            f'test -d {base_directory}{specie}/CDS_protein || mkdir -p {base_directory}{specie}/CDS_protein'
        ]
        for command in commands:
            output, error = execute_command_via_ssh(ssh, command)
            if error:
                print(f"Error creating directory: {error}")

    ssh.close()
    print("All directories created successfully.")

# Function to execute a command via SSH and return output and error
def execute_command_via_ssh(client, command):
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read()
    error = stderr.read()
    return output.decode(), error.decode()

# Function to execute MUSCLE alignment on CDS files on the server
def execute_muscle_on_cds(server, port, user, key_file):
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/"

    # Create SSH client to connect to the server
    client = create_ssh_client(server, port, user, key_file)

    try:
        # List the files in the base directory and subdirectories
        output, error = execute_command_via_ssh(client, f"find {base_directory} -type f -name '*.fasta'")
        if error:
            print(f"Error retrieving list from {base_directory}: {error}")
            return
        files_list = [file for file in output.strip().split()]

        # Function to run MUSCLE alignment for a specific file with enhanced debugging
        def run_muscle_alignment(file):
            output_path = file.replace(".fasta", ".afa").replace("CDS_nucleotide", "CDS_protein/aligned")
            # Check if the alignment already exists
            output, error = execute_command_via_ssh(client, f"test -f {output_path} && echo exists")
            if output.strip() == 'exists':
                print(f"Skipping {file}, already aligned.")
                return

            # Run MUSCLE alignment using Singularity container
            command = f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {file} -out {output_path}"
            print(f"Aligning {file} using MUSCLE...")
            output, error = execute_command_via_ssh(client, command)
            if error:
                print(f"Failed to align {file}: {error}")

        # Run MUSCLE alignment for each file
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            futures = [executor.submit(run_muscle_alignment, file) for file in files_list]
            concurrent.futures.wait(futures)

    finally:
        # Close the SSH client after completing all alignments
        client.close()

    print("Alignment process completed for all specified files.")

# Function to upload files from local to server using SFTP
def upload_files(ssh_client, local_directory, remote_directory):
    sftp = ssh_client.open_sftp()
    try:
        # Upload each .fasta file
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

# Function to download aligned files from server to local
def download_files(ssh_client, base_directory, local_directory):
    sftp = ssh_client.open_sftp()
    try:
        for species_dir in os.listdir(local_directory):
            local_species_dir = os.path.join(local_directory, species_dir, "CDS_protein", "aligned")
            remote_aligned_dir = f"{base_directory}{species_dir}/CDS_protein/aligned"
            os.makedirs(local_species_dir, exist_ok=True)

            try:
                files = sftp.listdir(remote_aligned_dir)
                afa_files = [file for file in files if file.endswith('.afa')]
                if not afa_files:
                    print(f"No .afa files to download in {remote_aligned_dir}")
                    continue
                print(f"Attempting to download .afa files from {remote_aligned_dir}: {afa_files}")
            except IOError as e:
                print(f"Failed to list files in {remote_aligned_dir}: {str(e)}")
                continue

            for file in afa_files:
                remote_file_path = os.path.join(remote_aligned_dir, file).replace("\\", "/")
                local_file_path = os.path.join(local_species_dir, file)
                temp_file_path = os.path.join(local_species_dir, "temp_" + file)
                print(f"Preparing to download: {remote_file_path} to {temp_file_path}")

                if os.path.exists(local_file_path):
                    print(f"{file} already exists and won't be downloaded.")
                else:
                    try:
                        sftp.get(remote_file_path, temp_file_path)
                        if os.path.getsize(temp_file_path) > 0:
                            shutil.move(temp_file_path, local_file_path)
                            print(f"Downloaded and verified {file} to {local_file_path}")
                        else:
                            print(f"Downloaded file {file} is empty, check server content.")
                            os.remove(temp_file_path)
                    except Exception as e:
                        print(f"Failed to download {file}: {str(e)}")
                        if os.path.exists(temp_file_path):
                            os.remove(temp_file_path)
    finally:
        sftp.close()

# Main logic to run the entire workflow
server = "203.101.229.234"
port = 22
user = "mfreeman"
key_file = "C:/Users/freem/OneDrive/Documents/USC/Honours/API keys/mfreeman-private-key.txt"
remote_dir = "/home/mfreeman/USCServer/sequences/CDS/"
local_dir = "sequences/CDS"

# Ensure required directories exist on the server
ensure_directories_on_server(server, port, user, key_file)

# Create SSH client
client = create_ssh_client(server, port, user, key_file)

# Upload CDS files from local to server
upload_files(client, local_dir, remote_dir)

# Execute MUSCLE alignment on CDS files on the server
execute_muscle_on_cds(server, port, user, key_file)

# Download CDS alignment files
client = create_ssh_client(server, port, user, key_file)
download_files(client, remote_dir, local_dir)

# Close SSH client
client.close()

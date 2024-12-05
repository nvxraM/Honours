import paramiko
import subprocess
import os
import shutil
import time
import concurrent.futures

# Function to run a shell command and return the output, error, and exit code
def run_command(command):
    # Using subprocess to execute the command
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode

# Function to create an SSH client and connect to the server
def create_ssh_client(server, port, user, key_file):
    client = paramiko.SSHClient()
    # Automatically add the server's host key if missing
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # Connect to the server with the given credentials
    client.connect(server, port, username=user, key_filename=key_file)
    return client

# Function to ensure necessary directories exist on the server
def ensure_directories_on_server(server, port, user, key_file):
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/"

    # Create SSH client to connect to the server
    ssh = create_ssh_client(server, port, user, key_file)

    # Define the species and the required directory paths
    specie = "Phocoena"
    commands = [
        f'test -d {base_directory}{specie} || mkdir -p {base_directory}{specie}',
        f'test -d {base_directory}{specie}/CDS_protein || mkdir -p {base_directory}{specie}/CDS_protein',
        f'test -d {base_directory}{specie}/CDS_protein/aligned || mkdir -p {base_directory}{specie}/CDS_protein/aligned'
    ]
    # Ensure each directory exists or create it if it doesn't
    for command in commands:
        output, error = execute_command_via_ssh(ssh, command)
        if error:
            print(f"Error creating directory: {error}")

    # Close the SSH client after creating directories
    ssh.close()
    print("All directories created successfully.")

# Function to execute a command via SSH and return output and error
def execute_command_via_ssh(client, command):
    # Execute command using the SSH client
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read()
    error = stderr.read()
    return output.decode(), error.decode()

# Function to execute MUSCLE alignment on CDS files on the server
def execute_muscle_on_cds(server, port, user, key_file):
    base_directory = "/home/mfreeman/USCServer/sequences/CDS/Phocoena/"
    aligned_directory = f"{base_directory}CDS_protein/aligned/"

    # Create SSH client to connect to the server
    client = create_ssh_client(server, port, user, key_file)

    try:
        # List the files in the base directory and its subdirectories that match the '*.fasta' pattern
        output, error = execute_command_via_ssh(client, f"find {base_directory}CDS_protein -type f -name '*.fasta'")
        if error:
            print(f"Error retrieving list from {base_directory}: {error}")
            return
        files_list = [file for file in output.strip().split()]

        # Function to run MUSCLE alignment for a specific file
        def run_muscle_alignment(file):
            # Define the output path for the alignment result
            output_filename = os.path.basename(file).replace(".fasta", ".afa")
            output_path = os.path.join(aligned_directory, output_filename)
            # Check if the alignment already exists
            output, error = execute_command_via_ssh(client, f"test -f {output_path} && echo exists")
            if output.strip() == 'exists':
                print(f"Skipping {file}, already aligned.")
                return

            # Run MUSCLE alignment using the Singularity container
            command = f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {file} -out {output_path}"
            print(f"Aligning {file} using MUSCLE...")
            output, error = execute_command_via_ssh(client, command)
            if error:
                print(f"Failed to align {file}: {error}")

        # Run MUSCLE alignment for each file in parallel using ThreadPoolExecutor
        with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor:
            futures = [executor.submit(run_muscle_alignment, file) for file in files_list]
            concurrent.futures.wait(futures)

    finally:
        # Close the SSH client after completing all alignments
        client.close()

    print("Alignment process completed for all specified files.")

# Function to upload protein CDS files from local to server using SFTP
def upload_files(ssh_client, local_directory, remote_directory):
    sftp = ssh_client.open_sftp()
    try:
        # Walk through the local directory to find and upload protein CDS files
        for root, _, files in os.walk(local_directory):
            for file in files:
                # Only upload .fasta files that are part of the Phocoena CDS_protein directory
                if file.endswith('.fasta') and "Phocoena" in root and "CDS_protein" in root:
                    local_file_path = os.path.join(root, file)
                    relative_path = os.path.relpath(local_file_path, local_directory)
                    remote_file_path = os.path.join(remote_directory, relative_path).replace("\\", "/")
                    remote_dir = os.path.dirname(remote_file_path)

                    # Ensure the remote directory exists
                    ssh_client.exec_command(f'mkdir -p {remote_dir}')

                    # Upload the file
                    print(f"Uploading {local_file_path} to {remote_file_path}")
                    sftp.put(local_file_path, remote_file_path)
    except Exception as e:
        print(f"Error while uploading files: {e}")
    finally:
        # Close the SFTP session
        sftp.close()

# Function to download aligned files from server to local
def download_files(ssh_client, base_directory, local_directory):
    sftp = ssh_client.open_sftp()
    try:
        species_dir = "Phocoena"
        # Define local and remote directories for aligned files
        local_species_dir = os.path.join(local_directory, species_dir, "CDS_protein", "aligned")
        remote_aligned_dir = f"{base_directory}{species_dir}/CDS_protein/aligned"
        os.makedirs(local_species_dir, exist_ok=True)

        try:
            # List all .afa files in the remote aligned directory
            files = sftp.listdir(remote_aligned_dir)
            afa_files = [file for file in files if file.endswith('.afa')]
            if not afa_files:
                print(f"No .afa files to download in {remote_aligned_dir}")
                return
            print(f"Attempting to download .afa files from {remote_aligned_dir}: {afa_files}")
        except IOError as e:
            print(f"Failed to list files in {remote_aligned_dir}: {str(e)}")
            return

        for file in afa_files:
            # Define remote and local paths for each .afa file
            remote_file_path = os.path.join(remote_aligned_dir, file).replace("\\", "/")
            local_file_path = os.path.join(local_species_dir, file)
            temp_file_path = os.path.join(local_species_dir, "temp_" + file)
            print(f"Preparing to download: {remote_file_path} to {temp_file_path}")

            if os.path.exists(local_file_path):
                # Skip the download if the file already exists locally
                print(f"{file} already exists and won't be downloaded.")
            else:
                try:
                    # Download the file to a temporary path
                    sftp.get(remote_file_path, temp_file_path)
                    # Verify that the downloaded file is not empty
                    if os.path.getsize(temp_file_path) > 0:
                        # Move the temporary file to the final location
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
        # Close the SFTP session
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

# Upload CDS protein files from local to server
upload_files(client, local_dir, remote_dir)

# Execute MUSCLE alignment on CDS files on the server
execute_muscle_on_cds(server, port, user, key_file)

# Download CDS alignment files
client = create_ssh_client(server, port, user, key_file)
download_files(client, remote_dir, local_dir)

# Close SSH client
client.close()

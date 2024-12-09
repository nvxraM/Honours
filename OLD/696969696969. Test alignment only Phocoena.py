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

# Function to execute a command via SSH and return output and error
def execute_command_via_ssh(client, command):
    # Execute command using the SSH client
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read()
    error = stderr.read()
    return output.decode(), error.decode()

# Function to create necessary directories on the server
def create_directories_on_server(server, port, user, key_file, remote_directories):
    client = create_ssh_client(server, port, user, key_file)
    try:
        for directory in remote_directories:
            command = f'mkdir -p {directory}'
            output, error = execute_command_via_ssh(client, command)
            if error:
                print(f"Error creating directory {directory}: {error}")
            else:
                print(f"Directory created or already exists: {directory}")
    finally:
        client.close()

# Function to upload the specified file from local to server
def upload_file(server, port, user, key_file, local_file_path, remote_file_path):
    # Check if the local file exists before attempting to upload
    absolute_local_path = os.path.abspath(local_file_path)
    print(f"Checking if local file exists at: {absolute_local_path}")
    if not os.path.exists(absolute_local_path):
        print(f"Error: Local file does not exist: {absolute_local_path}")
        return

    # Create SSH client to connect to the server
    client = create_ssh_client(server, port, user, key_file)
    sftp = client.open_sftp()
    try:
        print(f"Uploading {absolute_local_path} to {remote_file_path}")
        sftp.put(absolute_local_path, remote_file_path)
        print(f"Upload completed: {absolute_local_path} to {remote_file_path}")

        # List the server directory to verify if the file is present
        remote_dir = os.path.dirname(remote_file_path)
        try:
            files = sftp.listdir(remote_dir)
            if os.path.basename(remote_file_path) in files:
                print(f"Verified: File exists on server: {remote_file_path}")
            else:
                print(f"Error: Uploaded file not found in directory: {remote_dir}")
        except IOError as e:
            print(f"Failed to list directory {remote_dir}: {str(e)}")

    except Exception as e:
        print(f"Failed to upload {absolute_local_path}: {str(e)}")
    finally:
        # Close the SFTP session and SSH client
        sftp.close()
        client.close()

# Function to execute MUSCLE alignment on the specified CDS file on the server
def execute_muscle_on_cds(server, port, user, key_file):
    input_file = "/home/mfreeman/USCServer/sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.fasta"
    output_file = "/home/mfreeman/USCServer/sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.afa"

    # Create SSH client to connect to the server
    client = create_ssh_client(server, port, user, key_file)

    try:
        # Check if the alignment already exists
        output, error = execute_command_via_ssh(client, f"test -f {output_file} && echo exists")
        if output.strip() == 'exists':
            print(f"Skipping alignment, already exists: {output_file}")
            return

        # Run MUSCLE alignment using the Singularity container
        command = f"singularity exec /RDS/Q1233/singularity/muscle.sif muscle -in {input_file} -out {output_file}"
        print(f"Aligning {input_file} using MUSCLE...")
        output, error = execute_command_via_ssh(client, command)
        if error:
            print(f"Failed to align {input_file}: {error}")
        else:
            print(f"Alignment completed: {output_file}")

    finally:
        # Close the SSH client after completing the alignment
        client.close()

    print("Alignment process completed for the specified file.")

# Function to download the aligned file from server to local
def download_file(ssh_client, remote_file_path, local_file_path):
    sftp = ssh_client.open_sftp()
    try:
        temp_file_path = local_file_path + ".temp"
        print(f"Preparing to download: {remote_file_path} to {temp_file_path}")

        if os.path.exists(local_file_path):
            # Skip the download if the file already exists locally
            print(f"{local_file_path} already exists and won't be downloaded.")
        else:
            try:
                # Download the file to a temporary path
                sftp.get(remote_file_path, temp_file_path)
                # Verify that the downloaded file is not empty
                if os.path.getsize(temp_file_path) > 0:
                    # Move the temporary file to the final location
                    shutil.move(temp_file_path, local_file_path)
                    print(f"Downloaded and verified {local_file_path}")
                else:
                    print(f"Downloaded file is empty, check server content.")
                    os.remove(temp_file_path)
            except Exception as e:
                print(f"Failed to download {remote_file_path}: {str(e)}")
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
local_file_to_upload = "sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.fasta"
remote_file_to_upload = "/home/mfreeman/USCServer/sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.fasta"
remote_file = "/home/mfreeman/USCServer/sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.afa"
local_file = "sequences/CDS/Phocoena/CDS_nucleotide/Phocoena_concatenated_sequences_nucleotide.afa"

# Create necessary directories on the server
remote_directories = [
    "/home/mfreeman/USCServer/sequences/CDS/Phocoena/CDS_nucleotide"
]
create_directories_on_server(server, port, user, key_file, remote_directories)

# Upload the CDS nucleotide file from local to server
upload_file(server, port, user, key_file, local_file_to_upload, remote_file_to_upload)

# Execute MUSCLE alignment on the specified CDS file on the server
execute_muscle_on_cds(server, port, user, key_file)

# Create SSH client
client = create_ssh_client(server, port, user, key_file)

# Download the aligned CDS file
download_file(client, remote_file, local_file)

# Close SSH client
client.close()

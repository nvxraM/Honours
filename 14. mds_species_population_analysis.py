import paramiko
import subprocess
import os
import time
import concurrent.futures
from scp import SCPClient

########################################
# Utility Functions
########################################

def run_command(command):
    """
    Run a shell command locally and return (output, error, returncode).
    """
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True
    )
    output, error = process.communicate()
    return output, error, process.returncode

def create_ssh_client(server, port, user, key_file):
    """
    Create and return an SSH client connected to the remote server.
    """
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, username=user, key_filename=key_file)
    return client

########################################
# FASTA Parsing and Conversion Functions
########################################

def extract_accession_genus_species(header_line):
    """
    Extract accession, genus, and species from a header line.

    Expected format: Accession_Genus_species
    Example: NC_005268.1_Balaena_mysticetus
    """
    parts = header_line.split('_')
    if len(parts) < 3:
        # Handle unexpected format gracefully
        accession = parts[0] if len(parts) > 0 else "UnknownAcc"
        genus = parts[1] if len(parts) > 1 else "UnknownGenus"
        species = "UnknownSpecies"
        return accession, genus, species

    genus = parts[-2]
    species = parts[-1]
    accession = "_".join(parts[:-2])
    return accession, genus, species

def parse_fasta_headers(fasta_path):
    """
    Parse a concatenated FASTA file to split sequences by species.

    Input format:
    >Accession_Genus_species
    """
    species_data = {}
    current_header = None
    current_seq = []

    if not os.path.isfile(fasta_path):
        print(f"FASTA file does not exist: {fasta_path}")
        return species_data

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Store previous record if exists
                if current_header and current_seq:
                    accession, genus, species = extract_accession_genus_species(current_header[1:])
                    species_name = f"{genus}_{species}"

                    if species_name not in species_data:
                        species_data[species_name] = {"headers": [], "sequences": []}
                    species_data[species_name]["headers"].append(current_header)
                    species_data[species_name]["sequences"].append("".join(current_seq))

                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Handle the last record
        if current_header and current_seq:
            accession, genus, species = extract_accession_genus_species(current_header[1:])
            species_name = f"{genus}_{species}"
            if species_name not in species_data:
                species_data[species_name] = {"headers": [], "sequences": []}
            species_data[species_name]["headers"].append(current_header)
            species_data[species_name]["sequences"].append("".join(current_seq))

    return species_data

def read_fasta(fasta_file):
    """
    Read a FASTA file and return (headers, sequences).
    """
    sequences = []
    headers = []
    with open(fasta_file, "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
                headers.append(line.strip()[1:])
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    return headers, sequences

def fasta_to_vcf(input_fasta, output_vcf):
    """
    Convert a FASTA file to VCF format, using the first sequence as reference.
    Only biallelic SNPs are included.
    """
    headers, sequences = read_fasta(input_fasta)
    if not sequences:
        raise ValueError("No sequences found in the input FASTA file.")

    reference_sequence = sequences[0]
    sequence_length = min(len(seq) for seq in sequences)

    with open(output_vcf, "w") as vcf_file:
        # VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write("##source=CustomFastaToVCF\n")
        vcf_file.write(f"##reference={headers[0]}\n")
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(headers) + '\n')

        # Iterate over positions
        for pos in range(sequence_length):
            column_bases = [seq[pos] for seq in sequences if len(seq) > pos]
            unique_bases = list(set(column_bases))

            # Keep positions with 2 unique valid alleles (A,T,C,G)
            if len(unique_bases) == 2 and all(base in "ATCG" for base in unique_bases):
                ref = reference_sequence[pos]
                alt = unique_bases[1] if unique_bases[0] == ref else unique_bases[0]

                genotypes = []
                for base in column_bases:
                    if base == ref:
                        genotypes.append("0/0")
                    elif base == alt:
                        genotypes.append("1/1")
                    else:
                        genotypes.append("./.")

                vcf_file.write(f"1\t{pos+1}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t" + "\t".join(genotypes) + "\n")

    print(f"VCF file created: {output_vcf}")

########################################
# Processing Species FASTA and Converting to VCF
########################################

def write_species_fasta_files(genus_dir, genus):
    """
    For a given genus, split the genus-level concatenated gapped FASTA into species FASTA files.
    """
    input_fasta = os.path.join(genus_dir, "CDS_nucleotide_gapped", f"{genus}_concatenated_gapped_sequences.fasta")
    if not os.path.exists(input_fasta):
        print(f"No input FASTA found for genus {genus}: {input_fasta}")
        return

    species_data = parse_fasta_headers(input_fasta)

    # Create MDS directory if not exist
    mds_dir = os.path.join(genus_dir, "CDS_nucleotide_gapped", "MDS")
    if not os.path.exists(mds_dir):
        os.makedirs(mds_dir)

    # Write each species FASTA
    for species_name, data in species_data.items():
        species_subdir = os.path.join(mds_dir, species_name)
        if not os.path.exists(species_subdir):
            os.makedirs(species_subdir)

        output_fasta = os.path.join(species_subdir, f"{species_name}_concatenated_gapped_sequences.fasta")
        with open(output_fasta, 'w') as out_f:
            for header, seq in zip(data["headers"], data["sequences"]):
                out_f.write(header + "\n")
                out_f.write(seq + "\n")

def convert_all_fasta_to_vcf(local_base):
    """
    Convert all species FASTA files in MDS/[Species] directories into VCF files.
    Input FASTA: [Species]_concatenated_gapped_sequences.fasta
    Output VCF: [Species]_concatenated_sequences.vcf
    """
    genus_dirs = [g for g in os.listdir(local_base) if os.path.isdir(os.path.join(local_base, g))]
    for genus in genus_dirs:
        genus_mds_dir = os.path.join(local_base, genus, "CDS_nucleotide_gapped", "MDS")
        if not os.path.exists(genus_mds_dir):
            continue

        species_dirs = [s for s in os.listdir(genus_mds_dir) if os.path.isdir(os.path.join(genus_mds_dir, s))]
        for species_name in species_dirs:
            species_dir = os.path.join(genus_mds_dir, species_name)
            fasta_path = os.path.join(species_dir, f"{species_name}_concatenated_gapped_sequences.fasta")
            if not os.path.exists(fasta_path):
                continue

            # VCF file name
            output_vcf = os.path.join(species_dir, f"{species_name}_concatenated_sequences.vcf")

            if os.path.exists(output_vcf):
                print(f"VCF file already exists for {species_name}, skipping conversion.")
                continue

            try:
                fasta_to_vcf(fasta_path, output_vcf)
            except ValueError as e:
                print(f"Error processing {fasta_path}: {e}")

########################################
# File Transfer and MDS Analysis
########################################

def upload_vcf_files_to_server():
    """
    Upload VCF files to the server at ~/USCServer/MDS/[Genus]/[Species].
    """
    local_base = "sequences/CDS_Genus"
    remote_base = "~/USCServer/MDS"
    key_path = r"C:\Users\freem\Downloads\mfreeman-private-key.txt"

    genus_dirs = [g for g in os.listdir(local_base) if os.path.isdir(os.path.join(local_base, g))]

    for genus in genus_dirs:
        genus_mds_dir = os.path.join(local_base, genus, "CDS_nucleotide_gapped", "MDS")
        if not os.path.exists(genus_mds_dir):
            continue

        species_dirs = [s for s in os.listdir(genus_mds_dir) if os.path.isdir(os.path.join(genus_mds_dir, s))]
        for species_name in species_dirs:
            species_dir = os.path.join(genus_mds_dir, species_name)
            local_vcf_path = os.path.join(species_dir, f"{species_name}_concatenated_sequences.vcf")
            if not os.path.isfile(local_vcf_path):
                continue

            remote_genus_dir = f"{remote_base}/{genus}"
            remote_species_dir = f"{remote_genus_dir}/{species_name}"
            remote_vcf_path = f"{remote_species_dir}/{species_name}_concatenated_sequences.vcf"

            try:
                ssh = create_ssh_client("203.101.229.234", 22, "mfreeman", key_path)
                scp = SCPClient(ssh.get_transport())

                # Create directories on server
                stdin, stdout, stderr = ssh.exec_command(f"mkdir -p {remote_species_dir}")
                exit_status = stdout.channel.recv_exit_status()
                if exit_status != 0:
                    error_msg = stderr.read().decode()
                    print(f"Failed to create directory {remote_species_dir}: {error_msg}")
                    scp.close()
                    ssh.close()
                    continue

                print(f"Uploading {local_vcf_path} to {remote_vcf_path}...")
                scp.put(local_vcf_path, remote_vcf_path)
                print("Upload complete.")

                scp.close()
                ssh.close()

            except Exception as e:
                print(f"Error uploading {local_vcf_path}: {str(e)}")

def run_mds_analysis():
    """
    Run MDS analysis on the server using plink commands on the uploaded VCF files.
    We add '--double-id' to the first plink command to handle underscores in IDs.
    """
    key_path = r"C:\Users\freem\Downloads\mfreeman-private-key.txt"
    local_base = "sequences/CDS_Genus"
    remote_base = "~/USCServer/MDS"

    print("Establishing SSH connection for MDS analysis...")
    ssh = create_ssh_client("203.101.229.234", 22, "mfreeman", key_path)
    print("SSH connection established.")

    genus_dirs = [g for g in os.listdir(local_base) if os.path.isdir(os.path.join(local_base, g))]
    genus_species_pairs = []
    for genus in genus_dirs:
        genus_mds_dir = os.path.join(local_base, genus, "CDS_nucleotide_gapped", "MDS")
        if not os.path.exists(genus_mds_dir):
            continue
        species_dirs = [s for s in os.listdir(genus_mds_dir) if os.path.isdir(os.path.join(genus_mds_dir, s))]
        for species_name in species_dirs:
            genus_species_pairs.append((genus, species_name))

    def run_species_mds(pair):
        genus, species_name = pair
        remote_dir = f"{remote_base}/{genus}/{species_name}"

        # Ensure directory
        mkdir_cmd = f"mkdir -p {remote_dir}"
        stdin, stdout, stderr = ssh.exec_command(mkdir_cmd)
        exit_status = stdout.channel.recv_exit_status()
        if exit_status != 0:
            error_msg = stderr.read().decode()
            print(f"Failed to create directory {remote_dir}: {error_msg}")
            return

        vcf_file = f"{remote_dir}/{species_name}_concatenated_sequences.vcf"
        filtered_prefix = f"{remote_dir}/{species_name}_filtered"
        filtered_mind_geno_prefix = f"{remote_dir}/{species_name}_filtered_mind_geno"
        genome_prefix = f"{remote_dir}/{species_name}_genome"
        mds_prefix = f"{remote_dir}/{species_name}_mds"

        # Add '--double-id' to handle underscores in sample IDs
        commands = [
            f"singularity exec /RDS/Q1233/singularity/plink1-9.sif plink --vcf {vcf_file} --double-id --make-bed --out {filtered_prefix}",
            f"singularity exec /RDS/Q1233/singularity/plink1-9.sif plink --bfile {filtered_prefix} --mind 0.2 --geno 0.1 --make-bed --out {filtered_mind_geno_prefix}",
            f"singularity exec /RDS/Q1233/singularity/plink1-9.sif plink --bfile {filtered_mind_geno_prefix} --genome --out {genome_prefix}",
            f"singularity exec /RDS/Q1233/singularity/plink1-9.sif plink --bfile {filtered_mind_geno_prefix} --read-genome {genome_prefix}.genome --cluster --mds-plot 10 --out {mds_prefix}"
        ]

        # Execute each command
        for command in commands:
            print(f"Running command for {genus}/{species_name}: {command}")
            stdin, stdout, stderr = ssh.exec_command(command)

            # 5-minute timeout
            start_time = time.time()
            while not stdout.channel.exit_status_ready():
                if time.time() - start_time > 300:
                    print(f"Timeout reached for command: {command}")
                    return
                time.sleep(5)

            stdout_output = stdout.read().decode()
            stderr_output = stderr.read().decode()
            exit_status = stdout.channel.recv_exit_status()

            if exit_status == 0:
                print(f"Command succeeded for {genus}/{species_name}: {command}")
            else:
                print(f"Command failed for {genus}/{species_name}: {command}\n{stderr_output}")
                with open(f"error_log_{species_name}.txt", "w") as error_log:
                    error_log.write(f"Command: {command}\n")
                    error_log.write(f"Error: {stderr_output}\n")
                return

        # Check if MDS file exists
        mds_file = f"{mds_prefix}.mds"
        check_cmd = f"test -f {mds_file} && echo 'exists' || echo 'not_exists'"
        stdin, stdout, stderr = ssh.exec_command(check_cmd)
        output = stdout.read().decode().strip()
        if output == "exists":
            print(f"MDS file for {genus}/{species_name} successfully created at {mds_file}")
        else:
            print(f"Warning: MDS file for {genus}/{species_name} does not exist at {mds_file}")

    # Run MDS analysis in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        executor.map(run_species_mds, genus_species_pairs)

    ssh.close()
    print("SSH connection closed after MDS analysis.")

def download_mds_files_from_server():
    """
    Download MDS files (.mds) from the server back to local machine.
    """
    local_base = "sequences/CDS_Genus"
    remote_base = "~/USCServer/MDS"
    key_path = r"C:\Users\freem\Downloads\mfreeman-private-key.txt"

    genus_dirs = [g for g in os.listdir(local_base) if os.path.isdir(os.path.join(local_base, g))]
    for genus in genus_dirs:
        genus_mds_dir = os.path.join(local_base, genus, "CDS_nucleotide_gapped", "MDS")
        if not os.path.exists(genus_mds_dir):
            continue
        species_dirs = [s for s in os.listdir(genus_mds_dir) if os.path.isdir(os.path.join(genus_mds_dir, s))]
        for species_name in species_dirs:
            local_mds_path = os.path.join(genus_mds_dir, species_name, f"{species_name}_mds.mds")
            remote_mds_file = f"{remote_base}/{genus}/{species_name}/{species_name}_mds.mds"

            try:
                ssh = create_ssh_client("203.101.229.234", 22, "mfreeman", key_path)
                scp = SCPClient(ssh.get_transport())

                if not os.path.exists(os.path.dirname(local_mds_path)):
                    os.makedirs(os.path.dirname(local_mds_path))

                print(f"Downloading {remote_mds_file} to {local_mds_path}...")
                scp.get(remote_mds_file, local_mds_path)
                print("Download complete.")

                scp.close()
                ssh.close()

            except Exception as e:
                print(f"Error downloading {remote_mds_file}: {str(e)}")

########################################
# Main Execution
########################################

if __name__ == "__main__":
    local_base = "sequences/CDS_Genus"

    # Step 1: Split concatenated FASTA files into per-species FASTA files
    genus_list = [g for g in os.listdir(local_base) if os.path.isdir(os.path.join(local_base, g))]
    for genus in genus_list:
        genus_dir = os.path.join(local_base, genus)
        write_species_fasta_files(genus_dir, genus)

    # Step 2: Convert species FASTA to VCF
    convert_all_fasta_to_vcf(local_base)

    # Step 3: Upload VCF files to the server
    upload_vcf_files_to_server()

    # Step 4: Run MDS analysis on the server with --double-id
    run_mds_analysis()

    # Step 5: Download MDS results from the server
    download_mds_files_from_server()

    print("All steps completed.")

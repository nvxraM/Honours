import os
from pathlib import Path
from Bio import Entrez, SeqIO
import shelve
from tqdm import tqdm
import time
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Rate limit lock for controlling requests per second
rate_limit_lock = threading.Lock()
requests_per_second = 10
last_request_time = [0.0]  # Use a list to allow modification within threads

# Function to rate limit requests to NCBI
def rate_limited_request():
    with rate_limit_lock:
        current_time = time.time()
        elapsed = current_time - last_request_time[0]
        if elapsed < 1.0 / requests_per_second:
            time.sleep((1.0 / requests_per_second) - elapsed)
        last_request_time[0] = time.time()

# Function to read Entrez credentials from a file
def read_entrez_credentials(filepath):
    credentials = {}
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split(' = ')
            credentials[key] = value.strip('"')
    return credentials['email'], credentials['api_key']

# Function to setup Entrez configuration and count the available records for the given query
def setup(entrez_email, entrez_api_key, entrez_search_query, entrez_database):
    Entrez.email = entrez_email
    Entrez.api_key = entrez_api_key
    logging.info(f"Setup: Configuring NCBI Entrez with email: {entrez_email}")
    logging.info(f"Setup: Counting available records for query: '{entrez_search_query}' on NCBI...")
    rate_limited_request()
    with Entrez.esearch(db=entrez_database, term=entrez_search_query, retmax=0) as handle:
        total_records_count = int(Entrez.read(handle)['Count'])
    logging.info(f"Setup: Total records available on NCBI: {total_records_count}")
    return total_records_count

# Utility function to format organism name
def format_organism_name(organism_name):
    return organism_name.replace(" ", "_") if organism_name else 'Unknown'

# Utility function to execute a command with retry logic
def execute_with_retry(command, retries, wait=5):
    attempt = 0
    while attempt < retries:
        try:
            rate_limited_request()
            return command()
        except (RuntimeError, OSError, ValueError) as e:
            attempt += 1
            logging.warning(f"Error occurred: {e}. Retrying ({attempt}/{retries})...")
            time.sleep(wait)
    raise Exception("Failed to complete operation after multiple attempts.")

# Function to download a single batch of GenBank files
def download_gb_batch(start, batch_size, entrez_database, entrez_search_query, gb_directory):
    try:
        handle = execute_with_retry(lambda: Entrez.esearch(db=entrez_database, term=entrez_search_query, retstart=start, retmax=batch_size), 3)
        id_list = Entrez.read(handle)['IdList']

        if not id_list:
            logging.warning(f"No IDs found for batch starting at {start}")
            return

        fetch_handle = execute_with_retry(lambda: Entrez.efetch(db=entrez_database, id=",".join(id_list), rettype="gb", retmode="text"), 3)
        records = SeqIO.parse(fetch_handle, "gb")
        records_to_write = []

        for record in records:
            if not record.id:
                logging.warning(f"Record without ID found in batch starting at {start}")
                continue
            organism = format_organism_name(record.annotations.get('organism', 'Unknown'))
            organism_dir = os.path.join(gb_directory, organism)
            os.makedirs(organism_dir, exist_ok=True)
            filename = os.path.join(organism_dir, f"{record.id}.gb")

            if not os.path.exists(filename):
                records_to_write.append((filename, record))

        for filename, record in records_to_write:
            SeqIO.write([record], filename, "gb")

    except Exception as e:
        logging.warning(f"Failed to download batch starting at {start}. Error: {e}")

# Main function to parallelize downloading of GenBank files
def download_gb_files_parallel(email, api_key, search_query, database, gb_directory, total_records, batch_size=500):
    with ThreadPoolExecutor(max_workers=5) as executor:  # Limit to 5 workers to stay within 10 rps
        futures = [
            executor.submit(download_gb_batch, start, batch_size, database, search_query, gb_directory)
            for start in range(0, total_records, batch_size)
        ]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading GenBank files", unit="batch"):
            pass

# Function to create a mapping of sequence IDs to organism names from downloaded GenBank files
def refresh_id_to_organism_map(gb_directory):
    logging.info(f"Mapping: Refreshing ID to organism map from files in {gb_directory}")
    id_to_organism = {}
    for path in Path(gb_directory).rglob('*.gb'):
        with open(path, 'r') as gb_file:
            for record in SeqIO.parse(gb_file, 'genbank'):
                if not record.id:
                    logging.warning(f"Record without ID found in file: {path}")
                    continue
                organism_name = format_organism_name(record.annotations.get('organism', 'Unknown'))
                base_id = record.id.split('.')[0]
                id_to_organism[base_id] = organism_name
    logging.info(f"Mapping: Refreshed map with {len(id_to_organism)} entries")
    return id_to_organism

# Function to download FASTA files in batches
def download_fasta_files(entrez_email, entrez_api_key, entrez_search_query, entrez_database, fasta_directory, total_records, id_to_organism, batch_size=500):
    logging.info(f"Downloading FASTA: Checking and downloading FASTA files to {fasta_directory}")
    os.makedirs(fasta_directory, exist_ok=True)
    fasta_db = None
    try:
        fasta_db = shelve.open('fasta_records_db', writeback=True)
        for start in tqdm(range(0, total_records, batch_size), desc="Downloading FASTA files", unit="batch"):
            try:
                handle = execute_with_retry(lambda: Entrez.esearch(db=entrez_database, term=entrez_search_query, retstart=start, retmax=batch_size), 3)
                id_list = Entrez.read(handle)['IdList']

                if not id_list:
                    logging.warning(f"No IDs found for batch starting at {start}")
                    continue

                fetch_handle = execute_with_retry(lambda: Entrez.efetch(db=entrez_database, id=",".join(id_list), rettype="fasta", retmode="text"), 3)
                records = SeqIO.parse(fetch_handle, "fasta")

                records_to_write = []
                for record in records:
                    if not record.id:
                        logging.warning(f"Record without ID found in batch starting at {start}")
                        continue
                    base_id = record.id.split('.')[0]
                    organism = id_to_organism.get(base_id, 'Unknown')
                    organism_dir = os.path.join(fasta_directory, organism)
                    os.makedirs(organism_dir, exist_ok=True)
                    organism_file = os.path.join(organism_dir, f"{organism}.fasta")

                    if base_id not in fasta_db:
                        records_to_write.append((organism_file, record))
                        fasta_db[base_id] = True

                for organism_file, record in records_to_write:
                    with open(organism_file, 'a') as f:
                        SeqIO.write([record], f, "fasta")
            except Exception as e:
                tqdm.write(f"Failed to download batch starting at {start}. Error: {e}")
                continue
    except Exception as e:
        logging.error(f"An error occurred while accessing the shelve database: {e}")
    finally:
        if fasta_db is not None:
            fasta_db.close()
            # Delete shelve files if they exist
            for ext in ['bak', 'dat', 'dir']:
                file_path = f'fasta_records_db.{ext}'
                if os.path.exists(file_path):
                    os.remove(file_path)

# Utility function to count the number of sequence records in the given directory with the specified extension
def count_sequence_records(directory, extension):
    total_sequences = len(list(Path(directory).rglob(f'*{extension}')))
    return total_sequences

# Main function to initiate the download of GenBank and FASTA files
if __name__ == "__main__":
    credentials_file = r"C:\Users\freem\OneDrive\Documents\USC\Honours\API keys\Entrez.txt"
    email, api_key = read_entrez_credentials(credentials_file)
    search_query = "(Cetacea[Organism] OR (Cetacea[Organism] OR Cetacea[All Fields])) AND Mitochondrion[All Fields] "
    search_query += "AND (complete[All Fields] AND genome[All Fields]) AND (14000[SLEN] : 18000[SLEN])"
    database = "nucleotide"
    base_directory = "sequences/"
    total_records = setup(email, api_key, search_query, database)
    if total_records > 0:
        gb_directory = os.path.join(base_directory, 'gb')
        fasta_directory = os.path.join(base_directory, 'fasta')
        download_gb_files_parallel(email, api_key, search_query, database, gb_directory, total_records)
        id_to_organism = refresh_id_to_organism_map(gb_directory)
        download_fasta_files(email, api_key, search_query, database, fasta_directory, total_records, id_to_organism)

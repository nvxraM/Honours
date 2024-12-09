import os
import time
import logging
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import Entrez, SeqIO
import shelve
from tqdm import tqdm

# Configure logging for detailed runtime information
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Global variables for rate limiting
rate_limit_lock = threading.Lock()
requests_per_second = 10
last_request_time = [0.0]  # list used as mutable to store last request time

def rate_limited_request():
    """
    Enforce a rate limit on requests to NCBI Entrez servers.
    Ensures that no more than 'requests_per_second' requests are sent.
    """
    with rate_limit_lock:
        current_time = time.time()
        elapsed = current_time - last_request_time[0]
        # If too soon since last request, wait until next permissible request time
        if elapsed < 1.0 / requests_per_second:
            time.sleep((1.0 / requests_per_second) - elapsed)
        last_request_time[0] = time.time()

def read_entrez_credentials(filepath):
    """
    Read Entrez API credentials (email, api_key) from a file.
    The file is expected to have lines in the format: key = "value".
    """
    credentials = {}
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            key, value = line.strip().split(' = ')
            credentials[key] = value.strip('"')
    return credentials['email'], credentials['api_key']

def setup(entrez_email, entrez_api_key, entrez_search_query, entrez_database):
    """
    Set up Entrez parameters and obtain the total count of records matching the search query.
    """
    Entrez.email = entrez_email
    Entrez.api_key = entrez_api_key
    logging.info(f"Setup: Counting available records for query: '{entrez_search_query}'...")
    rate_limited_request()
    with Entrez.esearch(db=entrez_database, term=entrez_search_query, retmax=0) as handle:
        total_records_count = int(Entrez.read(handle)['Count'])
    logging.info(f"Setup: Total records available on NCBI: {total_records_count}")
    return total_records_count

def format_organism_name(organism_name):
    """
    Format an organism name by replacing spaces with underscores.
    Returns 'Unknown' if organism_name is None.
    """
    return organism_name.replace(" ", "_") if organism_name else 'Unknown'

def execute_with_retry(command, retries, wait=5):
    """
    Execute a given command (usually an Entrez request) with retries.
    If an error occurs, it will retry up to 'retries' times, waiting 'wait' seconds between attempts.
    """
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

def get_all_ids(entrez_database, entrez_search_query, total_records):
    """
    Fetch all matching IDs from NCBI using ESearch in a paginated manner.
    This helps when there are more records than can be returned in a single call.
    """
    all_ids = []
    batch_size = 10000  # Large batch size reduces the number of API calls
    for start in range(0, total_records, batch_size):
        handle = execute_with_retry(
            lambda: Entrez.esearch(db=entrez_database,
                                   term=entrez_search_query,
                                   retstart=start,
                                   retmax=batch_size),
            3
        )
        result = Entrez.read(handle)
        batch_ids = result['IdList']
        all_ids.extend(batch_ids)
    logging.info(f"Retrieved {len(all_ids)} IDs in total.")
    return all_ids

def post_ids(all_ids, entrez_database):
    """
    Use EPost to upload the list of IDs to the NCBI server.
    This returns a WebEnv and a query_key that can be used with EFetch without re-specifying IDs.
    """
    id_str = ",".join(all_ids)
    handle = execute_with_retry(lambda: Entrez.epost(db=entrez_database, id=id_str), 3)
    result = Entrez.read(handle)
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    logging.info("Successfully posted IDs to NCBI (EPost).")
    return webenv, query_key

def download_gb_batch(start, batch_size, entrez_database, gb_directory, webenv, query_key):
    """
    Download a batch of GenBank records using EFetch.
    Writes each record to a GenBank file in a subdirectory named after its organism.
    """
    try:
        fetch_handle = execute_with_retry(
            lambda: Entrez.efetch(db=entrez_database,
                                  query_key=query_key,
                                  WebEnv=webenv,
                                  rettype="gb",
                                  retmode="text",
                                  retstart=start,
                                  retmax=batch_size),
            3
        )
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

        # Write all new records
        for filename, record in records_to_write:
            SeqIO.write([record], filename, "gb")

    except Exception as e:
        logging.warning(f"Failed to download batch starting at {start}. Error: {e}")

def download_gb_files_parallel(email, api_key, database, gb_directory, total_records, webenv, query_key, batch_size=500):
    """
    Download all GenBank files in parallel.
    Uses a ThreadPoolExecutor to fetch multiple batches simultaneously.
    """
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = [
            executor.submit(download_gb_batch, start, batch_size, database, gb_directory, webenv, query_key)
            for start in range(0, total_records, batch_size)
        ]
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Downloading GenBank files", unit="batch"):
            pass

def refresh_id_to_organism_map(gb_directory):
    """
    Re-build a map from sequence ID (base ID) to organism name by scanning all downloaded GenBank files.
    """
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
    logging.info(f"Mapping: Refreshed map with {len(id_to_organism)} entries.")
    return id_to_organism

def download_fasta_files(entrez_database, fasta_directory, total_records, id_to_organism, webenv, query_key, batch_size=500):
    """
    Download FASTA files corresponding to previously fetched GenBank records.
    Uses the ID->organism map to place sequences in organism-specific directories.
    Uses a shelve database to keep track of previously downloaded FASTA records and avoid duplicates.
    """
    logging.info(f"Downloading FASTA: Checking and downloading FASTA files to {fasta_directory}")
    os.makedirs(fasta_directory, exist_ok=True)

    fasta_db = None
    try:
        # Use a shelve database to track downloaded FASTA records
        fasta_db = shelve.open('fasta_records_db', writeback=True)
        for start in tqdm(range(0, total_records, batch_size), desc="Downloading FASTA files", unit="batch"):
            try:
                # Fetch FASTA sequences in batches using the stored WebEnv and query_key
                fetch_handle = execute_with_retry(
                    lambda: Entrez.efetch(db=entrez_database,
                                          query_key=query_key,
                                          WebEnv=webenv,
                                          rettype="fasta",
                                          retmode="text",
                                          retstart=start,
                                          retmax=batch_size),
                    3
                )
                records = SeqIO.parse(fetch_handle, "fasta")

                records_to_write = []
                for record in records:
                    if not record.id:
                        logging.warning(f"Record without ID found in FASTA batch starting at {start}")
                        continue
                    base_id = record.id.split('.')[0]
                    organism = id_to_organism.get(base_id, 'Unknown')
                    organism_dir = os.path.join(fasta_directory, organism)
                    os.makedirs(organism_dir, exist_ok=True)
                    organism_file = os.path.join(organism_dir, f"{organism}.fasta")

                    # Only write records not already in the database
                    if base_id not in fasta_db:
                        records_to_write.append((organism_file, record))
                        fasta_db[base_id] = True

                # Append new records to the respective organism file
                for organism_file, record in records_to_write:
                    with open(organism_file, 'a') as f:
                        SeqIO.write([record], f, "fasta")

            except Exception as e:
                tqdm.write(f"Failed to download batch starting at {start}. Error: {e}")
                continue

    except Exception as e:
        logging.error(f"An error occurred while accessing the shelve database: {e}")
    finally:
        # Close and clean up the shelve database
        if fasta_db is not None:
            fasta_db.close()
            # Remove shelve auxiliary files after use
            for ext in ['bak', 'dat', 'dir']:
                file_path = f'fasta_records_db.{ext}'
                if os.path.exists(file_path):
                    os.remove(file_path)

def count_sequence_records(directory, extension):
    """
    Count the number of sequence files (with given extension) in a directory.
    Useful for verification after download.
    """
    total_sequences = len(list(Path(directory).rglob(f'*{extension}')))
    return total_sequences

if __name__ == "__main__":
    # Path to the file containing Entrez credentials
    credentials_file = r"C:\Users\freem\OneDrive\Documents\USC\Honours\API keys\Entrez.txt"
    email, api_key = read_entrez_credentials(credentials_file)

    # Define search query and database
    search_query = ("(Cetacea[Organism] OR (Cetacea[Organism] OR Cetacea[All Fields])) "
                    "AND Mitochondrion[All Fields] AND (complete[All Fields] AND genome[All Fields]) "
                    "AND (14000[SLEN] : 18000[SLEN])")
    database = "nucleotide"
    base_directory = "sequences/"

    # 1. Setup and get total count of records
    total_records = setup(email, api_key, search_query, database)

    if total_records > 0:
        gb_directory = os.path.join(base_directory, 'gb')
        fasta_directory = os.path.join(base_directory, 'fasta')

        # 2. Retrieve all IDs
        all_ids = get_all_ids(database, search_query, total_records)

        # 3. Post IDs to get WebEnv and query_key
        webenv, query_key = post_ids(all_ids, database)

        # 4. Download GenBank records in parallel
        download_gb_files_parallel(email, api_key, database, gb_directory, total_records, webenv, query_key)

        # 5. Refresh ID -> organism map from downloaded GenBank files
        id_to_organism = refresh_id_to_organism_map(gb_directory)

        # 6. Download corresponding FASTA sequences for each ID
        download_fasta_files(database, fasta_directory, total_records, id_to_organism, webenv, query_key)

import os
import time
import logging
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio import Entrez, SeqIO
import shelve
from tqdm import tqdm

# Configure logging: INFO level shows general progress and warnings
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Global rate-limiting variables
rate_limit_lock = threading.Lock()  # Lock to prevent race conditions on request timing
requests_per_second = 10  # Maximum requests per second
last_request_time = [0.0]  # Using a list for mutability in threads


def rate_limited_request():
    """
    Ensures that requests to the NCBI Entrez API do not exceed the allowed rate.
    This function enforces a maximum number of requests per second by introducing a delay if needed.
    """
    with rate_limit_lock:
        current_time = time.time()
        elapsed = current_time - last_request_time[0]
        if elapsed < 1.0 / requests_per_second:
            # If we haven't waited enough since the last request, wait the required extra time
            time.sleep((1.0 / requests_per_second) - elapsed)
        last_request_time[0] = time.time()


def read_entrez_credentials(filepath):
    """
    Reads Entrez credentials from a file and returns the email and API key.

    Parameters:
        filepath (str): Path to the file containing 'email' and 'api_key'.

    Returns:
        (str, str): Email and API key for the Entrez API.
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
    Configures Entrez with the provided credentials and counts the total records available on NCBI
    that match the given query.

    Parameters:
        entrez_email (str): Email address for Entrez.
        entrez_api_key (str): API key for Entrez.
        entrez_search_query (str): NCBI Entrez search query.
        entrez_database (str): NCBI database name (e.g., "nucleotide").

    Returns:
        int: Total count of records available for the search query.
    """
    Entrez.email = entrez_email
    Entrez.api_key = entrez_api_key
    logging.info(f"Setup: Configuring NCBI Entrez with email: {entrez_email}")
    logging.info(f"Setup: Counting available records for query: '{entrez_search_query}' on NCBI...")

    # Rate limit the request
    rate_limited_request()

    # Perform the search and read the count of matching records
    with Entrez.esearch(db=entrez_database, term=entrez_search_query, retmax=0) as handle:
        total_records_count = int(Entrez.read(handle)['Count'])
    logging.info(f"Setup: Total records available on NCBI: {total_records_count}")
    return total_records_count


def format_organism_name(organism_name):
    """
    Replaces spaces with underscores in organism names for use as directory names.

    Parameters:
        organism_name (str): The organism's name.

    Returns:
        str: Formatted organism name.
    """
    return organism_name.replace(" ", "_") if organism_name else 'Unknown'


def execute_with_retry(command, retries, wait=5):
    """
    Executes a given command (function) with retry logic.

    Parameters:
        command (callable): Function or lambda to call.
        retries (int): Number of attempts before raising an error.
        wait (int): Seconds to wait between retries.

    Returns:
        The result of the command if successful.

    Raises:
        Exception: If the command fails after the given number of retries.
    """
    attempt = 0
    while attempt < retries:
        try:
            # Ensure we respect the request rate limit
            rate_limited_request()
            return command()
        except (RuntimeError, OSError, ValueError) as e:
            attempt += 1
            logging.warning(f"Error occurred: {e}. Retrying ({attempt}/{retries})...")
            time.sleep(wait)
    raise Exception("Failed to complete operation after multiple attempts.")


def download_gb_batch(start, batch_size, entrez_database, entrez_search_query, gb_directory):
    """
    Downloads a batch of GenBank (GB) records starting at 'start' index.

    Parameters:
        start (int): The starting record index for this batch.
        batch_size (int): How many records to fetch in this batch.
        entrez_database (str): NCBI database name.
        entrez_search_query (str): Query used to search for records.
        gb_directory (str): Directory to store downloaded GB files.
    """
    try:
        # Search for a batch of record IDs
        handle = execute_with_retry(lambda: Entrez.esearch(db=entrez_database, term=entrez_search_query,
                                                           retstart=start, retmax=batch_size), 3)
        id_list = Entrez.read(handle)['IdList']

        if not id_list:
            logging.warning(f"No IDs found for batch starting at {start}")
            return

        # Fetch the GenBank records
        fetch_handle = execute_with_retry(lambda: Entrez.efetch(db=entrez_database, id=",".join(id_list),
                                                                rettype="gb", retmode="text"), 3)
        records = SeqIO.parse(fetch_handle, "gb")

        # We'll only write out records that don't already exist
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

        # Write the new GenBank files
        for filename, record in records_to_write:
            SeqIO.write([record], filename, "gb")

    except Exception as e:
        logging.warning(f"Failed to download batch starting at {start}. Error: {e}")


def download_gb_files_parallel(email, api_key, search_query, database, gb_directory, total_records, batch_size=500):
    """
    Parallelizes the downloading of GenBank files using multiple threads.

    Parameters:
        email (str): Entrez email address.
        api_key (str): Entrez API key.
        search_query (str): NCBI search query.
        database (str): NCBI database name.
        gb_directory (str): Directory to store GenBank files.
        total_records (int): Total number of records to download.
        batch_size (int): Number of records per batch.
    """
    # Limit the number of workers to 5 because we have 10 requests/sec and want to avoid exceeding that.
    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = [
            executor.submit(download_gb_batch, start, batch_size, database, search_query, gb_directory)
            for start in range(0, total_records, batch_size)
        ]
        # as_completed allows us to track progress of all submitted tasks
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Downloading GenBank files", unit="batch"):
            pass


def refresh_id_to_organism_map(gb_directory):
    """
    Scans the downloaded GenBank files to create a mapping from sequence ID to organism name.

    Parameters:
        gb_directory (str): Directory containing GenBank files organized by organism.

    Returns:
        dict: A dictionary mapping {base_id: organism_name}.
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
    logging.info(f"Mapping: Refreshed map with {len(id_to_organism)} entries")
    return id_to_organism


def download_fasta_files(entrez_email, entrez_api_key, entrez_search_query, entrez_database, fasta_directory,
                         total_records, id_to_organism, batch_size=500):
    """
    Downloads FASTA files for each record previously retrieved, grouping them by organism. Uses a shelve
    database to keep track of processed IDs and avoid duplicates.

    Parameters:
        entrez_email (str): Entrez email address.
        entrez_api_key (str): Entrez API key.
        entrez_search_query (str): NCBI search query.
        entrez_database (str): NCBI database name.
        fasta_directory (str): Directory to store FASTA files.
        total_records (int): Total number of records to download.
        id_to_organism (dict): Mapping from base_id to organism name.
        batch_size (int): Number of records per batch.
    """
    logging.info(f"Downloading FASTA: Checking and downloading FASTA files to {fasta_directory}")
    os.makedirs(fasta_directory, exist_ok=True)

    fasta_db = None
    try:
        # Using a shelve database to remember which records have been processed
        fasta_db = shelve.open('fasta_records_db', writeback=True)

        # Iterate over the total records in batches
        for start in tqdm(range(0, total_records, batch_size), desc="Downloading FASTA files", unit="batch"):
            try:
                # Get the list of IDs for this batch
                handle = execute_with_retry(lambda: Entrez.esearch(db=entrez_database, term=entrez_search_query,
                                                                   retstart=start, retmax=batch_size), 3)
                id_list = Entrez.read(handle)['IdList']

                if not id_list:
                    logging.warning(f"No IDs found for batch starting at {start}")
                    continue

                # Fetch FASTA data for these IDs
                fetch_handle = execute_with_retry(lambda: Entrez.efetch(db=entrez_database, id=",".join(id_list),
                                                                        rettype="fasta", retmode="text"), 3)
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

                    # Only write if we haven't processed this ID before
                    if base_id not in fasta_db:
                        records_to_write.append((organism_file, record))
                        fasta_db[base_id] = True

                # Append the new records to the organism-specific FASTA file
                for organism_file, record in records_to_write:
                    with open(organism_file, 'a') as f:
                        SeqIO.write([record], f, "fasta")

            except Exception as e:
                tqdm.write(f"Failed to download batch starting at {start}. Error: {e}")
                continue

    except Exception as e:
        logging.error(f"An error occurred while accessing the shelve database: {e}")
    finally:
        # Close the shelve database
        if fasta_db is not None:
            fasta_db.close()
            # Clean up shelve auxiliary files
            for ext in ['bak', 'dat', 'dir']:
                file_path = f'fasta_records_db.{ext}'
                if os.path.exists(file_path):
                    os.remove(file_path)


def count_sequence_records(directory, extension):
    """
    Counts the number of sequence files in the given directory that match a certain extension.

    Parameters:
        directory (str): Directory to search.
        extension (str): File extension to match (e.g. '.gb', '.fasta').

    Returns:
        int: Count of matching files.
    """
    total_sequences = len(list(Path(directory).rglob(f'*{extension}')))
    return total_sequences


if __name__ == "__main__":
    # File containing Entrez credentials in the format:
    # email = "your_email@example.com"
    # api_key = "your_api_key"
    credentials_file = r"C:\Users\freem\OneDrive\Documents\USC\Honours\API keys\Entrez.txt"

    email, api_key = read_entrez_credentials(credentials_file)

    # This search query looks for Cetacea mitochondrial genomes within a certain length range
    search_query = "(Cetacea[Organism] OR (Cetacea[Organism] OR Cetacea[All Fields])) AND Mitochondrion[All Fields] "
    search_query += "AND (complete[All Fields] AND genome[All Fields]) AND (14000[SLEN] : 18000[SLEN])"
    database = "nucleotide"
    base_directory = "sequences/"

    # Count how many records are available for the given query
    total_records = setup(email, api_key, search_query, database)

    if total_records > 0:
        gb_directory = os.path.join(base_directory, 'gb')
        fasta_directory = os.path.join(base_directory, 'fasta')

        # Download GenBank files in parallel
        download_gb_files_parallel(email, api_key, search_query, database, gb_directory, total_records)

        # Create a map from ID to organism name based on GenBank annotations
        id_to_organism = refresh_id_to_organism_map(gb_directory)

        # Download FASTA files using the ID-to-organism map
        download_fasta_files(email, api_key, search_query, database, fasta_directory, total_records, id_to_organism)

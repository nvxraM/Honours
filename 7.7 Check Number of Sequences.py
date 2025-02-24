import glob
import os


def count_fasta_headers_in_files(pattern):
    """
    Counts the total number of FASTA headers (lines beginning with '>')
    across all files matched by the given glob pattern.
    """
    total_headers = 0

    # Find all files that match the pattern
    file_paths = glob.glob(pattern)

    # Iterate over each file and count headers
    for file_path in file_paths:
        with open(file_path, 'r') as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    total_headers += 1

    return total_headers


if __name__ == "__main__":
    # Adjust the pattern as needed to match your directory structure
    # Example pattern: 'sequences/CDS_Genus/*/CDS_nucleotide_gapped/*_concatenated_gapped_sequences.fasta'
    pattern = 'sequences/CDS_Genus/*/CDS_nucleotide_gapped/*_concatenated_gapped_sequences.fasta'

    total_sequences = count_fasta_headers_in_files(pattern)
    print(f"Total number of sequences: {total_sequences}")

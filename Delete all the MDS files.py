import os
import shutil

# Base directory containing the genus folders
base_dir = "sequences/CDS_Genus"

# Iterate over each genus directory
for genus in os.listdir(base_dir):
    genus_dir = os.path.join(base_dir, genus)

    # Check if it's indeed a directory
    if os.path.isdir(genus_dir):
        # Construct the path to the MDS directory
        mds_dir = os.path.join(genus_dir, "CDS_nucleotide_gapped", "MDS")

        # If the MDS directory exists, remove it and all its contents
        if os.path.isdir(mds_dir):
            shutil.rmtree(mds_dir)
            print(f"Deleted MDS directory and all contents for genus: {genus}")

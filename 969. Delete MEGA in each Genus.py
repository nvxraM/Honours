import os

def delete_mega_folders(base_dir):
    """
    Delete all MEGA folders in each genus directory within the specified base directory.

    Args:
        base_dir (str): The base directory containing genus folders.
    """
    for genus_folder in os.listdir(base_dir):
        genus_dir = os.path.join(base_dir, genus_folder)
        mega_dir = os.path.join(genus_dir, "MEGA_Grouped_by_Species")

        # Check if MEGA directory exists
        if os.path.isdir(mega_dir):
            try:
                # Remove the MEGA directory and its contents
                for root, dirs, files in os.walk(mega_dir, topdown=False):
                    for file in files:
                        os.remove(os.path.join(root, file))
                    for dir in dirs:
                        os.rmdir(os.path.join(root, dir))
                os.rmdir(mega_dir)
                print(f"Deleted MEGA folder: {mega_dir}")
            except Exception as e:
                print(f"Error deleting MEGA folder {mega_dir}: {e}")
        else:
            print(f"MEGA folder does not exist for genus: {genus_folder}")

# Example usage
if __name__ == "__main__":
    base_directory = r"sequences\CDS"
    delete_mega_folders(base_directory)

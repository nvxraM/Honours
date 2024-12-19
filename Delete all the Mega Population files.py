import os
import shutil

base_dir = "sequences/Interactive_group_editor"

for root, dirs, files in os.walk(base_dir, topdown=False):
    for d in dirs:
        if d == "MEGA_Populations":
            dir_path = os.path.join(root, d)
            print(f"Deleting: {dir_path}")
            shutil.rmtree(dir_path)

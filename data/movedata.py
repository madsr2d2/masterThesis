"""
One-off utility that built data/data from the delivered archive.

data/data holds a flat copy of every instrument export (.txt) and recipe
workbook the pipeline reads; data/Mads holds the archive as delivered, in
the experimenter's own folder structure. This is the script that flattened
the first into the second, keeping the mapping reproducible rather than a
manual step nobody recorded.

Not part of the pipeline and not imported by anything. Kept because it
documents where data/data came from. See DATA_VERIFICATION.md for the
archive-to-dataset mapping it produced.
"""
import os
import shutil


def copy_txt_files(source_folder, destination_folder):
    """
    Recursively copies all .txt files from the source folder and its subfolders
    to the destination folder.

    Args:
        source_folder (str): Path to the source folder.
        destination_folder (str): Path to the destination folder.

    Returns:
        None
    """
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    for root, _, files in os.walk(source_folder):
        for file in files:
            if file.endswith(".txt"):
                source_path = os.path.join(root, file)
                destination_path = os.path.join(destination_folder, file)

                # Ensure unique filenames in the destination folder
                counter = 1
                while os.path.exists(destination_path):
                    base, ext = os.path.splitext(file)
                    destination_path = os.path.join(
                        destination_folder, f"{base}_{counter}{ext}"
                    )
                    counter += 1

                shutil.copy2(source_path, destination_path)
                print(f"Copied: {source_path} to {destination_path}")


# Example usage
source_folder = "Mads"
destination_folder = "data"
copy_txt_files(source_folder, destination_folder)

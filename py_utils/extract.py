import shutil
import os
import zipfile
import re
import sys

def extract_zip(zip_path):
    """
    Extract an existing ZIP file into a directory with the same name
    The original ZIP file is kept.
    """

    if not zip_path.endswith(".zip"):
        raise ValueError(f"Expected a .zip file: {zip_path}")

    if not os.path.exists(zip_path):
        raise FileNotFoundError(f"ZIP file not found: {zip_path}")

    extract_path = zip_path[:-4]

    if os.path.exists(extract_path):
        print(f"Already extracted: {extract_path}")
        return extract_path

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(path=extract_path)

    return extract_path


def flatten_product_directory(product_dir, out_dir):
    """
    Move only NetCDF files numbered 0021-0041 from an extracted
    EUMETSAT product directory into out_dir and remove the extracted
    directory afterwards.
    """

    moved = []

    for root, _, files in os.walk(product_dir):
        for filename in files:
            if not filename.endswith(".nc"):
                continue

            # Extract the trailing four-digit file number
            m = re.search(r"_(\d{4})\.nc$", filename)
            if m is None:
                continue

            number = int(m.group(1))

            # Keep only files 21-41 inclusive
            # northern hemisphere only
            if number < 21:
                continue

            src = os.path.join(root, filename)
            dst = os.path.join(out_dir, filename)

            shutil.move(src, dst)
            moved.append(dst)

    shutil.rmtree(product_dir)

    return moved


# Command-line argument
data_dir = os.path.abspath(sys.argv[1])

if not os.path.isdir(data_dir):
    raise FileNotFoundError(f"Demo data directory not found: {data_dir}")

# Find ZIP files
zip_files = [
    f for f in os.listdir(data_dir)
    if f.endswith(".zip")
]

if len(zip_files) == 0:
    raise FileNotFoundError(f"No ZIP file found in {data_dir}")

if len(zip_files) > 1:
    raise RuntimeError(
        f"Expected exactly one ZIP file in {data_dir}, "
        f"but found {len(zip_files)}"
    )

zip_path = os.path.join(data_dir, zip_files[0])

# Extract the local ZIP
product_dir = extract_zip(zip_path)

# Flatten the extracted product into the same directory
flatten_product_directory(product_dir, data_dir)
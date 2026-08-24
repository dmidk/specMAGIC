from pathlib import Path
from urllib.request import urlretrieve
import argparse
import sys

BASE_URL = "http://tds.webservice-energy.org/thredds/fileServer/ground-albedo"

FILES = [
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band1.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band2.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band3.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band4.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band5.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band6.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_Band7.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_nir.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_shortwave.nc",
    "Climatology_monthly_BRDF_parameters_MODIS_MCD43C1_vis.nc",
]


def progress(block_number, block_size, total_size):
    downloaded = block_number * block_size

    if total_size > 0:
        percent = min(100.0, downloaded * 100.0 / total_size)
        downloaded_mb = downloaded / 1024 / 1024
        total_mb = total_size / 1024 / 1024
        sys.stdout.write(
            f"\r    {percent:6.2f}%  {downloaded_mb:8.1f} MB / {total_mb:8.1f} MB"
        )
    else:
        downloaded_mb = downloaded / 1024 / 1024
        sys.stdout.write(f"\r    downloaded {downloaded_mb:8.1f} MB")

    sys.stdout.flush()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output_dir",
        help="Folder where the MODIS BRDF NetCDF files should be saved",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Saving MODIS BRDF files to: {output_dir}", flush=True)

    for filename in FILES:
        destination = output_dir / filename

        if destination.exists():
            print(f"Skipping existing file: {destination}", flush=True)
            continue

        url = f"{BASE_URL}/{filename}"
        print(f"Downloading {url}", flush=True)

        try:
            urlretrieve(url, destination, reporthook=progress)
            print(f"\nSaved to {destination}", flush=True)
        except Exception as error:
            print(f"\nERROR while downloading {filename}: {error}", flush=True)
            raise


if __name__ == "__main__":
    main()
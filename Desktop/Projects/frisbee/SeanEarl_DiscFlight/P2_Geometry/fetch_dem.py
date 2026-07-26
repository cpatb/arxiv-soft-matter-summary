"""
Fetch a USGS 3DEP 1-m DEM tile for the Sean Earl Field bbox.

USGS National Map TNMAccess REST:
  https://tnmaccess.nationalmap.gov/api/v1/products?datasets=Digital%20Elevation%20Model%20(DEM)%201%20meter
The Cook County 1-m DEM is part of the Illinois Height Modernization Program.

Run locally; the sandbox cannot reach USGS reliably.
"""

from __future__ import annotations
import argparse


def fetch_dem(bbox, out_path: str):
    minlon, minlat, maxlon, maxlat = bbox
    print(f"[stub] Would query TNM Access for DEM covering {bbox}")
    print(f"[stub] and download the 1-m DEM tile to {out_path}")
    print()
    print("Steps to do manually if API is being annoying:")
    print("  1. Go to https://apps.nationalmap.gov/downloader/")
    print("  2. Zoom to Loyola Lake Shore Campus")
    print("  3. Filter: 'Elevation Products (3DEP)' -> '1 meter DEM'")
    print("  4. Select the tile, download the .tif, place it at:", out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", nargs=4, type=float,
                    metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"),
                    required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    fetch_dem(tuple(args.bbox), args.out)


if __name__ == "__main__":
    main()

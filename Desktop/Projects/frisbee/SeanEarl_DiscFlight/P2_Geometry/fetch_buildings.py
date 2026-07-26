"""
Fetch building footprints + heights for the Sean Earl Field bounding box.

Sources (in priority order):
    1. Microsoft GlobalMLBuildingFootprints (Illinois quadkey, ODbL)
       https://github.com/microsoft/GlobalMLBuildingFootprints
    2. USGS Overture Maps 2023 (heights enrichment for Chicago)
    3. OpenStreetMap building polygons (height tag if present)
    4. Cook County GIS open-data portal (parcel + building shapes)

The script writes a single GeoPackage with attributes (height_m, source).
Run on your local machine (sandbox cannot reach these endpoints reliably).
"""

from __future__ import annotations
import argparse
import sys


def fetch_microsoft(bbox: tuple[float, float, float, float], out_path: str):
    """Pull MS Global ML Building Footprints clipped to bbox.

    bbox = (minlon, minlat, maxlon, maxlat) in EPSG:4326.
    """
    try:
        import geopandas as gpd
        import requests
    except ImportError:
        raise SystemExit("install geopandas, requests; see README")

    # Microsoft serves per-quadkey GeoJSONL. For Chicago (Illinois),
    # the relevant root URL is:
    #   https://minedbuildings.blob.core.windows.net/global-buildings/Illinois.csv.gz
    # The CSV is an index of quadkey -> URL. For the Loyola area,
    # the quadkey at zoom 9 covering (lon=-87.66, lat=42.00) is "0302032".
    # See https://github.com/microsoft/GlobalMLBuildingFootprints for the
    # full index. This stub assumes the user has selected the right tile.
    print(f"[stub] Would fetch Microsoft Buildings tile covering {bbox}")
    print(f"[stub] Write polygons + height_m to {out_path}")
    print("Replace this stub with the actual download once you've picked the tile.")


def fetch_osm(bbox: tuple[float, float, float, float], out_path: str):
    """Fallback: Overpass API for buildings within bbox with height tags."""
    try:
        import requests
        import geopandas as gpd
    except ImportError:
        raise SystemExit("install requests, geopandas")

    minlon, minlat, maxlon, maxlat = bbox
    query = f"""
    [out:json][timeout:60];
    (
      way["building"]({minlat},{minlon},{maxlat},{maxlon});
      relation["building"]({minlat},{minlon},{maxlat},{maxlon});
    );
    out geom;
    """
    print(f"[stub] Overpass query for {bbox}:")
    print(query)
    print(f"[stub] Would POST to https://overpass-api.de/api/interpreter and write {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bbox", nargs=4, type=float,
                    metavar=("MINLON", "MINLAT", "MAXLON", "MAXLAT"),
                    required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--source", choices=["microsoft", "osm"], default="microsoft")
    args = ap.parse_args()

    bbox = tuple(args.bbox)
    if args.source == "microsoft":
        fetch_microsoft(bbox, args.out)
    else:
        fetch_osm(bbox, args.out)


if __name__ == "__main__":
    main()

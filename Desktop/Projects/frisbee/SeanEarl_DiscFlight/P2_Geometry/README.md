# P2 — Geometry Pipeline for Sean Earl Field

Builds the STL surface model + DEM that P3's OpenFOAM `snappyHexMesh` will
mesh. Three data layers:

1. **Ground**: USGS 3DEP LiDAR DEM → STL ground patch
2. **Buildings**: Microsoft GlobalMLBuildingFootprints + USGS 2023 heights → extruded STL
3. **Tree canopy**: hand-traced polygons over the lakefront strip → porous-zone STL

## Site bounding box (Sean Earl Field)

Centered at approximately:

```
lat  = 42.0006   N
lon  = -87.6580  W   (Loyola Lake Shore Campus, Norville Center area)
```

Recommended domain: 500 m × 500 m horizontal, 0–80 m vertical. The
field center is the origin of the simulation coordinate system; +x East,
+y North, +z up.

Verify the lat/lon before running — Loyola's published campus map
PDF (`luc.edu/media/lucedu/pdfs-campusmaps/2025-map-lsc-8.5x11.pdf`)
gives the field's location relative to Mertz Hall and the Norville Center.

## Files

- `fetch_buildings.py` — pulls Microsoft GlobalMLBuildingFootprints (Illinois) + Cook County GIS polygon + height data for the bounding box.
- `fetch_dem.py` — pulls USGS 3DEP 1 m DEM tile from The National Map.
- `buildings_to_stl.py` — extrudes each polygon to its height, writes STL.
- `dem_to_stl.py` — meshes DEM raster to a triangulated STL ground patch.
- `tree_canopy.py` — defines lakefront tree polygons; writes a porous-zone-region STL.
- `compose_domain.py` — concatenates the three STLs into the final `seanEarl.stl` for `snappyHexMesh`.

## Workflow

```bash
# requires: requests, rasterio, shapely, trimesh, geopandas, pyproj
pip install requests rasterio shapely trimesh geopandas pyproj

python fetch_buildings.py \
    --bbox -87.661 41.998 -87.655 42.003 \
    --out raw_buildings.gpkg

python fetch_dem.py \
    --bbox -87.661 41.998 -87.655 42.003 \
    --out raw_dem.tif

python buildings_to_stl.py raw_buildings.gpkg --out buildings.stl
python dem_to_stl.py raw_dem.tif --out ground.stl --simplify 1.0
python tree_canopy.py --out trees.stl
python compose_domain.py --buildings buildings.stl --ground ground.stl \
                          --trees trees.stl --out seanEarl.stl
```

The output `seanEarl.stl` is fed to P3's `snappyHexMeshDict` as the surface
to refine around.

## Key buildings to include (Lake Shore Campus, east of field)

The following structures are within the snappyHexMesh refinement zone and
shape the local wind:

- **Mertz Residence Hall** (immediately west of field, ~22 m tall) — defines the eastern wake under westerly winds and the lee zone under easterlies.
- **Norville Center for Intercollegiate Athletics** (immediately south, ~12 m).
- **Damen Student Center** (~150 m W, ~25 m).
- **Cudahy Library** (~200 m WNW, ~30 m).
- **Information Commons** (~250 m W, ~20 m).
- Plus 2–3 nearest residential blocks west of N. Sheridan Rd.

## What's stubbed and what needs doing

The Python files in this directory are **scaffolds** with the structure,
argument parsing, and library imports laid out — you need to:

1. **Run them on your own machine.** The sandbox cannot reach USGS or Microsoft endpoints reliably; geospatial data download is intentionally not run here.
2. **Confirm the bounding box** against the Loyola map PDF before downloading data.
3. **Spot-check 3–5 building heights** with a handheld laser rangefinder during a campus visit, especially Mertz and Norville which dominate the field's near wake.
4. **Hand-trace the tree-canopy polygons** in QGIS over a recent aerial image; the script provides a template.

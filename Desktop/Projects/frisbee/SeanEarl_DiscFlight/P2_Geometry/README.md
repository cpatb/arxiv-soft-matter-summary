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

## Status (2026-07-26) — BUILT from real public data

```bash
pip install trimesh shapely mapbox_earcut
python fetch_buildings.py --half-size 350 --out site_geometry.json
python buildings_to_stl.py site_geometry.json --out buildings.stl --half-size 320
python compose_domain.py --buildings buildings.stl --out seanEarl.stl
python verify_site.py site_geometry.json          # -> site_verification.png
```

Output: `seanEarl.stl`, 210 buildings, 10 292 triangles, **210/210 watertight**,
zero inverted normals, bounds x[−337, 236] y[−333, 334] z[0, 83.2] m.

Source is **OpenStreetMap via Overpass**, not Microsoft GlobalMLBuildingFootprints:
OSM has the named Loyola buildings (needed for the campus-map cross-check) and
carries height/levels tags, whereas the MS footprints are geometry-only.
171/303 fetched buildings (56%) have an OSM `building:levels` tag; the rest use
documented type defaults. Every building records its `height_source`.

### Site verification findings — four corrections to the manuscript

1. **Origin was wrong.** The README's assumed centre (42.0006, −87.6580) is
   ~75 m south of the true field centroid. The pipeline now resolves the origin
   from OSM way 210532242 at run time: **42.001267, −87.658214**.
2. **Mertz Hall is 19 storeys ≈ 61 m, not 22 m.** Confirmed by OSM
   `building:levels=19` and by Loyola's own housing page. The manuscript's
   `H ≈ 22 m` under-states the tallest near-field building by ~3×, which
   invalidates its COST 732 domain sizing (see below).
3. **The field's long axis runs E–W, not N–S.** The OSM polygon is 78 m (E–W)
   × 25 m (N–S). Note this is also smaller than a full Ultimate pitch
   (100 × 37 m), so the OSM outline is likely a partial/rough trace — flagged
   for the author.
4. **Cudahy Library and the Information Commons are ESE of the field, not
   west.** The manuscript places them "~200 m WNW" and "~250 m W"; OSM puts them
   at 129 m ESE and 185 m ESE — i.e. on the *lakefront* side.

### The finding that matters most for the physics

**Campus Towers**, a 14-storey (44.8 m) apartment block, occupies
x ∈ [+55.8, +79.3], y ∈ [−15.9, +21.5] — its footprint begins **9 m east of the
field's eastern edge** and spans the full width of the field. For easterly
inflow the field therefore sits in the *near wake* of a 45 m bluff body, not in
a clean marine internal boundary layer. This directly contradicts §3.3 of the
draft ("easterly inflow ... passes over ~80 km of unbroken Lake Michigan fetch
... *this is the dominant local physical feature*"). Reported in §3.7/§3.8.

### Domain sizing consequence

With the true tallest building in the modelled box (H = 83.2 m, Caroline Hedger
Apartments) the COST 732 / Tominaga (2008) guidelines require ≥5H = 416 m
upstream/lateral/above and ≥15H = 1248 m downstream. The manuscript's proposed
500 × 500 × 80 m domain puts the top only **1.0 H** above the tallest building
and gives **3.5% blockage** from that building alone, against a <3% guideline.
P3 records which criteria are met and which are relaxed for compute.

## Scoped down (stated in the paper)

- **Ground is flat**, not a DEM triangulation. USGS 3DEP point queries give
  179.1–181.4 m MSL across the field and near-field buildings (±1 m), against
  building heights of 12–83 m. The ground is the blockMesh bottom patch with a
  rough-wall function. The one real feature, a ~4 m step down to the lake east
  of x ≈ +150 m (3DEP reads 175.2–175.6 m there), is represented as the inlet
  roughness transition rather than as topography.
- **No tree-canopy porous zone.** OSM maps 453 individual trees but no
  contiguous canopy polygon east of the field; the two mapped woodland polygons
  are ~200 m SW and the lakefront parks are NE. Rather than hand-trace an
  unvalidatable canopy, it is omitted — defensible because the dominant
  obstruction between field and lake is Campus Towers, not vegetation.

## Original scaffold notes (superseded)

## What's stubbed and what needs doing

The Python files in this directory are **scaffolds** with the structure,
argument parsing, and library imports laid out — you need to:

1. **Run them on your own machine.** The sandbox cannot reach USGS or Microsoft endpoints reliably; geospatial data download is intentionally not run here.
2. **Confirm the bounding box** against the Loyola map PDF before downloading data.
3. **Spot-check 3–5 building heights** with a handheld laser rangefinder during a campus visit, especially Mertz and Norville which dominate the field's near wake.
4. **Hand-trace the tree-canopy polygons** in QGIS over a recent aerial image; the script provides a template.

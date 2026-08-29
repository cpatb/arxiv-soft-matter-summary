"""
data_chicago.py
================
Curated geographic dataset for the Chicago 'L' Physarum network-growth study.

All coordinates are real-world (longitude, latitude) in decimal degrees, drawn
from CTA station data, City of Chicago community-area centroids, and published
transit-desert analyses (see references in the accompanying LaTeX document).

The module exposes:
  * EXISTING_NODES   : anchor nodes of the current 8-line 'L' system (terminals,
                       major junctions, the Loop core).
  * EXISTING_LINES   : polylines (ordered node keys) describing existing corridors,
                       used to grant a build-cost discount for right-of-way reuse.
  * DEMAND_NODES     : candidate destinations the *new* network should serve, each
                       tagged with population, activity-center, and underserved
                       (transit-desert) sub-scores.
  * LAND_POLYGON     : (lon,lat) vertices of the buildable land area; everything
                       east of the Lake Michigan shoreline is excluded (hard no-build).
  * project()        : equirectangular projection (lon,lat) -> (x,y) kilometres.

The demand weight used to drive the slime-mold flows is a convex combination
    w = a*pop + b*activity + c*underserved
with (a,b,c) configurable in run_experiments.py (default 0.4/0.3/0.3), realising
the user's requested "weighted mix" of population density, activity centres, and
underserved regions.
"""
import numpy as np

# Projection origin (central Chicago)
LON0, LAT0 = -87.68, 41.85
_KM_PER_DEG_LAT = 110.574
_KM_PER_DEG_LON = 111.320 * np.cos(np.radians(LAT0))


def project(lon, lat):
    """Equirectangular projection to local kilometres about (LON0, LAT0)."""
    x = (np.asarray(lon) - LON0) * _KM_PER_DEG_LON
    y = (np.asarray(lat) - LAT0) * _KM_PER_DEG_LAT
    return x, y


# ---------------------------------------------------------------------------
# Existing 'L' system anchor nodes  (key -> (lon, lat))
# Terminals, the downtown Loop core, and a few structural junctions.
# ---------------------------------------------------------------------------
EXISTING_NODES = {
    "Loop":            (-87.6300, 41.8786),  # downtown core (all lines converge)
    "OHare":           (-87.9040, 41.9803),  # Blue NW terminus (airport)
    "ForestPark":      (-87.8170, 41.8743),  # Blue W terminus
    "Howard":          (-87.6731, 42.0192),  # Red/Purple/Yellow N terminus
    "Linden":          (-87.6906, 42.0731),  # Purple terminus (Wilmette)
    "DempsterSkokie":  (-87.7517, 42.0383),  # Yellow terminus
    "Dan95th":         (-87.6244, 41.7224),  # Red S terminus (95th/Dan Ryan)
    "Midway":          (-87.7378, 41.7868),  # Orange terminus (airport)
    "HarlemLake":      (-87.8030, 41.8869),  # Green W terminus
    "Ashland63":       (-87.6638, 41.7794),  # Green SW branch terminus
    "CottageGrove":    (-87.6058, 41.7800),  # Green E (63rd) terminus
    "Cermak54":        (-87.7565, 41.8587),  # Pink terminus
    "Kimball":         (-87.7113, 41.9679),  # Brown terminus
    "RooseveltJct":    (-87.6270, 41.8674),  # south Loop junction
    "BelmontJct":      (-87.6530, 41.9395),  # north-side junction (Red/Brown/Purple)
    "ClarkLake":       (-87.6309, 41.8856),  # Loop (Blue) hub
}

# Ordered polylines approximating existing corridors (for ROW-reuse discount).
EXISTING_LINES = [
    ["OHare", "ForestPark"],                       # (proxy) NW approach
    ["OHare", "BelmontJct", "ClarkLake", "Loop"],  # Blue-ish trunk to core
    ["Linden", "Howard", "BelmontJct", "Loop"],    # Red/Purple N trunk
    ["DempsterSkokie", "Howard"],                  # Yellow
    ["Loop", "RooseveltJct", "Dan95th"],           # Red S (Dan Ryan)
    ["Loop", "RooseveltJct", "Midway"],            # Orange
    ["Loop", "HarlemLake"],                        # Green W
    ["RooseveltJct", "Ashland63"],                 # Green SW
    ["RooseveltJct", "CottageGrove"],              # Green E
    ["Loop", "Cermak54"],                          # Pink
    ["BelmontJct", "Kimball"],                      # Brown
    ["ForestPark", "Loop"],                        # Blue/Forest Park branch
]

# ---------------------------------------------------------------------------
# Candidate demand nodes for the NEW network.
# Each entry: key -> dict(lonlat, pop, activity, underserved)
#   pop        : relative residential population density (0-1)
#   activity   : relative activity-centre / employment draw (0-1)
#   underserved: relative transit-desert severity (0-1)
# Values are normalised scores assembled from community-area population,
# major trip generators, and published Chicago transit-desert maps.
# ---------------------------------------------------------------------------
DEMAND_NODES = {
    # --- Activity centres / major trip generators ---
    "HydePark_UChicago": dict(lonlat=(-87.5987, 41.7886), pop=0.80, activity=0.95, underserved=0.75),
    "McCormickPlace":    dict(lonlat=(-87.6168, 41.8519), pop=0.30, activity=0.85, underserved=0.55),
    "NavyPier":          dict(lonlat=(-87.6086, 41.8917), pop=0.40, activity=0.80, underserved=0.45),
    "IllMedDistrict":    dict(lonlat=(-87.6730, 41.8690), pop=0.55, activity=0.75, underserved=0.20),
    "MuseumCampus":      dict(lonlat=(-87.6170, 41.8650), pop=0.20, activity=0.65, underserved=0.40),
    # --- Dense residential population centroids ---
    "Lakeview":          dict(lonlat=(-87.6540, 41.9400), pop=0.95, activity=0.55, underserved=0.10),
    "LoganSquare":       dict(lonlat=(-87.7080, 41.9230), pop=0.85, activity=0.45, underserved=0.30),
    "Pilsen":            dict(lonlat=(-87.6570, 41.8570), pop=0.75, activity=0.40, underserved=0.45),
    "Bridgeport":        dict(lonlat=(-87.6500, 41.8380), pop=0.70, activity=0.35, underserved=0.55),
    "AlbanyPark":        dict(lonlat=(-87.7160, 41.9680), pop=0.80, activity=0.30, underserved=0.55),
    "BrightonPark":      dict(lonlat=(-87.6990, 41.8170), pop=0.75, activity=0.25, underserved=0.70),
    # --- Transit deserts (far south / west underserved) ---
    "Englewood":         dict(lonlat=(-87.6440, 41.7800), pop=0.60, activity=0.20, underserved=0.95),
    "WestEnglewood":     dict(lonlat=(-87.6670, 41.7754), pop=0.58, activity=0.15, underserved=0.95),
    "Austin":            dict(lonlat=(-87.7615, 41.8916), pop=0.85, activity=0.25, underserved=0.90),
    "HumboldtPark":      dict(lonlat=(-87.7015, 41.9020), pop=0.78, activity=0.25, underserved=0.80),
    "Chatham":           dict(lonlat=(-87.6090, 41.7400), pop=0.65, activity=0.20, underserved=0.85),
    "SouthChicago":      dict(lonlat=(-87.5540, 41.7400), pop=0.55, activity=0.25, underserved=0.90),
    "SouthDeering":      dict(lonlat=(-87.5560, 41.7100), pop=0.45, activity=0.20, underserved=0.95),
    "Roseland":          dict(lonlat=(-87.6200, 41.6900), pop=0.65, activity=0.20, underserved=0.95),
    "Hegewisch":         dict(lonlat=(-87.5460, 41.6560), pop=0.35, activity=0.20, underserved=0.95),
    "BelmontCragin":     dict(lonlat=(-87.7700, 41.9300), pop=0.82, activity=0.25, underserved=0.85),
    "AshburnMidwaySW":   dict(lonlat=(-87.7100, 41.7470), pop=0.60, activity=0.20, underserved=0.80),
}

# ---------------------------------------------------------------------------
# Buildable land polygon.  Vertices traced clockwise from the NW corner along
# the north edge, DOWN the Lake Michigan shoreline (the binding eastern no-build
# boundary), across the south edge, and back up the western edge.  Any candidate
# point east of the shoreline falls outside this polygon and is rejected.
# ---------------------------------------------------------------------------
LAND_POLYGON = [
    (-87.95, 42.08),   # NW corner (near O'Hare / Evanston latitude)
    (-87.66, 42.08),   # NE corner at the north lakeshore
    # --- Lake Michigan shoreline, north -> south (eastern boundary) ---
    (-87.655, 41.99),
    (-87.640, 41.95),
    (-87.620, 41.91),
    (-87.610, 41.88),  # Navy Pier / Loop lakefront
    (-87.605, 41.85),
    (-87.585, 41.81),
    (-87.560, 41.77),
    (-87.535, 41.73),
    (-87.525, 41.70),
    (-87.525, 41.66),  # SE corner near the Indiana state line
    # --- South edge ---
    (-87.95, 41.64),   # SW corner (far south / west)
    # implicit close back to NW corner up the western edge
]


def demand_weight(node_attr, a=0.4, b=0.3, c=0.3):
    """Convex combination of population, activity, and underserved sub-scores."""
    return a * node_attr["pop"] + b * node_attr["activity"] + c * node_attr["underserved"]


if __name__ == "__main__":
    print(f"Existing anchor nodes : {len(EXISTING_NODES)}")
    print(f"Demand (target) nodes : {len(DEMAND_NODES)}")
    print(f"Land polygon vertices : {len(LAND_POLYGON)}")
    xs, ys = project([v[0] for v in EXISTING_NODES.values()],
                     [v[1] for v in EXISTING_NODES.values()])
    print(f"Existing-node extent  : x[{xs.min():.1f},{xs.max():.1f}] "
          f"y[{ys.min():.1f},{ys.max():.1f}] km")

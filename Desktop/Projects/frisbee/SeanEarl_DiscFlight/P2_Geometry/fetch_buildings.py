"""
Fetch building footprints with heights around Sean Earl Field.

SOURCE CHOICE
-------------
The original scaffold proposed Microsoft GlobalMLBuildingFootprints. We use
OpenStreetMap via the Overpass API instead, because:

  * Loyola's Lake Shore Campus is well mapped in OSM, including NAMED campus
    buildings -- and the named buildings are exactly what is needed to
    sanity-check the model against the published campus map;
  * OSM carries `height` / `building:levels` tags, whereas the Microsoft
    footprints are geometry-only and would need a separate height join;
  * the Microsoft dataset ships as continent-scale quadkey shards, a large
    download for a 500 m box.

Buildings without a usable height tag fall back to a storey count, and failing
that to a type-dependent default. Every building records which rule was used in
`height_source`, so the paper can state how many heights are tagged versus
assumed rather than implying all are surveyed.

SITE ORIGIN
-----------
The local Cartesian origin is the OSM centroid of Sean Earl Field itself
(way 210532242), resolved at run time rather than hardcoded. Note that the
coordinate in the P2 README (42.0006, -87.6580) is ~75 m SOUTH of the true
field centroid; using it would have offset the whole CFD domain relative to the
launch grid.

Output: JSON with polygons in LOCAL METRES (x East, y North, origin at the
field centroid), heights in metres, plus the field polygon and metadata.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import urllib.parse
import urllib.request
from collections import Counter

OVERPASS_URLS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
]

# Sean Earl Field, OSM way id. Defines the local origin.
FIELD_WAY_ID = 210532242

# Storey height used where only `building:levels` is tagged. 3.2 m is typical
# for institutional / residential mid-rise and matches OSM rendering defaults.
STOREY_HEIGHT_M = 3.2

# Fallbacks by building type, used only when neither `height` nor
# `building:levels` is present. Deliberately conservative.
DEFAULT_HEIGHT_M = {
    "university": 12.0,
    "college": 12.0,
    "dormitory": 15.0,
    "residential": 9.0,
    "apartments": 12.0,
    "house": 7.0,
    "garage": 3.0,
    "roof": 3.0,
    "shed": 3.0,
}
FALLBACK_HEIGHT_M = 8.0


def overpass(query: str, timeout: int = 180) -> dict:
    """POST an Overpass QL query, trying mirrors in order."""
    data = urllib.parse.urlencode({"data": query}).encode()
    last = None
    for url in OVERPASS_URLS:
        try:
            req = urllib.request.Request(
                url, data=data,
                headers={"User-Agent": "SeanEarlDiscFlight/1.0 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode())
        except Exception as exc:                       # noqa: BLE001
            print(f"  [{url}] failed: {exc}", file=sys.stderr)
            last = exc
    raise SystemExit(f"all Overpass mirrors failed; last error: {last}")


def parse_height(tags: dict) -> tuple[float, str]:
    """Return (height_m, source_rule)."""
    h = tags.get("height")
    if h:
        try:
            hs = str(h).strip().lower().replace("m", "").strip()
            if hs.endswith("'"):
                return float(hs[:-1]) * 0.3048, "height_ft"
            return float(hs), "height_tag"
        except ValueError:
            pass

    lv = tags.get("building:levels")
    if lv:
        try:
            return float(str(lv).split(";")[0]) * STOREY_HEIGHT_M, "levels_tag"
        except ValueError:
            pass

    btype = tags.get("building", "yes")
    if btype in DEFAULT_HEIGHT_M:
        return DEFAULT_HEIGHT_M[btype], f"default:{btype}"
    return FALLBACK_HEIGHT_M, "default:fallback"


def lonlat_to_local(lon, lat, lon0, lat0):
    """Equirectangular projection about (lon0, lat0); accurate over ~1 km."""
    R = 6378137.0
    x = math.radians(lon - lon0) * R * math.cos(math.radians(lat0))
    y = math.radians(lat - lat0) * R
    return x, y


def _ring_from_element(e):
    """Local-frame ring for a way, or the outer ring of a multipolygon."""
    geom = e.get("geometry")
    if not geom or len(geom) < 4:
        outer = [m for m in e.get("members", [])
                 if m.get("role") == "outer" and m.get("geometry")]
        if not outer:
            return None
        geom = outer[0]["geometry"]
        if len(geom) < 4:
            return None
    return geom


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--half-size", type=float, default=350.0,
                    help="half-width of the fetch box in metres")
    ap.add_argument("--out", default="site_geometry.json")
    args = ap.parse_args()

    # --- 1. Resolve the field centroid to define the local origin ----------
    print("resolving Sean Earl Field centroid from OSM ...")
    fd = overpass(f"[out:json][timeout:60];\nway({FIELD_WAY_ID});\nout geom tags;")
    field_el = next((e for e in fd["elements"] if e["type"] == "way"), None)
    if field_el is None:
        raise SystemExit(f"OSM way {FIELD_WAY_ID} (Sean Earl Field) not found")

    fgeom = field_el["geometry"]
    lat0 = sum(p["lat"] for p in fgeom) / len(fgeom)
    lon0 = sum(p["lon"] for p in fgeom) / len(fgeom)
    print(f"  origin = ({lat0:.6f}, {lon0:.6f})  "
          f"[{field_el.get('tags', {}).get('name', '?')}]")

    # --- 2. Fetch every building in the box --------------------------------
    dlat = args.half_size / 111320.0
    dlon = args.half_size / (111320.0 * math.cos(math.radians(lat0)))
    bbox = (lat0 - dlat, lon0 - dlon, lat0 + dlat, lon0 + dlon)
    print(f"fetching buildings in a {2*args.half_size:.0f} m box ...")

    bd = overpass(f"""[out:json][timeout:120];
(
  way["building"]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
  relation["building"]({bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]});
);
out geom tags;""")

    buildings = []
    for e in bd["elements"]:
        geom = _ring_from_element(e)
        if geom is None:
            continue
        tags = e.get("tags", {})
        h, src = parse_height(tags)
        ring = [lonlat_to_local(p["lon"], p["lat"], lon0, lat0) for p in geom]
        if len(ring) > 1 and abs(ring[0][0] - ring[-1][0]) < 1e-9 \
                and abs(ring[0][1] - ring[-1][1]) < 1e-9:
            ring = ring[:-1]
        if len(ring) < 3:
            continue
        buildings.append(dict(
            osm_id=e.get("id"), osm_type=e.get("type"),
            name=tags.get("name"),
            building=tags.get("building", "yes"),
            height_m=h, height_source=src,
            ring=ring,
        ))

    # --- 3. The field polygon itself, for figures and the launch grid ------
    field_ring = [lonlat_to_local(p["lon"], p["lat"], lon0, lat0)
                  for p in fgeom]
    if len(field_ring) > 1 and abs(field_ring[0][0] - field_ring[-1][0]) < 1e-9:
        field_ring = field_ring[:-1]

    out = dict(
        origin_lat=lat0, origin_lon=lon0,
        half_size_m=args.half_size,
        source="OpenStreetMap via Overpass API",
        storey_height_m=STOREY_HEIGHT_M,
        field=dict(name=field_el.get("tags", {}).get("name"),
                   osm_id=FIELD_WAY_ID, ring=field_ring),
        buildings=buildings,
    )
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)

    # --- 4. Report ---------------------------------------------------------
    print(f"\nwrote {args.out}: {len(buildings)} buildings")
    src_counts = Counter(b["height_source"] for b in buildings)
    print("  height provenance:")
    for k, v in src_counts.most_common():
        print(f"    {k:22s} {v:4d}")
    measured = sum(v for k, v in src_counts.items()
                   if not k.startswith("default"))
    print(f"  {measured}/{len(buildings)} heights from OSM tags "
          f"({100*measured/max(len(buildings),1):.0f}%)")

    named = [b for b in buildings if b["name"]]
    named.sort(key=lambda b: -b["height_m"])
    print(f"\n  tallest named buildings ({len(named)} named):")
    for b in named[:14]:
        cx = sum(p[0] for p in b["ring"]) / len(b["ring"])
        cy = sum(p[1] for p in b["ring"]) / len(b["ring"])
        print(f"    {b['name'][:38]:40s} h={b['height_m']:5.1f} m "
              f"({b['height_source']:16s}) at ({cx:+7.1f}, {cy:+7.1f}) m")

    xs = [p[0] for b in buildings for p in b["ring"]]
    ys = [p[1] for b in buildings for p in b["ring"]]
    if xs:
        print(f"\n  extent: x [{min(xs):.0f}, {max(xs):.0f}] m, "
              f"y [{min(ys):.0f}, {max(ys):.0f}] m")
    fx = [p[0] for p in field_ring]
    fy = [p[1] for p in field_ring]
    print(f"  field:  x [{min(fx):.0f}, {max(fx):.0f}] m, "
          f"y [{min(fy):.0f}, {max(fy):.0f}] m "
          f"({max(fx)-min(fx):.0f} x {max(fy)-min(fy):.0f} m)")


if __name__ == "__main__":
    main()

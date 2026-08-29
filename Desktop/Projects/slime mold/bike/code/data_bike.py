"""
data_bike.py -- Chicago bike-network instance for the Physarum gap analysis.
Real inputs encoded here:
  * Facility mileage by type and the named protected corridors come from the City
    of Chicago "Bike Routes" dataset (hvv9-38ut), retrieved via the Socrata API.
  * PROTECTED_CORRIDORS are the highest-mileage protected-lane streets from that
    dataset, encoded as approximate polylines along the real streets.
  * Demand nodes blend residential population, activity/employment draw, and an
    equity/crash-severity score informed by documented West/South-side crash
    hotspots (Austin, Garfield Park, North Lawndale, Greater Grand Crossing).
Projection and demand_weight match data_chicago.py so physarum.CandidateGraph can
be reused by monkeypatching physarum.dc = data_bike.
"""
import numpy as np

LON0, LAT0 = -87.68, 41.85
_KM_PER_DEG_LAT = 110.574
_KM_PER_DEG_LON = 111.320 * np.cos(np.radians(LAT0))


def project(lon, lat):
    x = (np.asarray(lon) - LON0) * _KM_PER_DEG_LON
    y = (np.asarray(lat) - LAT0) * _KM_PER_DEG_LAT
    return x, y


# Real facility mileage by type (City of Chicago Bike Routes hvv9-38ut, Dec 2025)
FACILITY_MILES = {
    "Protected Bike Lane": 68.74,
    "Buffered Bike Lane": 106.46,
    "Bike Lane": 138.49,
    "Neighborhood Greenway": 85.23,
    "Marked Shared Lane": 46.98,
}

# Highest-mileage real PROTECTED corridors (street, miles) and encoded polylines.
# Miles are from the dataset; polylines trace the real streets (approximate).
PROTECTED_CORRIDORS = {
    "Milwaukee Ave":  ([(-87.6430, 41.8930), (-87.6670, 41.9090), (-87.7010, 41.9300), (-87.7170, 41.9470)], 5.09),
    "Lake St":        ([(-87.6320, 41.8855), (-87.7050, 41.8855), (-87.7770, 41.8855)], 4.47),
    "Halsted St":     ([(-87.6470, 41.9500), (-87.6470, 41.8600), (-87.6470, 41.7800)], 4.12),
    "Indiana Ave":    ([(-87.6220, 41.8500), (-87.6220, 41.7900), (-87.6220, 41.7400)], 3.01),
    "Kedzie Ave":     ([(-87.7050, 41.9500), (-87.7050, 41.8700), (-87.7050, 41.8000)], 2.87),
    "Clark St":       ([(-87.6320, 41.8800), (-87.6520, 41.9300), (-87.6680, 41.9700), (-87.6720, 42.0050)], 2.79),
    "Harrison St":    ([(-87.6300, 41.8740), (-87.6800, 41.8740), (-87.7200, 41.8740)], 2.28),
    "Jackson Blvd":   ([(-87.6300, 41.8780), (-87.6900, 41.8770), (-87.7300, 41.8770)], 2.21),
    "Belmont Ave":    ([(-87.6400, 41.9395), (-87.6800, 41.9395), (-87.7150, 41.9395)], 1.86),
    "Dearborn St":    ([(-87.6295, 41.8740), (-87.6295, 41.8970)], 1.85),
    "Vincennes Ave":  ([(-87.6300, 41.7800), (-87.6470, 41.7500), (-87.6620, 41.7250)], 1.84),
    "Damen Ave":      ([(-87.6770, 41.9300), (-87.6770, 41.8800), (-87.6770, 41.8400)], 1.51),
    "Augusta Blvd":   ([(-87.6600, 41.8990), (-87.7000, 41.8990), (-87.7350, 41.8990)], 1.34),
    "31st St":        ([(-87.6250, 41.8380), (-87.6700, 41.8380), (-87.7050, 41.8380)], 1.19),
}

# Build EXISTING_NODES / EXISTING_LINES (keys) programmatically from the corridors.
EXISTING_NODES = {"Loop": (-87.6300, 41.8786)}  # 'Loop' is the solver ground node
EXISTING_LINES = []
for _name, (_pts, _mi) in PROTECTED_CORRIDORS.items():
    keys = []
    for _k, (_lo, _la) in enumerate(_pts):
        key = f"{_name}#{_k}"
        EXISTING_NODES[key] = (_lo, _la)
        keys.append(key)
    EXISTING_LINES.append(keys)

# Demand nodes: pop = population density, activity = employment/trip draw,
# underserved = equity + crash-severity score (West/South side emphasis).
DEMAND_NODES = {
    # Activity / employment centres
    "Loop_Job":          dict(lonlat=(-87.6290, 41.8800), pop=0.55, activity=0.95, underserved=0.20),
    "IllMedDistrict":    dict(lonlat=(-87.6730, 41.8690), pop=0.55, activity=0.80, underserved=0.35),
    "UIC":               dict(lonlat=(-87.6500, 41.8710), pop=0.65, activity=0.80, underserved=0.30),
    "Streeterville_Hosp":dict(lonlat=(-87.6210, 41.8950), pop=0.50, activity=0.80, underserved=0.20),
    "UChicago_HydePark": dict(lonlat=(-87.5987, 41.7886), pop=0.75, activity=0.80, underserved=0.55),
    "McCormick":         dict(lonlat=(-87.6168, 41.8519), pop=0.30, activity=0.60, underserved=0.40),
    # Dense bike-friendly residential
    "Lakeview":          dict(lonlat=(-87.6540, 41.9400), pop=0.95, activity=0.55, underserved=0.10),
    "LoganSquare":       dict(lonlat=(-87.7080, 41.9230), pop=0.90, activity=0.45, underserved=0.30),
    "WickerPark":        dict(lonlat=(-87.6770, 41.9090), pop=0.90, activity=0.55, underserved=0.15),
    "Avondale":          dict(lonlat=(-87.7110, 41.9390), pop=0.80, activity=0.35, underserved=0.45),
    "Uptown":            dict(lonlat=(-87.6550, 41.9660), pop=0.85, activity=0.40, underserved=0.30),
    "Edgewater":         dict(lonlat=(-87.6610, 41.9840), pop=0.85, activity=0.35, underserved=0.25),
    "RogersPark":        dict(lonlat=(-87.6670, 42.0100), pop=0.85, activity=0.35, underserved=0.45),
    "AlbanyPark":        dict(lonlat=(-87.7160, 41.9680), pop=0.80, activity=0.30, underserved=0.55),
    "Pilsen":            dict(lonlat=(-87.6570, 41.8570), pop=0.80, activity=0.40, underserved=0.55),
    "Bridgeport":        dict(lonlat=(-87.6500, 41.8380), pop=0.70, activity=0.35, underserved=0.55),
    "Bronzeville":       dict(lonlat=(-87.6190, 41.8120), pop=0.70, activity=0.40, underserved=0.60),
    # Equity / high-crash West and South side
    "Austin":            dict(lonlat=(-87.7615, 41.8916), pop=0.85, activity=0.25, underserved=0.95),
    "EGarfieldPark":     dict(lonlat=(-87.7000, 41.8810), pop=0.70, activity=0.30, underserved=0.95),
    "NorthLawndale":     dict(lonlat=(-87.7180, 41.8580), pop=0.72, activity=0.25, underserved=0.95),
    "LittleVillage":     dict(lonlat=(-87.7140, 41.8430), pop=0.85, activity=0.30, underserved=0.85),
    "HumboldtPark":      dict(lonlat=(-87.7015, 41.9020), pop=0.80, activity=0.25, underserved=0.80),
    "BelmontCragin":     dict(lonlat=(-87.7700, 41.9300), pop=0.82, activity=0.25, underserved=0.85),
    "Englewood":         dict(lonlat=(-87.6440, 41.7800), pop=0.62, activity=0.20, underserved=0.95),
    "GreaterGrandXing":  dict(lonlat=(-87.6150, 41.7600), pop=0.60, activity=0.20, underserved=0.95),
    "AuburnGresham":     dict(lonlat=(-87.6560, 41.7430), pop=0.65, activity=0.20, underserved=0.90),
    "Chatham":           dict(lonlat=(-87.6090, 41.7400), pop=0.65, activity=0.25, underserved=0.85),
    "SouthShore":        dict(lonlat=(-87.5760, 41.7600), pop=0.65, activity=0.25, underserved=0.80),
    "Roseland":          dict(lonlat=(-87.6200, 41.6900), pop=0.65, activity=0.20, underserved=0.90),
}

# Buildable land polygon (same Lake Michigan no-build boundary as the rail model)
LAND_POLYGON = [
    (-87.95, 42.08), (-87.66, 42.08),
    (-87.655, 41.99), (-87.640, 41.95), (-87.620, 41.91), (-87.610, 41.88),
    (-87.605, 41.85), (-87.585, 41.81), (-87.560, 41.77), (-87.535, 41.73),
    (-87.525, 41.70), (-87.525, 41.66), (-87.95, 41.64),
]


def demand_weight(node_attr, a=0.35, b=0.30, c=0.35):
    return a * node_attr["pop"] + b * node_attr["activity"] + c * node_attr["underserved"]


if __name__ == "__main__":
    print("protected corridors:", len(PROTECTED_CORRIDORS),
          "total encoded mi:", round(sum(m for _, m in PROTECTED_CORRIDORS.values()), 1))
    print("demand nodes:", len(DEMAND_NODES))

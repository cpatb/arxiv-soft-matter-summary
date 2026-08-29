# Bike Gap-Analysis Brief - Verification Notes

## Data provenance (all measured numbers traced to source)
- Facility mileage by type: City of Chicago "Bike Routes" dataset hvv9-38ut (Dec 2025),
  pulled via Socrata aggregate query. Values: Protected 68.74, Buffered 106.46, Bike Lane
  138.49, Neighborhood Greenway 85.23, Marked Shared Lane 46.98 mi. Total 445.9 mi.
  Protected share = 15.4% (verified by recomputation).
- Named protected corridors (Milwaukee 5.1, Lake 4.5, Halsted 4.1, Indiana 3.0, Kedzie 2.9,
  Clark 2.8 ...): same dataset, grouped by street.
- Cyclist crashes: "Traffic Crashes - Crashes" dataset 85ca-t3if, where first_crash_type=
  'PEDALCYCLIST'. Count 17,677; fatalities 40; incapacitating injuries 1,575;
  non-incapacitating 8,760 (confirmed in a single combined query).

## Arithmetic checks (all pass)
- Total mapped miles sum = 445.9; protected = 15.4%; low-stress (protected+greenway) = 34.5%.
- Ideal network TL = 228.7 km; aligned 144.2 + gap 84.5 = 228.7; gap share = 36.9% ~ 37%.
- Ideal-network quality: fault tolerance 0.88, mean trip 12.1 km (model output).

## Modeled vs. measured (stated plainly in the brief, Sec. 6)
Measured: facility mileage, named corridors, crash counts. Modeled: the ideal network,
corridor alignment classification (quarter-mile proximity threshold), and the gap ranking
(demand sub-scores are informed estimates, not calibrated OD ridership; protected corridors
encoded as simplified real-street polylines).

## Attribution
Safety-benefit claim cited to Marshall & Ferenchak (Injury Prevention, 2019). The slime-mold
model is the Tero et al. (2010) framework documented in the companion rail technical report.

## Residual limitations
Demand and protection-status detection would sharpen with Divvy OD flows, ACS/LEHD
demographics, and CDOT's full geocoded facility layer (noted in the brief).

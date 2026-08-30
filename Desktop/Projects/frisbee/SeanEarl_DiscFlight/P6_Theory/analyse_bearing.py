"""Throw direction relative to the wind -- and why the obvious version is wrong.

TWO TRAPS, BOTH LIVE AT THIS SITE.

1. THE RELATIVE ANGLE MUST BE MEASURED AGAINST THE LOCAL DISC-ALTITUDE WIND,
   NOT THE SYNOPTIC DIRECTION. At the release point the two differ by

       theta_wind = 90 deg  (easterly)   local flow differs by -173 deg
       theta_wind = 210 deg (SSW)                              +55 deg
       theta_wind = 270 deg (westerly)                        -120 deg

   The easterly case is a near-perfect REVERSAL: the release point sits in the
   Campus Towers wake, where the near-ground return flow runs back toward the
   building. So a throw that is upwind of the synoptic easterly is downwind of
   the air the disc is actually in. Any relative-angle analysis built on the
   met-station direction would have the sign backwards for a third of the
   headings here.

2. FROM A FIXED LAUNCH POINT THE BEARING AXIS IS CONFOUNDED WITH THE SITE.
   The thrower stands at the west end, so a westward throw flies straight into
   the built-up campus: the fraction of throws landing inside a building runs
   from 0% to 62% across the eight bearings. Those cells are measuring
   architecture, not wind. Cells with more than 2% blocked are dropped here and
   the count is reported, rather than a curve being drawn through all of them.

   Fixing this properly needs the sweep repeated from several launch points --
   which is what the production 12-point grid is for. That is not done here.

WITHIN-HEADING, NOT POOLED. Each heading has its own local wind SPEED (0.59 to
2.23 m/s at the release point), so pooling all three mixes an amplitude axis
into an angle axis. The comparison is made inside each heading.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P4_Coupling"))
sys.path.insert(0, str(ROOT / "P5_Synthesis"))

BLOCKED_MAX = 0.02


def local_wind(field, heading, launch):
    from cfd_field import CFDField
    f = CFDField(str(field), heading=int(heading))
    U = f.U(np.asarray(launch, float), 0.0)
    return (float(np.degrees(np.arctan2(U[1], U[0])) % 360.0),
            float(np.hypot(U[0], U[1])))


def main():
    import site_view as sv
    src = HERE / "out" / "sweep_bearing.json"
    d = json.loads(src.read_text())
    rows = d["rows"]
    field = d["meta"]["field"]

    loc = {h: local_wind(field, h, sv.LAUNCH)
           for h in sorted({r["heading"] for r in rows})}
    print("release-point wind direction, local vs synoptic:")
    for h, (ang, spd) in loc.items():
        syn = (90.0 - ((h + 180.0) % 360.0)) % 360.0
        print(f"  theta_wind {h:>3}: local blows toward {ang:6.1f} deg "
              f"(|U| {spd:.2f} m/s), synoptic {syn:6.1f} deg, "
              f"difference {((ang - syn + 180) % 360) - 180:+7.1f} deg")

    blocked = np.array([r["full__blocked"] for r in rows])
    keep = [r for r in rows if r["full__blocked"] <= BLOCKED_MAX]
    print(f"\nblocked fraction across {len(rows)} cells: "
          f"min {blocked.min():.0%}, median {np.median(blocked):.0%}, "
          f"max {blocked.max():.0%}")
    print(f"{len(keep)} cells kept at <= {BLOCKED_MAX:.0%} blocked")

    out = []
    for h, (ang, spd) in loc.items():
        rs = [r for r in keep if r["heading"] == h]
        if not rs:
            continue
        print(f"\n  theta_wind {h}, local |U| = {spd:.2f} m/s "
              f"({len(rs)} unblocked cells)")
        print(f"    {'Uref':>5} {'rel':>6} {'bias':>7} {'sig_wind':>8} "
              f"{'sig_gust':>8} {'total':>7} {'range':>7}")
        for r in sorted(rs, key=lambda r: (r["u_ref"],
                                           abs(((r["bearing"] - ang + 180)
                                                % 360) - 180))):
            rel = ((r["bearing"] - ang + 180) % 360) - 180
            tw = float(np.hypot(np.hypot(r["wind__bias"], r["wind__sigma"]),
                                r["gust__sigma"]))
            print(f"    {r['u_ref']:5.1f} {rel:6.0f} {r['wind__bias']:7.2f} "
                  f"{r['wind__sigma']:8.2f} {r['gust__sigma']:8.2f} "
                  f"{tw:7.2f} {r['full__down_mean']:7.1f}")
            out.append(dict(heading=h, u_ref=r["u_ref"], bearing=r["bearing"],
                            rel_local_deg=rel, local_wind_toward_deg=ang,
                            local_wind_speed=spd, bias=r["wind__bias"],
                            sigma_wind=r["wind__sigma"],
                            sigma_gust=r["gust__sigma"], total_wind=tw,
                            range_m=r["full__down_mean"]))

    # The one ordering claim, made within each heading at the higher speed.
    print("\n  downwind vs upwind, within each heading at Uref = 10 m/s:")
    for h in loc:
        rs = [o for o in out if o["heading"] == h and o["u_ref"] == 10.0]
        if len(rs) < 2:
            continue
        near = min(rs, key=lambda o: abs(o["rel_local_deg"]))
        far = max(rs, key=lambda o: abs(o["rel_local_deg"]))
        print(f"    {h:>3}: |rel| {abs(near['rel_local_deg']):3.0f} deg -> "
              f"{near['total_wind']:5.2f} m   |rel| "
              f"{abs(far['rel_local_deg']):3.0f} deg -> {far['total_wind']:5.2f} m"
              f"   ratio {far['total_wind']/near['total_wind']:.2f}x")

    p = HERE / "out" / "bearing_results.json"
    p.write_text(json.dumps(dict(local_wind={str(k): v for k, v in loc.items()},
                                 blocked_max=BLOCKED_MAX,
                                 n_cells=len(rows), n_kept=len(keep),
                                 rows=out), indent=2))
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()

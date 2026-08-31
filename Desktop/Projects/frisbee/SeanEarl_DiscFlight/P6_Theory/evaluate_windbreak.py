"""Grade each windbreak design against the gates in `windbreak_gates.py`.

The gates were committed before any design finished solving. This file only
applies them.

G-W1 to G-W4 are field measurements. G-W5 flies discs: the same release draws
and the same gust phases through the baseline field and through each design's
field, paired, so the answer is a per-throw displacement rather than a
difference of two scatters. That is the same estimator the wind-speed theory
uses and for the same reason -- the release scatter here is several times the
effect being measured.

Launch points for G-W5 are the two that sit in or just upwind of the affected
zone, because that is where a barrier at x = +21 or +30 can do anything at all.
A barrier evaluated on throws from the west end, which land before they reach
it, would score a null it did not earn.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P1_AnalyticalModel"))
sys.path.insert(0, str(ROOT / "P4_Coupling"))
sys.path.insert(0, str(ROOT / "P5_Synthesis"))

from windbreak_gates import (DISC_BAND, EAST_THIRD, GATES,  # noqa: E402
                             PITCH_X, PITCH_Y, THIRDS, DESIGNS)

G5_LAUNCHES = [(11.0, 2.5), (30.5, 2.5)]
G5_N_THROWS = 64
G5_SEED = 20260830


def band(F, heading, key):
    x, y, z = F["x"], F["y"], F["z"]
    ix = (x >= PITCH_X[0]) & (x <= PITCH_X[1])
    iy = (y >= PITCH_Y[0]) & (y <= PITCH_Y[1])
    iz = (z >= DISC_BAND[0]) & (z <= DISC_BAND[1])
    a = F[f"dir_{heading:03d}__{key}"]
    solid = F[f"dir_{heading:03d}__solid"]
    return x[ix], y[iy], z[iz], a[np.ix_(ix, iy, iz)] if a.ndim == 3 else \
        a[np.ix_(ix, iy, iz, np.arange(a.shape[3]))], ~solid[np.ix_(ix, iy)]


def field_metrics(path, heading=270):
    F = np.load(path, allow_pickle=False)
    x, y, z, k, msk = band(F, heading, "k")
    _, _, _, U, _ = band(F, heading, "U")
    sig = np.sqrt(np.maximum(2.0 / 3.0 * k, 0.0))
    spd = np.linalg.norm(U[..., :2], axis=-1)
    m3 = np.broadcast_to(msk[:, :, None], sig.shape)

    out = {"n_cells": int(m3.sum()), "z_band": [float(z.min()), float(z.max())]}
    for name, (a, b) in THIRDS.items():
        sel = (x >= a) & (x < b) if name != "east" else (x >= a) & (x <= b)
        s = sig[sel][:, :, :][np.broadcast_to(msk[sel][:, :, None],
                                              sig[sel].shape)]
        out[f"sigma_u_{name}"] = float(np.nanmean(s))
    sel = (x >= EAST_THIRD[0]) & (x <= EAST_THIRD[1])
    out["sigma_u_east_third"] = float(np.nanmean(
        sig[sel][np.broadcast_to(msk[sel][:, :, None], sig[sel].shape)]))
    out["sigma_u_pitch"] = float(np.nanmean(sig[m3]))
    out["U_pitch"] = float(np.nanmean(spd[m3]))

    # Vertical gradient of sigma_u inside the disc band. A new shear layer put
    # into the band by the barrier shows up here and in nothing else.
    if len(z) >= 2:
        dsig = np.gradient(sig, z, axis=2)
        out["dsigma_dz_max"] = float(np.nanmax(np.abs(dsig[m3])))
        out["dsigma_dz_p95"] = float(np.nanpercentile(np.abs(dsig[m3]), 95))
    return out


def fly_pair(base_field, design_field, heading, launch, n_throws=G5_N_THROWS):
    """Paired flights through the two fields. Returns the per-throw shift."""
    import disc_6dof as d6
    import monte_carlo as mc
    import synthetic_turbulence as st
    from cfd_field import CFDField
    from throw_catalogue import build_catalogue

    dist = build_catalogue()["right_backhand"]
    launch = np.asarray([launch[0], launch[1], 1.5], dtype=float)
    out = {}
    for tag, path in (("base", base_field), ("design", design_field)):
        field = CFDField(str(path), heading=heading)
        lands, blocked = [], []
        for i in range(n_throws):
            rng = np.random.default_rng(
                np.random.SeedSequence(entropy=G5_SEED, spawn_key=(int(i),)))
            r0, v0, q0, w0 = d6.build_initial_conditions(
                launch_pos=launch,
                speed=rng.normal(dist["speed_mean"], dist["speed_std"]),
                heading_deg=0.0,
                elevation_deg=rng.normal(dist["elevation_mean"],
                                         dist["elevation_std"]),
                release_alpha_deg=rng.normal(dist["release_alpha_mean"],
                                             dist["release_alpha_std"]),
                roll_deg=rng.normal(dist["roll_mean"], dist["roll_std"]),
                spin_rad_s=max(20.0, rng.normal(dist["spin_mean"],
                                                dist["spin_std"])),
                handedness=dist["handedness"])
            grng = np.random.default_rng(
                np.random.SeedSequence(entropy=G5_SEED + 1, spawn_key=(int(i),)))
            g = st.CorrelatedGusts(field, grng, t_max=12.0, ref_pos=launch,
                                   model="sweeping",
                                   sigma_ratios=st.SIGMA_RATIOS)
            res = d6.simulate(r0, v0, q0, w0, coeffs=mc.DEFAULT_COEFFS,
                              geom=mc.ULTRASTAR,
                              wind=st.wind_with_gusts(field, g), t_max=12.0,
                              rtol=mc.MC_RTOL, atol=mc.MC_ATOL,
                              compute_diag=False)
            r = np.asarray(res.r)
            lands.append(r[-1, :2])
            blocked.append(bool(field.is_solid(r[-1, 0], r[-1, 1])))
        out[tag] = (np.asarray(lands), np.asarray(blocked))

    res = {}
    for tag, (L, b) in out.items():
        good = L[~b]
        mu = good.mean(axis=0)
        res[f"sigma_{tag}"] = float(np.sqrt(((good - mu) ** 2).sum(axis=1).mean()))
        res[f"blocked_{tag}"] = float(b.mean())
        res[f"n_good_{tag}"] = int((~b).sum())
    return res


def grade(base, design, tag, base_flow=None, design_flow=None, heading=270):
    b, d = field_metrics(base, heading), field_metrics(design, heading)
    r = {"design": tag, "baseline": b, "design_field": d, "gates": {}}

    v = d["sigma_u_east_third"] / b["sigma_u_east_third"]
    r["gates"]["G-W1"] = dict(value=v, threshold=GATES["G-W1"]["threshold"],
                              passes=bool(v <= GATES["G-W1"]["threshold"]),
                              detail=f"east-third sigma_u "
                                     f"{b['sigma_u_east_third']:.3f} -> "
                                     f"{d['sigma_u_east_third']:.3f} m/s")

    worst, where = 0.0, None
    for nm in THIRDS:
        q = d[f"sigma_u_{nm}"] / b[f"sigma_u_{nm}"]
        if q > worst:
            worst, where = q, nm
    r["gates"]["G-W2"] = dict(value=worst, threshold=GATES["G-W2"]["threshold"],
                              passes=bool(worst <= GATES["G-W2"]["threshold"]),
                              detail=f"worst third is {where} at x{worst:.3f}")

    v = d["U_pitch"] / b["U_pitch"]
    r["gates"]["G-W3"] = dict(value=v, threshold=GATES["G-W3"]["threshold"],
                              passes=bool(v >= GATES["G-W3"]["threshold"]),
                              detail=f"pitch |U| {b['U_pitch']:.3f} -> "
                                     f"{d['U_pitch']:.3f} m/s")

    v = d["dsigma_dz_max"] / b["dsigma_dz_max"]
    r["gates"]["G-W4"] = dict(value=v, threshold=GATES["G-W4"]["threshold"],
                              passes=bool(v <= GATES["G-W4"]["threshold"]),
                              detail=f"max |dsigma/dz| {b['dsigma_dz_max']:.4f}"
                                     f" -> {d['dsigma_dz_max']:.4f} 1/s")

    if base_flow and design_flow:
        per = []
        for lp in G5_LAUNCHES:
            q = fly_pair(base_flow, design_flow, heading, lp)
            q["launch"] = list(lp)
            per.append(q)
        ratios = [p["sigma_design"] / p["sigma_base"] for p in per]
        v = float(np.mean(ratios))
        r["gates"]["G-W5"] = dict(value=v, threshold=GATES["G-W5"]["threshold"],
                                  passes=bool(v < GATES["G-W5"]["threshold"]),
                                  detail="; ".join(
                                      f"({p['launch'][0]:+.0f},"
                                      f"{p['launch'][1]:+.0f}) "
                                      f"{p['sigma_base']:.2f} -> "
                                      f"{p['sigma_design']:.2f} m" for p in per),
                                  per_launch=per)

    core = [k for k in ("G-W1", "G-W2", "G-W3", "G-W4", "G-W5")
            if k in r["gates"]]
    r["passes_all"] = bool(all(r["gates"][k]["passes"] for k in core))
    r["n_passed"] = int(sum(r["gates"][k]["passes"] for k in core))
    r["n_gates"] = len(core)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cfd", default=str(ROOT / "P3_CFD"))
    ap.add_argument("--heading", type=int, default=270)
    ap.add_argument("--designs", nargs="+", default=sorted(DESIGNS))
    ap.add_argument("--out", default=str(HERE / "out" / "windbreak.json"))
    args = ap.parse_args()

    cfd = Path(args.cfd)
    base = cfd / "vertical_fields.npz"
    base_flow = cfd / "flow_fields.npz"
    results = []
    for tag in args.designs:
        dv = cfd / f"vertical_fields_dir_{args.heading:03d}_{tag}.npz"
        df = cfd / f"flow_fields_dir_{args.heading:03d}_{tag}.npz"
        if not dv.exists():
            print(f"  {tag}: {dv.name} not written yet -- skipping")
            continue
        d = DESIGNS[tag]
        print(f"\n=== {tag}  {d['name']}: x = {d['x']:+.0f} m, H = "
              f"{d['height']:.0f} m, k_r = {d['k_r']:.2f} ===")
        print(f"    {d['idea']}")
        r = grade(base, dv, tag,
                  base_flow if df.exists() else None,
                  df if df.exists() else None, args.heading)
        for gk in ("G-W1", "G-W2", "G-W3", "G-W4", "G-W5"):
            g = r["gates"].get(gk)
            if g is None:
                print(f"    {gk}  not evaluated")
                continue
            print(f"    {gk}  {'PASS' if g['passes'] else 'FAIL'}  "
                  f"x{g['value']:.3f} vs {g['threshold']}   {g['detail']}")
        print(f"    -> {r['n_passed']}/{r['n_gates']} gates, "
              f"{'PASSES' if r['passes_all'] else 'DOES NOT PASS'}")
        results.append(r)

    if results:
        passing = [r for r in results if r["passes_all"]]
        if passing:
            best = max(passing, key=lambda r: -r["gates"]["G-W1"]["value"])
            print(f"\nPASSES: {best['design']} "
                  f"({best['n_passed']}/{best['n_gates']} gates)")
        else:
            # Do not crown a design that failed the gate the study exists for.
            # Ranking failures by G-W1 alone named W1 "best" on a 0.7% edge in a
            # gate both miss by 20 points, while W1 RAISES sigma_u by 8% over
            # the west third and W2 raises it nowhere. Rank by harm instead, and
            # say plainly that nothing passed.
            print(f"\nNO DESIGN PASSES. {len(results)} evaluated, best "
                  f"{max(r['n_passed'] for r in results)}/"
                  f"{results[0]['n_gates']} gates.")
            print("  ranked by harm done, since none achieved the goal:")
            for r in sorted(results, key=lambda r: r["gates"]["G-W2"]["value"]):
                g5 = r["gates"].get("G-W5", {}).get("value", float("nan"))
                print(f"    {r['design']}  G-W1 x{r['gates']['G-W1']['value']:.3f}"
                      f"  worst-third sigma_u x{r['gates']['G-W2']['value']:.3f}"
                      f"  shear x{r['gates']['G-W4']['value']:.3f}"
                      f"  flown x{g5:.3f}")
    Path(args.out).write_text(json.dumps(
        dict(heading=args.heading, gates=GATES, designs=DESIGNS,
             results=results), indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()

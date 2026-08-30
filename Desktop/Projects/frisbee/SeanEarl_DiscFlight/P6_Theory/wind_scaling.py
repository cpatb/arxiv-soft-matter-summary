"""Put the CFD field at any reference wind speed, and test whether that is legal.

THE PROBLEM
-----------
Everything this project has computed sits at ONE wind speed: `Uref = 5.0` m/s at
`Zref = 10` m. A theory of how throwable a site is "as a function of wind speed"
needs the field at every speed, and one OpenFOAM solve per speed per direction
is not affordable on 11 cores.

THE ESCAPE, AND WHY IT IS AN ASSUMPTION
---------------------------------------
Bluff-body flow at high Reynolds number is conventionally taken to be
Reynolds-INDEPENDENT: separation is fixed by sharp edges rather than by the
boundary layer, so the whole solution is self-similar in the inflow speed. In
that limit the steady RANS solution obeys

    U(x; s*Uref)     = s   * U(x; Uref)
    k(x; s*Uref)     = s^2 * k(x; Uref)
    omega(x; s*Uref) = s   * omega(x; Uref)
    eps = Cmu*k*omega            ->  s^3 * eps(x; Uref)

Three consequences follow immediately and they are the backbone of the theory:

    sigma_u    = sqrt(2k/3)   ~  s          -- amplitude scales
    T_L        = k/eps        ~  1/s        -- eddies turn over faster
    L = sigma_u * T_L         ~  s^0        -- THE INTEGRAL LENGTH IS INVARIANT

So under this scaling the site's turbulence *structure* -- its length scales,
its turbulence intensity sigma_u/U, the spatial pattern of shelter -- does not
depend on wind speed at all. Only the amplitude does. Every wind-speed
dependence in the throwability of a site therefore comes from the DISC, not
from the site. That is a strong claim and it is worth being sure of.

WHY IT IS NOT AN IDENTITY HERE
------------------------------
Two things in the case setup carry an explicit Reynolds number:

  * molecular viscosity is finite, `nu = 1.5e-5` m^2/s (`constant/transportProperties`);
  * the `buildings` patch uses `nutkWallFunction` and `omegaWallFunction` -- the
    SMOOTH-wall forms, whose blending depends on y+ and therefore on speed. The
    `ground` patch uses the roughness-based `atmNutkWallFunction`, which is
    scale-free; the buildings' is not.

So the scaling is an approximation whose error must be measured, not asserted.
`P3_CFD/runScaleTest.sh` solves the identical gridL1 mesh at Uref = 2.5, 5 and
10 m/s and `scale_invariance_report()` below reports the residual after
dividing out the expected powers. Nothing in this module claims the scaling is
exact; `ScaledField` records the measured residual so any result built on it
carries its own error bar.

WHAT SCALING DOES *NOT* COVER
-----------------------------
Buoyancy. The inlet is a neutral Richards & Hoxey ABL, so the solution has no
stratification and the scaled family is a family of NEUTRAL boundary layers at
different speeds. A real 1 m/s day is usually not neutral. The theory is
therefore about the mechanical part of the problem, and says so.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "P4_Coupling"))

from cfd_field import CFDField                             # noqa: E402

# Every case template in P3_CFD/cases is generated at this reference speed.
BASE_UREF = 5.0
Z_REF = 10.0

# Exponents of the reference speed for each field, from the self-similar limit.
SCALE_POWERS = {"U": 1, "k": 2, "eps": 3}


@dataclass
class ScaledField:
    """A `CFDField` re-expressed at a different reference wind speed.

    Wraps rather than copies: the underlying interpolator is untouched and the
    scale factors are applied on the way out, so the memory cost of a wind-speed
    sweep is one field, not one per speed.

    `residual` is the measured departure from exact self-similarity, in the same
    units as the quantity it applies to (a fraction). It is carried so a caller
    can propagate it; it is NOT applied as a correction, because a single-number
    correction to a field-valued error would be worse than reporting it.
    """

    base: CFDField
    u_ref: float
    residual: dict | None = None

    @property
    def s(self) -> float:
        return self.u_ref / BASE_UREF

    def U(self, r, t=0.0):
        return self.s * self.base.U(r, t)

    def U_and_k(self, r, t=0.0):
        U, k = self.base.U_and_k(r, t)
        return self.s * U, self.s ** 2 * k

    def k(self, r, t=0.0):
        return self.s ** 2 * self.base.k(r, t)

    def eps(self, r, t=0.0):
        return self.s ** 3 * self.base.eps(r, t)

    def is_solid(self, x, y):
        return self.base.is_solid(x, y)

    @property
    def heading(self):
        return self.base.heading

    @property
    def fraction_clamped(self):
        return self.base.fraction_clamped

    def bounds(self):
        return self.base.bounds()


def scaled_field(npz_path, heading, u_ref, residual=None):
    return ScaledField(CFDField(str(npz_path), heading=heading), float(u_ref),
                       residual)


# ---------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------

def _load(npz_path, heading):
    F = np.load(npz_path, allow_pickle=False)
    key = f"dir_{int(heading):03d}"
    return dict(x=F["x"], y=F["y"], z=F["z"], U=F[f"{key}__U"],
                k=F[f"{key}__k"], eps=F[f"{key}__eps"],
                solid=F[f"{key}__solid"])


def scale_invariance_report(cases, heading=90, pitch_x=(-31.0, 47.0),
                            pitch_y=(-10.0, 15.0)):
    """Compare solutions at several Uref after dividing out the expected powers.

    `cases` is {u_ref: npz_path}. The case at `BASE_UREF` is the reference.

    Reported per field: the median and 95th-percentile absolute relative
    departure from the reference, over cells that are fluid in BOTH solutions
    and inside the playing surface. Restricting to the pitch matters -- the
    sampling box includes courtyards and near-wall cells where the smooth-wall
    building function does most of its damage and where no disc ever flies.
    A box-wide figure is reported alongside so the restriction is visible rather
    than silent.
    """
    if BASE_UREF not in cases:
        raise ValueError(f"no reference case at Uref = {BASE_UREF}")
    ref = _load(cases[BASE_UREF], heading)
    x, y = ref["x"], ref["y"]
    in_pitch = np.ix_((x >= pitch_x[0]) & (x <= pitch_x[1]),
                      (y >= pitch_y[0]) & (y <= pitch_y[1]))

    out = {"base_u_ref": BASE_UREF, "heading": int(heading), "cases": {}}
    for u_ref, path in sorted(cases.items()):
        if u_ref == BASE_UREF:
            continue
        cur = _load(path, heading)
        for nm in ("x", "y", "z"):
            if not np.array_equal(cur[nm], ref[nm]):
                raise ValueError(f"{path} is not on the reference {nm} grid")
        s = u_ref / BASE_UREF
        fluid = (~ref["solid"]) & (~cur["solid"])
        entry = {"u_ref": u_ref, "s": s, "fields": {}}
        for nm, power in SCALE_POWERS.items():
            a = ref[nm]
            b = cur[nm] / s ** power           # rescaled to the reference speed
            if a.ndim == 4:                    # U: compare the vector magnitude
                a = np.linalg.norm(a, axis=-1)
                b = np.linalg.norm(b, axis=-1)
            m3 = np.broadcast_to(fluid[:, :, None], a.shape)
            # Relative to the field's own RMS, not to the local value: dividing
            # by a near-zero local speed in a wake makes the metric diverge
            # exactly where the flow is least resolved. That mistake is what
            # produced the withdrawn "grid independence FAIL, 57.1%".
            denom = float(np.sqrt(np.mean(a[m3] ** 2)))
            d = np.abs(b - a) / max(denom, 1e-12)
            dp = d[in_pitch[0], in_pitch[1], :]
            fp = np.broadcast_to(fluid[in_pitch][:, :, None], dp.shape)
            entry["fields"][nm] = dict(
                pitch_median=float(np.median(dp[fp])),
                pitch_p95=float(np.percentile(dp[fp], 95)),
                pitch_max=float(np.max(dp[fp])),
                box_median=float(np.median(d[m3])),
                box_p95=float(np.percentile(d[m3], 95)),
                n_cells=int(fp.sum()))
        out["cases"][str(u_ref)] = entry
    return out


def derived_invariants(cases, heading=90, z_target=2.0,
                       pitch_x=(-31.0, 47.0), pitch_y=(-10.0, 15.0)):
    """The three quantities the theory rests on, measured at each wind speed.

    sigma_u ~ s, T_L ~ 1/s and L ~ s^0 are what make a single CFD solve enough.
    They are DERIVED from U, k and eps, so they can fail even if each raw field
    passes -- L is a ratio of two scaled quantities and a 2% error in each can
    land anywhere between 0 and 4% in the ratio. Measuring them directly is the
    check that matters for the theory, as opposed to for the solver.
    """
    rows = []
    for u_ref, path in sorted(cases.items()):
        F = _load(path, heading)
        x, y, z = F["x"], F["y"], F["z"]
        kz = int(np.argmin(np.abs(z - z_target)))
        ix = (x >= pitch_x[0]) & (x <= pitch_x[1])
        iy = (y >= pitch_y[0]) & (y <= pitch_y[1])
        m = ~F["solid"][np.ix_(ix, iy)]
        k = F["k"][np.ix_(ix, iy)][:, :, kz][m]
        eps = np.maximum(F["eps"][np.ix_(ix, iy)][:, :, kz][m], 1e-12)
        U = np.linalg.norm(F["U"][np.ix_(ix, iy)][:, :, kz, :2], axis=-1)[m]
        sig = np.sqrt(2.0 / 3.0 * np.maximum(k, 0.0))
        T_L = k / eps
        rows.append(dict(u_ref=float(u_ref),
                         U_mean=float(np.mean(U)),
                         sigma_u=float(np.mean(sig)),
                         T_L=float(np.mean(T_L)),
                         L_int=float(np.mean(sig * T_L)),
                         TI=float(np.mean(sig) / max(np.mean(U), 1e-12)),
                         n=int(m.sum())))
    return rows


def _fit_power(u, q):
    """Least-squares exponent p in q ~ u^p, plus the r^2 of the log-log fit."""
    lu, lq = np.log(np.asarray(u, float)), np.log(np.asarray(q, float))
    p, c = np.polyfit(lu, lq, 1)
    resid = lq - (p * lu + c)
    ss = np.sum((lq - lq.mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss if ss > 0 else float("nan")
    return float(p), float(r2)


def invariance_verdict(rows, tol_exponent=0.10):
    """Grade the derived invariants against their predicted exponents.

    The band is on the EXPONENT, not on the value, because the exponent is what
    the theory uses: it is what lets one solve stand in for a whole family. The
    tolerance is pre-registered here at +/- 0.10 in the exponent, which at a
    4x speed range is a factor of 1.15 across the range.
    """
    u = [r["u_ref"] for r in rows]
    predicted = {"sigma_u": 1.0, "T_L": -1.0, "L_int": 0.0, "U_mean": 1.0,
                 "TI": 0.0}
    verdict = {"tol_exponent": tol_exponent, "u_ref": u, "quantities": {}}
    allpass = True
    for q, p_pred in predicted.items():
        vals = [r[q] for r in rows]
        p_fit, r2 = _fit_power(u, vals)
        ok = abs(p_fit - p_pred) <= tol_exponent
        allpass &= ok
        # r^2 is reported but is NOT part of the verdict, and for the two
        # predicted-invariant quantities it is close to meaningless: a
        # quantity that does not vary has almost no variance to explain, so
        # r^2 collapses toward whatever the round-off does. L_int's 0.97 and
        # TI's 0.59 are that, not evidence of a poor fit -- both are flat to
        # four decimals. Only the EXPONENT is graded.
        spread = (max(vals) - min(vals)) / max(abs(np.mean(vals)), 1e-30)
        verdict["quantities"][q] = dict(predicted=p_pred, fitted=round(p_fit, 4),
                                        r2=round(r2, 5),
                                        relative_spread=round(float(spread), 6),
                                        r2_meaningful=bool(spread > 0.05),
                                        passes=bool(ok))
    verdict["passes"] = bool(allpass)
    verdict["n_speeds"] = len(rows)
    if len(rows) < 3:
        verdict["warning"] = ("only %d speeds: a two-point log-log fit is exact "
                              "by construction and its r^2 = 1 means nothing"
                              % len(rows))
    return verdict


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cfd", default=str(ROOT / "P3_CFD"))
    ap.add_argument("--heading", type=int, default=90)
    ap.add_argument("--out", default=str(HERE / "out" / "scale_invariance.json"))
    args = ap.parse_args()

    cfd = Path(args.cfd)
    cases = {}
    for u_ref, name in ((2.5, "flow_fields_dir_090_gridL1_u2p5.npz"),
                        (5.0, "flow_fields_dir_090_gridL1.npz"),
                        (10.0, "flow_fields_dir_090_gridL1_u10.npz")):
        p = cfd / name
        if p.exists():
            cases[u_ref] = str(p)
        else:
            print(f"  missing {p.name} -- skipping Uref = {u_ref}")
    if len(cases) < 2:
        raise SystemExit("need at least two solved speeds; run "
                         "P3_CFD/runScaleTest.sh first")

    print(f"speeds solved: {sorted(cases)}  (heading {args.heading})")
    rep = scale_invariance_report(cases, args.heading)
    for u, e in rep["cases"].items():
        print(f"\n  Uref = {u} m/s   (s = {e['s']:g})")
        for nm, f in e["fields"].items():
            print(f"    {nm:>4}: pitch median {100*f['pitch_median']:5.2f}%  "
                  f"p95 {100*f['pitch_p95']:6.2f}%  max {100*f['pitch_max']:6.2f}%"
                  f"   | box median {100*f['box_median']:5.2f}%  "
                  f"p95 {100*f['box_p95']:6.2f}%")

    rows = derived_invariants(cases, args.heading)
    print("\n  derived quantities over the pitch at z = 2 m:")
    print(f"    {'Uref':>6} {'U_mean':>8} {'sigma_u':>8} {'T_L':>8} "
          f"{'L_int':>8} {'TI':>7}")
    for r in rows:
        print(f"    {r['u_ref']:6.2f} {r['U_mean']:8.3f} {r['sigma_u']:8.3f} "
              f"{r['T_L']:8.3f} {r['L_int']:8.3f} {r['TI']:7.3f}")

    v = invariance_verdict(rows)
    print(f"\n  exponent test (predicted vs fitted, tol +/- {v['tol_exponent']}):")
    for q, d in v["quantities"].items():
        print(f"    {q:>8}: predicted {d['predicted']:+.2f}  fitted "
              f"{d['fitted']:+.4f}  r2 {d['r2']:.4f}   "
              f"{'PASS' if d['passes'] else 'FAIL'}")
    print(f"\n  VERDICT: {'PASS' if v['passes'] else 'FAIL'}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(report=rep, derived=rows, verdict=v),
                              indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

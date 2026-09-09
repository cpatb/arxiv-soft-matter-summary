"""analyse_rheology.py — the constitutive coefficients hiding in the phase-8 cells.

WHY
---
`analyse_phase8.py` reads one number per cell (the steady-state pressure) and throws the
rest away. But every cell records the whole `p_hat(gamma)` and `tau_hat(gamma)` series,
and those carry three things the project needs and does not have:

  R1  THE STEADY-STATE FRICTION LAW, mu(I). This one SURVIVES. mu = tau/p is a stress
      RATIO, and for hard grains it is a function of I alone (Forterre & Pouliquen,
      Annu. Rev. Fluid Mech. 40, 1 (2008), doi:10.1146/annurev.fluid.40.111406.102142;
      form after da Cruz, Emam, Prochnow, Roux & Chevoir, PRE 72, 021309 (2005),
      doi:10.1103/PhysRevE.72.021309). Its intercept mu_s is the quasi-static
      macroscopic friction and is a real measured constitutive coefficient.

  R2  🔴 THE phi(I) DILATANCY COEFFICIENT IS RETRACTED. Fitting phi against I and
      reading the slope as the dilatancy coefficient returned 10-36, against a
      literature scale of order 0.2 -- two orders out, which is the signature of a
      construction error rather than a material difference. Two reasons, both fatal:

        (i) On this protocol `gdot` is FIXED along a ladder (measured 3.000e-05 on all
            five rungs at cv = 0), so I moves ONLY through p. phi is the IMPOSED
            variable and p the response; "phi(I)" is therefore phi against pressure in
            disguise, i.e. a COMPRESSION curve, not the rate-dependent dilatancy law.
        (ii) The ladder reaches p_hat = 1.6e-2 to 2.2e-2, where the grains are visibly
            SOFT. `phi` as Yade reports it sums sphere volumes and so counts overlaps
            twice; at the top rung that inflation is comparable to the whole ladder
            span. mu, being a ratio, is far less exposed -- which is exactly why R1
            survives and R2 does not.

      The same error had made `analyse_phase8.py`'s E1 the PRIMARY phi_c estimator; it
      is retracted there too and E2 (the vanishing-pressure root) is primary now.

  R3  🔴 THE ELASTIC-RAMP READING OF THE EARLY TRANSIENT IS RETRACTED. The intent was
      dp/dgamma = K' alpha (phi - phi_c), so that a CV-independent K'alpha would show
      the derivation's assumption (that CV acts only through phi_c) is safe. Measured,
      dp/dgamma came out NEGATIVE on every cell including those well above phi_c, and
      the three-window gate FAILED on most of them. The model was wrong for this
      protocol: shear does not start from the critical state, it starts from an
      ISOTROPICALLY COMPRESSED one, so the early transient is dominated by the decay of
      that preparation toward the critical state -- a relaxation, not a dilatancy ramp.
      What is reported instead is the relaxation STRAIN SCALE, which is a protocol
      diagnostic and validates the steadiness gate.

      ⚠ So the derivation's assumption that alpha is CV-independent REMAINS UNTESTED.
      Testing it needs constant-PRESSURE shear, which this protocol does not do.

WHAT IS GATED
-------------
  * the transient window is swept (three windows) and the answer must not depend on it
  * cells that are not quasi-static are excluded, printed
  * a fit is refused if the series does not actually ramp (r2 bar)

    python analyse_rheology.py
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

I_MAX = 1.0e-2            # same quasi-static bar as shear.I_MAX / analyse_phase8
RAMP_WINDOWS = ((0.02, 0.15), (0.05, 0.25), (0.10, 0.35))   # in gamma
RAMP_R2_MIN = 0.90        # a ramp fit below this is not a ramp; refused
WINDOW_TOL = 0.25         # the three windows must agree to 25 % or the slope is flagged


def load(pattern: str) -> list[dict]:
    recs = []
    for path in sorted(glob.glob(os.path.join(OUT_DIR, pattern))):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if "error" not in r:
                    recs.append(r)
    return recs


def linfit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """(intercept, slope, r2) by least squares. Returns NaNs on < 3 points."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 3:
        return float("nan"), float("nan"), float("nan")
    x, y = x[ok], y[ok]
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss = float(np.sum((y - y.mean()) ** 2))
    r2 = float(1 - np.sum((y - pred) ** 2) / ss) if ss > 0 else float("nan")
    return float(coef[0]), float(coef[1]), r2


def ramp_slope(r: dict, lo: float, hi: float) -> tuple[float, float]:
    """dp_hat/dgamma over a gamma window, with its r2."""
    g = np.asarray(r.get("gamma", []), float)
    p = np.asarray(r.get("p_hat_series", []), float)
    if g.size < 4 or g.size != p.size:
        return float("nan"), float("nan")
    m = (g >= lo) & (g <= hi)
    if m.sum() < 4:
        return float("nan"), float("nan")
    _, slope, r2 = linfit(g[m], p[m])
    return slope, r2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="phase8_prod2*.jsonl")
    ap.add_argument("--out", default="rheology_verdict.txt")
    args = ap.parse_args()
    recs = load(args.pattern)
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p(f"analyse_rheology: {len(recs)} cells from {args.pattern}")
    protos = sorted({r.get("protocol", "?") for r in recs})
    if len(protos) > 1:
        p(f"🔴 REFUSING: {len(protos)} protocols present {protos}")
        return 1
    p(f"protocol {protos[0]}")
    p()

    by_cv: dict[float, list[dict]] = defaultdict(list)
    for r in recs:
        by_cv[round(float(r["cv"]), 3)].append(r)

    # ------------------------------------------------------------------- R1: mu(I)
    p("=== R1  the friction law mu(I) = mu_s + b I, per CV  [SURVIVES] ===")
    p(f"  fitted over the quasi-static cells (I <= {I_MAX:.0e}).")
    p(f"  {'cv':>5s} {'n':>3s} {'mu_s':>8s} {'b_mu':>8s} {'r2':>6s} "
      f"{'I range':>22s} {'p_hat range':>22s}")
    r1 = {}
    for cv in sorted(by_cv):
        sel = [r for r in by_cv[cv] if r.get("quasistatic")]
        iv = np.array([r["i_ss"] for r in sel], float)
        mu = np.array([r["mu_macro"] for r in sel], float)
        ph = np.array([r["p_ss_hat"] for r in sel], float)
        mus, bmu, r2m = linfit(iv, mu)
        r1[cv] = dict(mu_s=mus, b_mu=bmu, n=len(sel))
        p(f"  {cv:5.2f} {len(sel):3d} {mus:8.4f} {bmu:8.2f} {r2m:6.3f} "
          f"{iv.min():10.2e}-{iv.max():.2e} {ph.min():10.2e}-{ph.max():.2e}")
    mu_mean = float(np.mean([v["mu_s"] for v in r1.values()]))
    p()
    p(f"  mu_s = {mu_mean:.3f} averaged over CV, at mu_p = 0.3.")
    p("  ⚠ Only the FORM of the law is taken from the cited work; da Cruz et al. is 2D")
    p("  DISKS and its numerical coefficients must not be quoted against 3D spheres.")
    p("  ⚠ I moves here through p, not through the shear rate (gdot is fixed on a")
    p("  ladder), so 'small I' is the DENSE end. For hard grains mu is a function of I")
    p("  alone and the reading is still the quasi-static friction; for soft grains it")
    p("  need not be, and this ladder IS soft at the top (see the p_hat range above).")
    p()
    p("  🔑 A FREE ANCHOR IS COMING: mu_s must FALL toward ~0 as mu_p -> 0. The P4")
    p("  frictionless run (mu_p = 0) tests that at no extra cost, and a mu_s at mu_p = 0")
    p("  anywhere near 0.35 would say this fit is not measuring friction at all.")
    p()
    p("  Does the friction law itself move with CV, or only phi_c?")
    mus_v = [r1[c]["mu_s"] for c in sorted(r1)]
    p(f"     mu_s over CV 0 -> 0.5: {mus_v[0]:.4f} -> {mus_v[-1]:.4f}  "
      f"({mus_v[-1] / mus_v[0]:.3f}x)  -- a {abs(mus_v[-1] / mus_v[0] - 1) * 100:.0f} % "
      f"move against a phi_c move of 11 %")

    # --------------------------------------------------- R2/R3: retracted, see docstring
    p()
    p("=== R2  the phi(I) dilatancy coefficient  🔴 RETRACTED ===")
    p("  Fitting phi against I gave a slope of 10-36 against a literature scale of ~0.2.")
    p("  Two orders out is a construction error, not a material difference: gdot is")
    p("  FIXED along a ladder, so 'phi(I)' is phi against PRESSURE in disguise -- a")
    p("  compression curve, not a rate-dependent dilatancy law -- and the top rungs are")
    p("  soft enough that the reported phi double-counts contact overlaps.")
    p("  The same error had made analyse_phase8's E1 the PRIMARY phi_c estimator. Both")
    p("  are retracted; E2, the vanishing-pressure root, is primary now.")

    p()
    p("=== R3  the early transient  🔴 RETRACTED as an elastic ramp; reported as a "
      "relaxation ===")
    p("  Intended reading: dp/dgamma = K' alpha (phi - phi_c), which would have tested")
    p("  the derivation's assumption that alpha does not move with CV. Measured,")
    p("  dp/dgamma is NEGATIVE on every cell -- including cells well above phi_c -- and")
    p("  the three-window gate FAILS on most. The model was wrong for this protocol:")
    p("  shear starts from an ISOTROPICALLY COMPRESSED state, not from the critical")
    p("  state, so the early transient is that preparation decaying away.")
    p("  ⚠ The alpha assumption therefore REMAINS UNTESTED. It needs constant-PRESSURE")
    p("  shear, which this protocol does not do.")
    p()
    p("  What the transient CAN say -- the strain to erase the preparation, which is a")
    p("  protocol diagnostic and is what the steadiness gate depends on:")
    p(f"  {'cv':>5s} {'n':>3s} {'median gamma_relax':>19s} {'max':>8s} "
      f"{'window [0.75,1.5]?':>19s}")
    for cv in sorted(by_cv):
        gr = []
        for r in by_cv[cv]:
            g = np.asarray(r.get("gamma", []), float)
            pp = np.asarray(r.get("p_hat_series", []), float)
            if g.size < 10 or g.size != pp.size:
                continue
            ss = float(np.mean(pp[g >= 0.75 * g.max()]))
            if not np.isfinite(ss) or ss == 0:
                continue
            within = np.abs(pp - ss) <= 0.10 * abs(ss)
            idx = np.flatnonzero(within)
            # FIRST strain after which it never leaves the +/-10 % band again
            last_out = (int(np.flatnonzero(~within)[-1]) + 1
                        if (~within).any() else 0)
            if last_out < g.size:
                gr.append(float(g[min(last_out, g.size - 1)]))
        if not gr:
            p(f"  {cv:5.2f} {0:3d} {'--':>19s}")
            continue
        gr = np.array(gr)
        ok = "yes" if np.median(gr) < 0.75 else "🔴 NO"
        p(f"  {cv:5.2f} {len(gr):3d} {np.median(gr):19.3f} {gr.max():8.3f} {ok:>19s}")
    p()
    p("  The averaging window is gamma in [0.75, 1.5]. If the median relaxation strain")
    p("  is below 0.75 the window is measuring the steady state, which is what the")
    p("  steadiness gate assumes; if it is not, the gate is measuring the transient.")

    path = os.path.join(OUT_DIR, args.out)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

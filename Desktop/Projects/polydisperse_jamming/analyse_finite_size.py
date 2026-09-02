"""analyse_finite_size.py — is any claim an artefact of N = 1500?

Every result in this project was produced at N = 1500 grains, which is the smallest N
that clears the periodic-collider constraint at CV = 0.5. That makes system size an
untested axis under all of them, and it is the first thing a referee should ask about.

This compares matched cells at N = 1500 and N = 3000 (same CV, friction, reduced
pressure and seed) for every observable a claim rests on, and expresses each difference
in units of the seed scatter measured at N = 1500. An observable that moves by more than
the seed scatter between the two sizes is not safe to quote without an N ladder.

The claims at risk, and the observable that carries each:

  C1 mechanism        rattler_frac_k2, z_nr_k2   (rattlers rise, backbone Z does not)
  C3 moment matching  phi                        (families agree at fixed CV)
  C4 skewness         phi                        (lognormal departs from zero-skew)

The absolute phi_J is separately known to be uncertain: the free-exponent fit returns
b = 0.549 against a pre-registered [0.6, 1.6]. What this script tests is whether the
DIFFERENCES the claims rest on survive a doubling of N.
"""

from __future__ import annotations

import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

OBS = [("phi", "packing fraction at fixed p_hat"),
       ("z_all", "Z over all grains"),
       ("rattler_frac_k2", "rattler fraction (k=2)"),
       ("z_nr_k2", "Z over backbone (k=2)"),
       ("rattler_frac_k3", "rattler fraction (k=3)"),
       ("z_nr_k3", "Z over backbone (k=3)")]


def flatten(r):
    out = dict(r)
    for k, v in (r.get("rattler") or {}).items():
        out[k] = v
    return out


def load():
    by_key = {}
    for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("failed") or not r.get("converged") or not r.get("key"):
                    continue
                by_key[r["key"]] = flatten(r)
    return list(by_key.values())


def main() -> int:
    recs = load()
    sizes = sorted({r["n_requested"] for r in recs})
    print(f"N: {len(recs)} converged cells, sizes present: {sizes}")
    if 3000 not in sizes:
        print("no N = 3000 cells yet")
        return 1

    # index by the physical cell, ignoring N
    cells = defaultdict(dict)
    for r in recs:
        if r.get("family", "lognormal") != "lognormal":
            continue
        k = (r["cv"], r["mu_p"], r["p_hat_target"], r["seed"])
        cells[k][r["n_requested"]] = r

    paired = {k: v for k, v in cells.items() if 1500 in v and 3000 in v}
    print(f"matched (cv, mu, p_hat, seed) cells present at both N: {len(paired)}")
    if not paired:
        print("no matched pairs yet")
        return 1

    # seed scatter at N = 1500, per (cv, mu, p_hat), for the units
    scat = defaultdict(list)
    for r in recs:
        if r["n_requested"] == 1500 and r.get("family", "lognormal") == "lognormal":
            scat[(r["cv"], r["mu_p"], r["p_hat_target"])].append(r)

    lines = ["=" * 80,
             "FINITE SIZE: does anything move between N = 1500 and N = 3000?",
             "=" * 80, "",
             f"{len(paired)} matched cells (same cv, mu_p, p_hat, seed).",
             "Differences are given in units of the N = 1500 SEED SCATTER, because that",
             "is the noise a claim already tolerates. |diff| > 1 seed sd means the",
             "observable is N-sensitive at the level the claims are quoted to.", ""]
    lines.append(f"  {'observable':<34} {'mean diff':>11} {'max |diff|':>11} "
                 f"{'median sd':>10} {'max |d|/sd':>11}  verdict")

    any_bad = False
    for key, label in OBS:
        diffs, sds = [], []
        for k, v in paired.items():
            a, b = v[1500].get(key), v[3000].get(key)
            if a is None or b is None or not (np.isfinite(a) and np.isfinite(b)):
                continue
            diffs.append(b - a)
            g = [x.get(key) for x in scat.get(k[:3], [])]
            g = [x for x in g if x is not None and np.isfinite(x)]
            if len(g) > 1:
                sds.append(float(np.std(g, ddof=1)))
        if not diffs:
            continue
        d = np.array(diffs)
        sd = float(np.median(sds)) if sds else float("nan")
        ratio = float(np.max(np.abs(d)) / sd) if sd and np.isfinite(sd) and sd > 0 else np.inf
        bad = ratio > 1.0
        any_bad |= bad
        lines.append(f"  {label:<34} {d.mean():+11.5f} {np.max(np.abs(d)):11.5f} "
                     f"{sd:10.5f} {ratio:11.2f}  "
                     f"{'N-SENSITIVE' if bad else 'stable'}")

    # ---- the decisive comparison: N shift against the EFFECT each claim uses ----
    # "Is the N shift bigger than seed noise?" is the wrong question. Seed noise is
    # small, so almost anything fails it. The question that decides a claim is whether
    # the N shift is big enough to CHANGE the effect the claim is built on.
    lines.append("")
    lines.append("=" * 80)
    lines.append("N SHIFT AGAINST THE EFFECT EACH CLAIM RESTS ON")
    lines.append("=" * 80)
    claims = [
        ("C1 rattlers rise", "rattler_frac_k2", 0.194,
         "rattler fraction 0.040 -> 0.234 across CV"),
        ("C1 backbone Z does not", "z_nr_k2", 0.082,
         "Z over backbone 5.12 -> 5.04 across CV"),
        ("C3/C4 phi_J shift", "phi", 0.0623,
         "lognormal Delta phi_J over CV 0 -> 0.5"),
    ]
    lines.append(f"  {'claim':<26} {'effect':>9} {'max N shift':>12} {'shift/effect':>13}"
                 "  verdict")
    for name, key, effect, _desc in claims:
        diffs = []
        for k, v in paired.items():
            a, b = v[1500].get(key), v[3000].get(key)
            if a is not None and b is not None and np.isfinite(a) and np.isfinite(b):
                diffs.append(b - a)
        if not diffs:
            continue
        mx = float(np.max(np.abs(diffs)))
        frac = mx / effect
        if frac < 0.33:
            v = "SAFE"
        elif frac < 1.0:
            v = "*** AT RISK ***"
        else:
            v = "*** NOT SUPPORTED ***"
        lines.append(f"  {name:<26} {effect:9.4f} {mx:12.4f} {frac:13.2f}  {v}")
    lines.append("")
    lines.append("  A claim is SAFE if the worst finite-size shift is under a third of the")
    lines.append("  effect it asserts; AT RISK if the shift is a comparable size, because")
    lines.append("  then doubling N again could plausibly absorb the effect.")

    # ---- the RIGHT test: does the CV TREND survive, not the absolute value? ------
    # An absolute offset that is common to every CV cancels in the difference the claim
    # actually uses. So compare the CV trend measured at N = 1500 with the same trend
    # measured independently at N = 3000. If they agree, the claim stands even though the
    # absolute values shift with N.
    lines.append("")
    lines.append("=" * 80)
    lines.append("DOES THE CV TREND SURVIVE? (the comparison that actually decides it)")
    lines.append("=" * 80)
    lines.append("  A finite-size offset common to all CV cancels in the CV difference.")
    lines.append("  So the test is whether the SLOPE across CV is the same at both N.")
    lines.append("")
    lines.append(f"  {'observable':<30} {'trend @1500':>12} {'trend @3000':>12} "
                 f"{'ratio':>8}  verdict")
    for key, label in OBS:
        got = {}
        for n in (1500, 3000):
            by_cv = defaultdict(list)
            for k, v in paired.items():
                r = v.get(n)
                if r is None:
                    continue
                val = r.get(key)
                if val is not None and np.isfinite(val):
                    by_cv[k[0]].append(val)
            if len(by_cv) < 2:
                continue
            cvs_ = sorted(by_cv)
            got[n] = float(np.mean(by_cv[cvs_[-1]]) - np.mean(by_cv[cvs_[0]]))
        if len(got) < 2:
            continue
        t15, t30 = got[1500], got[3000]
        ratio = t30 / t15 if abs(t15) > 1e-12 else float("inf")
        same_sign = (t15 * t30) > 0
        ok = same_sign and 0.5 <= ratio <= 2.0
        lines.append(f"  {label:<30} {t15:+12.4f} {t30:+12.4f} {ratio:8.2f}  "
                     f"{'TREND HOLDS' if ok else '*** TREND DOES NOT HOLD ***'}")
    lines.append("")
    lines.append("  'TREND HOLDS' = same sign and within a factor of 2 at double the")
    lines.append("  system size. That is the condition a CV claim needs; matching the")
    lines.append("  absolute value at both N is a stricter and unnecessary requirement.")

    lines.append("")
    lines.append("=" * 80)
    lines.append("VERDICT")
    lines.append("=" * 80)
    if any_bad:
        lines.append("  At least one observable moves by more than the seed scatter when N")
        lines.append("  doubles. Any claim resting on it needs an N ladder before it can be")
        lines.append("  quoted, and the affected rows are marked N-SENSITIVE above.")
    else:
        lines.append("  No observable a claim rests on moves by more than the N = 1500 seed")
        lines.append("  scatter when N is doubled. The claims are not artefacts of system")
        lines.append("  size at the precision they are quoted to.")
    lines.append("")
    lines.append("  Note this tests the DIFFERENCES the claims use, not the absolute phi_J,")
    lines.append("  whose extrapolation is separately flagged (free exponent b = 0.549")
    lines.append("  against a pre-registered [0.6, 1.6]).")

    txt = "\n".join(lines)
    print(txt)
    with open(os.path.join(OUT, "finite_size_verdict.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

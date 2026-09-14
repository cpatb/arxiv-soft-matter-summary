"""grade_mixture_mu0.py -- grade the FRICTIONLESS block, PREREG_MIXTURE_2026-09-13.md 4c.

This is a separate script from `grade_mixture.py` on purpose. M0-M6 are CLOSED on the
mu_p = 0.3 grid; F1-F3 are a different question asked of a different block, and keeping
them in one file is how a closed gate gets silently re-decided.

The frictionless block is the only genuinely BLIND test in this campaign: no binary
packing at mu_p = 0 existed anywhere in this project when the predictions in section 4c
were written, and the model's claim is that the interaction functions are GEOMETRIC, so
changing beta from 0.5947 to the measured 0.6413 and touching nothing else predicts the
whole composition curve.
"""

from __future__ import annotations

import collections
import glob
import json
import os
import sys

import numpy as np

import mixture_closure as mc
import mixture_model as mm
from phase12_mixture import MIN_PHASE_COUNT, N_PROD, X_ANCHOR
from psd import BinaryPSD

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU0 = 0.0
MU_REF = 0.3

F1_PASS_SCALES = 5.0
F1_FAIL_SCALES = 2.0
F2_RATIO_BAR = 3.0            # per-point chi2 at mu 0 must be within this x the mu 0.3 one
PHI_DENSE_BOUND = 0.72        # packing.PHI_DENSE_BOUND, the collider pre-flight assumption


def load(mu: float):
    g: dict[tuple[float, float], list[dict]] = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(OUT, "phase12_mixr*_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("error"):
                    g[(float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))].append(r)
                    continue
                if r.get("n_requested") != N_PROD:
                    continue
                if abs(float(r.get("mu_p", -1)) - mu) > 1e-9:
                    continue
                g[(float(r["mix_ratio"]), round(float(r["x_small_vol"]), 6))].append(r)
    return g


def gradable(key: tuple[float, float]) -> bool:
    if key[1] in (0.0, 1.0):
        return True
    p = BinaryPSD(d10=1.0, ratio=key[0], x_small=key[1])
    return min(p.n_small, 1.0 - p.n_small) * N_PROD >= MIN_PHASE_COUNT


def main() -> int:
    lines: list[str] = []
    failures: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    g0 = load(MU0)
    g3 = load(MU_REF)
    errs = [c for v in g0.values() for c in v if "error" in c]
    nonconv = [c for v in g0.values() for c in v
               if "error" not in c and not c.get("converged")]
    say("=== grade_mixture_mu0.py -- PREREG_MIXTURE_2026-09-13.md section 4c ===")
    say(f"   {sum(len(v) for v in g0.values())} frictionless cells, {len(errs)} errors, "
        f"{len(nonconv)} non-converged. scale = {SCALE:.1e}")
    for c in errs + nonconv:
        say(f"   PROBLEM {c.get('key')}: {str(c.get('error') or c.get('notes'))[:140]}")
    say()

    meas3 = mm.measured(mu_p=MU_REF)
    meas0 = mm.measured(mu_p=MU0)
    beta3 = meas3[("lognormal", 0.0)][0]
    beta0 = meas0[("lognormal", 0.0)][0]
    say("=== beta, the ONLY thing that changes between the two blocks ===")
    say(f"   beta(mu_p = {MU_REF}) = {beta3:.4f}   beta(mu_p = {MU0}) = {beta0:.4f}   "
        f"ratio {beta0/beta3:.4f}")
    ends = [float(np.mean([c["phi"] for c in g0[(r, x)]]))
            for r in (2.0, 3.0) for x in (0.0, 1.0) if (r, x) in g0]
    say("   the block's own monodisperse endpoints: "
        + ", ".join(f"{e:.4f}" for e in ends))
    say(f"   -> they agree with beta(mu_p = 0) to "
        f"{max(abs(e - beta0) for e in ends)/SCALE:.1f} scale, which is the control that")
    say("      the frictionless block is at the friction it says it is.")
    say()

    pa, pb, _ = mc.fit_on_binaries(mc.binary_points(mu_p=MU_REF), beta3)
    say(f"=== the prediction: the ALL-RATIO calibration (PA = {pa:.3f}, PB = {pb:.3f}),")
    say(f"    fitted on mu_p = {MU_REF} binaries, evaluated with beta = {beta0:.4f} ===")
    say(f"   {'r':>4}{'X_s':>8}{'phi_meas':>10}{'sem':>8}{'pred':>9}{'/scale':>8}"
        f"{'n_min':>7}{'grad':>6}")
    chi0 = 0.0
    n0 = 0
    over = []
    for key in sorted(g0):
        cells = [c for c in g0[key] if "error" not in c]
        if not cells:
            continue
        v = np.array([float(c["phi"]) for c in cells])
        sem = float(v.std(ddof=1) / np.sqrt(v.size)) if v.size > 1 else float("nan")
        p = BinaryPSD(d10=1.0, ratio=key[0], x_small=key[1])
        nmin = min(p.n_small, 1.0 - p.n_small) * N_PROD
        pred = mm.binary_density(key[0], key[1], beta0, pa, pb)
        ok = gradable(key)
        if ok and key[1] not in (0.0, 1.0):
            chi0 += ((pred - v.mean()) / SCALE) ** 2
            n0 += 1
        if v.mean() > PHI_DENSE_BOUND:
            over.append((key, float(v.mean())))
        say(f"   {key[0]:4.1f}{key[1]:8.4f}{v.mean():10.4f}{sem:8.4f}{pred:9.4f}"
            f"{(pred-v.mean())/SCALE:+8.1f}{nmin:7.0f}{'y' if ok else 'NO':>6}")
    say()

    # ---- F1 --------------------------------------------------------------------------
    say("=== F1 -- the matched pair, at ZERO friction ===")
    a_key, b_key = (3.0, round(X_ANCHOR, 6)), (3.0, 0.25)
    if a_key in g0 and b_key in g0:
        va = np.array([float(c["phi"]) for c in g0[a_key] if "error" not in c])
        vb = np.array([float(c["phi"]) for c in g0[b_key] if "error" not in c])
        d = float(vb.mean() - va.mean())
        sd = float(np.hypot(va.std(ddof=1) / np.sqrt(va.size),
                            vb.std(ddof=1) / np.sqrt(vb.size)))
        pred = (mm.binary_density(3.0, 0.25, beta0, pa, pb)
                - mm.binary_density(3.0, X_ANCHOR, beta0, pa, pb))
        cva = float([c for c in g0[a_key] if "error" not in c][0]["cv_realised"])
        cvb = float([c for c in g0[b_key] if "error" not in c][0]["cv_realised"])
        say(f"   both at size ratio 3 and number CV {cva:.6f} / {cvb:.6f} "
            f"(agreeing to {abs(cva-cvb):.1e})")
        say(f"   phi = {va.mean():.4f} and {vb.mean():.4f}")
        say(f"   MEASURED dphi = {d:+.4f} +- {sd:.4f}  ({d/SCALE:+.1f} scale, {d/sd:+.1f} se)")
        say(f"   REGISTERED PREDICTION (4c) = +0.0684; this run's calibration gives {pred:+.4f}")
        say("   descriptor account = +0.0000")
        if abs(d) > F1_PASS_SCALES * SCALE and d > 0:
            say(f"   -> F1 PASS ({abs(d)/SCALE:.1f} scale against a bar of {F1_PASS_SCALES}); "
                f"the prediction is {abs(pred/d - 1)*100:.0f}% off")
        elif abs(d) < F1_FAIL_SCALES * SCALE:
            say("   -> F1 FAIL")
            failures.append("F1")
        else:
            say("   -> F1 UNDETERMINED")
            failures.append("F1(undetermined)")
    say()

    # ---- F2 --------------------------------------------------------------------------
    say("=== F2 -- are the interaction functions GEOMETRIC? ===")
    # 🔴 The reference MUST be restricted to the compositions the frictionless block
    # actually has. The first version of this looped over every mu_p = 0.3 key, which
    # includes the ratio 1.25 and 1.5 ladder -- 32 points against the blind set's 16 --
    # while the comment beside it claimed "same compositions". Two arms of a comparison
    # taken on different grids is not a comparison (lessons.md, 2026-09-05). Both
    # numbers are printed so the difference the restriction makes is visible.
    blind_keys = {k for k in g0 if gradable(k) and k[1] not in (0.0, 1.0)}

    def ref_chi(keys):
        s, n = 0.0, 0
        for key in sorted(keys):
            cells = [c for c in g3.get(key, []) if "error" not in c]
            if not cells:
                continue
            m = float(np.mean([c["phi"] for c in cells]))
            s += ((mm.binary_density(key[0], key[1], beta3, pa, pb) - m) / SCALE) ** 2
            n += 1
        return s, n

    chi3, n3 = ref_chi(blind_keys)
    chi3_all, n3_all = ref_chi({k for k in g3 if gradable(k) and k[1] not in (0.0, 1.0)})
    say("   per-point chi2, same calibration, MATCHED compositions:")
    say(f"      mu_p = {MU_REF}: {chi3/max(n3,1):8.2f}  over {n3} points "
        f"(this is the fit's OWN set, so it is a floor, not a comparison)")
    say(f"         [over all {n3_all} mu_p = {MU_REF} compositions instead it would be "
        f"{chi3_all/max(n3_all,1):.2f}, i.e. the restriction matters]")
    say(f"      mu_p = {MU0}: {chi0/max(n0,1):8.2f}  over {n0} points  -- BLIND")
    ratio = (chi0 / max(n0, 1)) / max(chi3 / max(n3, 1), 1e-12)
    say(f"   ratio = {ratio:.2f} against a registered bar of {F2_RATIO_BAR}")
    if ratio <= F2_RATIO_BAR:
        say("   -> F2 PASS: only beta carried the friction, at this bar. The interaction")
        say("      functions transfer across a friction change that moves beta by "
            f"{(beta0-beta3)/SCALE:.0f} scale.")
    else:
        say("   -> F2 FAIL: the interaction functions are friction-DEPENDENT, so the")
        say("      closure is a mu_p = 0.3 object and must be quoted as one.")
        failures.append("F2")
    say("   /!\\ the mu_p = 0.3 column is the calibration set, so it is a floor rather")
    say("       than a fair comparison; the bar is 3x precisely to allow for that, and")
    say("       the blind number is the one that means anything.")
    say()

    # ---- F3 --------------------------------------------------------------------------
    say("=== F3 -- the collider guard, registered in advance ===")
    if over:
        say(f"   {len(over)} cells exceed packing.PHI_DENSE_BOUND = {PHI_DENSE_BOUND}:")
        for k, v in over:
            say(f"      r = {k[0]:.0f}, X_s = {k[1]:.2f}: phi = {v:.4f} "
                f"(+{v-PHI_DENSE_BOUND:.4f})")
        say("   The guard sizes the cell for the collider at the DENSEST state, so it was")
        say("   optimistic by this much. d_max/L scales as phi^(1/3), so the worst cell's")
        worst = max(v for _, v in over)
        say(f"   d_max/L is understated by {(worst/PHI_DENSE_BOUND)**(1/3)-1:.2%} -- it")
        say("   cannot invalidate a cell, and it was written down before the run.")
    else:
        say(f"   no cell exceeds {PHI_DENSE_BOUND}. The guard held.")
    say()

    say(f"FAILURES: {failures}")
    with open(os.path.join(OUT, "mixture_mu0.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

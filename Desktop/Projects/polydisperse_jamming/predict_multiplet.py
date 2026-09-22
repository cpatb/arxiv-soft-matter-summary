"""predict_multiplet.py -- BLIND predictions for the multiplet DEM cells.

Run BEFORE any multiplet cell exists. Emits, for each registered composition:
  * the family string the runner takes,
  * the realised CV, size ratio and Lambda (all closed form),
  * the calibrated closure's predicted phi at the campaign rung,
  * the same under the TEXTBOOK exponents, so the prediction's fit-dependence is visible.

THE DESIGN. Three blocks, and the third is the one that can fail informatively.

  A  LADDER IN d_f at fixed composition. Opens the aggregates up without changing how
     many are aggregated. Tests whether the closure tracks a PSD whose size ratio it has
     never seen from this direction.

  B  LADDER IN <k> at fixed d_f. Aggregates more without changing compactness.

  C  MATCHED-CV PAIRS. Two aggregation states with the SAME number CV to 3 decimals and
     DIFFERENT d_f. This project's manuscript already shows that two packs matched in
     width descriptors can differ in phi_J; block C asks the same question along the
     AGGREGATION axis, where the pair is reachable by a real bench protocol (salt
     changes d_f, time changes <k>). A CV-only closure predicts ZERO difference for
     every pair in block C, which is what makes it a discriminator rather than a fit.

Run:  python3 predict_multiplet.py     (writes out/predict_multiplet.txt)
"""

from __future__ import annotations

import itertools
import os
import sys

import numpy as np

import mixture_model as mm
from multiplet_levers import BETA, D10, PA_REFIT, PA_TEXT, PB_REFIT, PB_TEXT
from psd import MultipletPSD

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SCALE = mm.PHI_SCALE

# The campaign's own grid. CV is DERIVED for this family, so the axes are (d_f, n_k).
BLOCK_A = [(1.8, (0.70, 0.20, 0.10)), (2.1, (0.70, 0.20, 0.10)),
           (2.5, (0.70, 0.20, 0.10)), (3.0, (0.70, 0.20, 0.10))]
BLOCK_B = [(2.1, (0.90, 0.10)), (2.1, (0.70, 0.20, 0.10)),
           (2.1, (0.50, 0.30, 0.20)), (2.1, (0.30, 0.30, 0.40))]


def fam_string(d_f: float, n_k) -> str:
    return "multiplet:f" + f"{d_f:g}" + ":n" + "/".join(f"{x:g}" for x in n_k)


def describe(d_f: float, n_k) -> dict:
    p = MultipletPSD(d10=D10, d_f=d_f, n_k=tuple(n_k))
    return dict(
        fam=fam_string(d_f, n_k), d_f=d_f, n_k=tuple(n_k),
        cv=p.cv, ratio=p.size_ratio, lam=p.lambda_trap, k_bar=p.k_bar,
        phi_refit=mm.packing_density(p.diameters, p.volume_fractions,
                                     BETA, PA_REFIT, PB_REFIT),
        phi_text=mm.packing_density(p.diameters, p.volume_fractions,
                                    BETA, PA_TEXT, PB_TEXT),
    )


def matched_cv_pairs(tol: float = 5e-4, min_df_gap: float = 0.5,
                     min_dphi_scales: float = 2.0) -> list[tuple[dict, dict]]:
    """Find (d_f, n_k) pairs with the SAME CV and DIFFERENT d_f.

    Selection is on CV agreement and on `d_f` separation ONLY. The predicted `phi`
    difference is then REPORTED, and pairs whose predicted difference is below
    `min_dphi_scales * PHI_SCALE` are dropped as untestable -- not because they are
    uninteresting but because the DEM protocol could not resolve them, so a cell there
    could neither confirm nor refute anything.

    ⚠ The pair is chosen by CV, never by the size of the predicted difference; choosing
    on the response would be choosing the design from the answer.
    """
    grid = []
    for d_f in (1.6, 1.8, 2.1, 2.5, 3.0):
        for n2 in np.arange(0.05, 0.55, 0.05):
            for n3 in np.arange(0.0, 0.35, 0.05):
                n1 = 1.0 - n2 - n3
                if n1 < 0.30:
                    continue
                nk = (round(n1, 3), round(float(n2), 3), round(float(n3), 3))
                grid.append(describe(d_f, nk))
    out = []
    for a, b in itertools.combinations(grid, 2):
        if abs(a["d_f"] - b["d_f"]) < min_df_gap:
            continue
        if abs(a["cv"] - b["cv"]) > tol:
            continue
        if abs(a["phi_refit"] - b["phi_refit"]) < min_dphi_scales * SCALE:
            continue
        out.append((a, b))
    # keep the widest d_f separation per CV bucket, deterministically
    out.sort(key=lambda p: (-abs(p[0]["d_f"] - p[1]["d_f"]),
                            -abs(p[0]["phi_refit"] - p[1]["phi_refit"])))
    seen, keep = set(), []
    for a, b in out:
        bucket = round(a["cv"], 2)
        if bucket in seen:
            continue
        seen.add(bucket)
        keep.append((a, b))
    return keep[:4]


def main() -> int:
    L: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        L.append(s)

    say("=" * 100)
    say("predict_multiplet.py -- BLIND predictions, written before any multiplet cell exists")
    say("=" * 100)
    say(f"   closure: beta = {BETA} (measured), refit PA = {PA_REFIT}, PB = {PB_REFIT};")
    say(f"   textbook PA = {PA_TEXT}, PB = {PB_TEXT}. Rung p_hat = {mm.P_HAT}, "
        f"mu_p = {mm.MU_P}, N = {mm.N_REQ}.")
    say(f"   PHI_SCALE = {SCALE}. A predicted difference below 1 scale is NOT testable.")

    for name, block in (("A  ladder in d_f, composition fixed at 70/20/10", BLOCK_A),
                        ("B  ladder in <k>, d_f fixed at 2.1", BLOCK_B)):
        say(f"\n=== BLOCK {name} ===")
        say(f"   {'family string':>42}{'<k>':>7}{'CV':>8}{'ratio':>8}{'Lambda':>9}"
            f"{'phi refit':>11}{'phi text':>10}")
        for d_f, n_k in block:
            d = describe(d_f, n_k)
            say(f"   {d['fam']:>42}{d['k_bar']:7.3f}{d['cv']:8.4f}{d['ratio']:8.3f}"
                f"{d['lam']:9.4f}{d['phi_refit']:11.4f}{d['phi_text']:10.4f}")
        ph = [describe(f, n)["phi_refit"] for f, n in block]
        say(f"   span of the predicted phi across this block: "
            f"{max(ph)-min(ph):.4f} = {(max(ph)-min(ph))/SCALE:.1f} x PHI_SCALE")

    say("\n=== BLOCK C  MATCHED-CV PAIRS: same width, different compactness ===")
    say("   A closure written in CV alone predicts ZERO difference for every row here.")
    say("")
    pairs = matched_cv_pairs()
    if not pairs:
        say("   🔴 NO PAIR FOUND that is both CV-matched and predicted-resolvable.")
        say("   That is itself a result: it would mean the aggregation axis cannot")
        say("   produce a matched-CV contrast this protocol can see.")
    for i, (a, b) in enumerate(pairs, 1):
        d_phi = a["phi_refit"] - b["phi_refit"]
        say(f"   PAIR {i}:  CV = {a['cv']:.4f} vs {b['cv']:.4f}  "
            f"(matched to {abs(a['cv']-b['cv']):.1e})")
        for d in (a, b):
            say(f"      {d['fam']:>42}  d_f = {d['d_f']:.1f}  <k> = {d['k_bar']:.3f}  "
                f"ratio = {d['ratio']:.3f}  Lambda = {d['lam']:.4f}  "
                f"phi = {d['phi_refit']:.4f}")
        say(f"      PREDICTED d(phi) = {d_phi:+.4f} = {d_phi/SCALE:+.1f} x PHI_SCALE"
            f"   (textbook: {(a['phi_text']-b['phi_text'])/SCALE:+.1f})")
        say("")

    say("=== WHAT EACH BLOCK CAN AND CANNOT DECIDE ===")
    say("   A and B test the closure's TRANSFER to a PSD family no point of which is in")
    say("   its calibration set (16 binary compositions at size ratios 2 and 3). They")
    say("   cannot test Lambda: the DEM packs effective spheres, so it never sees the")
    say("   solid/excluded volume distinction Lambda is built from.")
    say("   C tests whether CV suffices. It cannot distinguish 'the closure is right'")
    say("   from 'some other non-CV descriptor is doing the work' -- only that a")
    say("   CV-only description is incomplete along the aggregation axis.")
    say("   NOTHING here tests assumption A1 (the effective-sphere replacement). That")
    say("   needs CLUMP cells at the same n_k, and it is registered separately.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "predict_multiplet.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    say("   written: out/predict_multiplet.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())

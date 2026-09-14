"""grade_ternary_N.py -- grade NA0-NA4 of PREREG_TERNARY_N_2026-09-14.md.

The discriminator T1 could not run: does the ternary under-prediction track COMPOSITION
(non-additivity) or the number of COARSE GRAINS (finite phase size)? Those two are
confounded by construction at fixed N, because number fraction goes as y_i/d_i^3.

This block holds every composition FIXED and doubles N, which moves n_L by 2x and moves
nothing else.

It loads `phase12_tern3k_shard*.jsonl` for the N = 3000 arm and
`phase12_tern_shard*.jsonl` for the N = 1500 arm. It cannot touch `mixr*`/`extr*`, and
`grade_ternary.py`'s own gates are not recomputed here -- T1 keeps its verdict at
N = 1500 whatever this says; NA1 decides what that verdict MEANS.

    python grade_ternary_N.py      -> out/grade_ternary_N.txt
"""

from __future__ import annotations

import collections
import glob
import json
import os

import numpy as np
from scipy.stats import spearmanr

import grade_mixture_ext as gx
import mixture_closure as mc
import mixture_model as mm
from phase12_ternary import COMPOSITIONS, MONO, make_ternary

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU = 0.3

# --- registered bars, PREREG_TERNARY_N_2026-09-14.md section 3 -----------------------
NA1_NONADD_BAR = 4.7     # scale; 2/3 of the N=1500 mean residual +7.53
NA1_FINITE_BAR = 2.8     # scale; the T1 bar -- at or below this, T1's verdict is withdrawn
NA2_K = 2.0              # multiples of the se of the difference of two 3-seed means
NA3_RHO_BAR = 0.6        # |pooled spearman| above which n_L is the better explanation
N_SMALL, N_BIG = 1500, 3000

LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s)


def load(pattern: str, n_req: int):
    g: dict[tuple, list[dict]] = collections.defaultdict(list)
    nraw = 0
    for f in sorted(glob.glob(os.path.join(OUT, pattern))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                nraw += 1
                r = json.loads(line)
                if r.get("error") or not r.get("converged"):
                    continue
                if r.get("n_requested") != n_req:
                    continue
                if abs(float(r.get("mu_p", -1)) - MU) > 1e-9:
                    continue
                g[(round(float(r["tern_y_L"]), 4), round(float(r["tern_y_M"]), 4),
                   round(float(r["tern_y_S"]), 4))].append(r)
    return g, nraw


def block_stats(g, beta, fa, fb, n_req):
    """(residual, mean phi, sd, n, n_L) per composition, against the block's OWN beta."""
    out = {}
    for comp in COMPOSITIONS:
        cs = g.get(tuple(round(v, 4) for v in comp), [])
        if not cs:
            continue
        v = np.array([c["phi"] for c in cs], dtype=float)
        p = make_ternary(comp)
        pred = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        out[comp] = {"res": (float(v.mean()) - pred) / SCALE,
                     "phi": float(v.mean()),
                     "sd": float(v.std(ddof=1)) if v.size > 1 else float("nan"),
                     "n": int(v.size),
                     "n_L": int(p.class_counts(n_req)[2]),
                     "pred": pred}
    return out


def main() -> int:
    meas = mm.measured(mu_p=MU)
    beta1500_iso = meas[("lognormal", 0.0)][0]
    allpts = gx.ext_points(beta1500_iso)
    tab, fa, fb = gx.tabulated_closure(beta1500_iso, allpts)

    say("=== grade_ternary_N.py -- PREREG_TERNARY_N_2026-09-14.md ===")
    say("    The closure (a, b) is the SAME tabulated one, calibrated on N = 1500")
    say("    binaries. Only beta is re-read per block. Nothing is refitted.")
    say()

    g15, nraw15 = load("phase12_tern_shard*.jsonl", N_SMALL)
    g30, nraw30 = load("phase12_tern3k_shard*.jsonl", N_BIG)
    say("=== data ===")
    say(f"   N = {N_SMALL}: {nraw15} records read, "
        f"{sum(len(v) for v in g15.values())} admitted, {len(g15)} compositions")
    say(f"   N = {N_BIG}: {nraw30} records read, "
        f"{sum(len(v) for v in g30.values())} admitted, {len(g30)} compositions")
    if not g30:
        say("   *** the N = 3000 block is absent; nothing below can be graded ***")
        with open(os.path.join(OUT, "grade_ternary_N.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(LINES) + "\n")
        return 1
    say()

    # ---------------------------------------------------------------- NA0
    say("=== NA0 -- beta at each N, re-measured from each block's OWN control ===")
    betas = {}
    for lbl, g, nreq in ((N_SMALL, g15, N_SMALL), (N_BIG, g30, N_BIG)):
        cs = g.get(tuple(round(v, 4) for v in MONO), [])
        v = np.array([c["phi"] for c in cs], dtype=float)
        betas[lbl] = float(v.mean())
        say(f"   N = {lbl}: monodisperse control phi = {betas[lbl]:.6f} "
            f"+- {float(v.std(ddof=1)):.6f} (n {v.size})")
    say(f"   isotropic lognormal CV 0 at N = 1500 (the closure's calibration beta): "
        f"{beta1500_iso:.6f}")
    say(f"   finite-size shift in beta itself, N {N_SMALL} -> {N_BIG}: "
        f"{(betas[N_BIG] - betas[N_SMALL]) / SCALE:+.2f} scale. Using the N = 1500 beta")
    say("   for the N = 3000 arm would fold this into the quantity under test, so each")
    say("   arm is graded against its own.")
    say()

    s15 = block_stats(g15, betas[N_SMALL], fa, fb, N_SMALL)
    s30 = block_stats(g30, betas[N_BIG], fa, fb, N_BIG)

    # ---------------------------------------------------------------- NA1/NA2
    say("=== NA1 / NA2 -- the residual at doubled N, composition held fixed ===")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}"
        f"{'n_L@1500':>10}{'res@1500':>10}{'n_L@3000':>10}{'res@3000':>10}"
        f"{'d(res)':>9}{'se(d)':>8}{'in bar?':>9}")
    na2_pass, na2_tot = 0, 0
    r15, r30 = [], []
    for comp in COMPOSITIONS:
        if comp not in s15 or comp not in s30:
            say(f"   {comp} missing from one arm")
            continue
        a, b = s15[comp], s30[comp]
        r15.append(a["res"]); r30.append(b["res"])
        # se of the difference of two independent three-seed means, in scale
        se = float(np.sqrt((a["sd"] / np.sqrt(a["n"])) ** 2
                           + (b["sd"] / np.sqrt(b["n"])) ** 2)) / SCALE
        d = b["res"] - a["res"]
        ok = abs(d) <= NA2_K * se
        na2_tot += 1
        na2_pass += int(ok)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}"
            f"{a['n_L']:10d}{a['res']:+10.2f}{b['n_L']:10d}{b['res']:+10.2f}"
            f"{d:+9.2f}{se:8.2f}{'yes' if ok else 'NO':>9}")
    m15, m30 = float(np.mean(r15)), float(np.mean(r30))
    say(f"   mean residual: N {N_SMALL} {m15:+.2f} scale;  N {N_BIG} {m30:+.2f} scale")
    say(f"   rms  residual: N {N_SMALL} {float(np.sqrt(np.mean(np.square(r15)))):.2f};  "
        f"N {N_BIG} {float(np.sqrt(np.mean(np.square(r30)))):.2f}")
    say(f"   sign structure at N {N_BIG}: {int(np.sum(np.array(r30) > 0))}/{len(r30)} positive")
    say()
    if m30 >= NA1_NONADD_BAR:
        na1 = "NON-ADDITIVE"
    elif m30 <= NA1_FINITE_BAR:
        na1 = "FINITE-SIZE ARTIFACT -- T1's verdict is WITHDRAWN"
    else:
        na1 = "UNDECIDED"
    say(f"   -> NA1 verdict: **{na1}**  (mean {m30:+.2f} against "
        f"non-additive >= +{NA1_NONADD_BAR}, finite-size <= +{NA1_FINITE_BAR})")
    say(f"   -> NA2: {na2_pass}/{na2_tot} compositions agree between the two N within "
        f"{NA2_K} se of the difference")
    say()

    # ---------------------------------------------------------------- NA3
    say("=== NA3 -- the correlation that raised the alarm, recomputed POOLED ===")
    say("   If the residual is set by n_L, doubling n_L must move it and the pooled")
    say("   correlation survives. If it is set by composition, the N = 3000 points sit")
    say("   at the SAME residual with DOUBLE the n_L, which destroys the correlation.")
    say("   `d_max/L` is read from the SOLVER's own record, not recomputed: it is the")
    say("   second face of the same confound (spearman +1.0000 with the residual at")
    say("   N = 1500) and doubling N multiplies it by 2^(-1/3) = 0.794.")
    nl_all, res_all, lbl_all, dm_all = [], [], [], []
    for comp in COMPOSITIONS:
        for tag, s, g in ((N_SMALL, s15, g15), (N_BIG, s30, g30)):
            if comp in s:
                nl_all.append(s[comp]["n_L"]); res_all.append(s[comp]["res"])
                lbl_all.append(tag)
                cs = g[tuple(round(v, 4) for v in comp)]
                dm_all.append(float(np.mean([c.get("d_max_over_L", np.nan) for c in cs])))
    rho_dm, p_dm = spearmanr(dm_all, res_all)
    rho_pool, p_pool = spearmanr(nl_all, res_all)
    rho15, p15 = spearmanr([s15[c]["n_L"] for c in COMPOSITIONS if c in s15],
                           [s15[c]["res"] for c in COMPOSITIONS if c in s15])
    rho30, p30 = spearmanr([s30[c]["n_L"] for c in COMPOSITIONS if c in s30],
                           [s30[c]["res"] for c in COMPOSITIONS if c in s30])
    say(f"   within N = {N_SMALL}: spearman(res, n_L) = {rho15:+.4f} (p {p15:.4f}, n 6)")
    say(f"   within N = {N_BIG}: spearman(res, n_L) = {rho30:+.4f} (p {p30:.4f}, n 6)")
    say(f"   POOLED over {len(nl_all)} points:      = {rho_pool:+.4f} (p {p_pool:.4f})")
    say(f"   POOLED spearman(res, d_max/L) over {len(dm_all)} points = {rho_dm:+.4f} "
        f"(p {p_dm:.4f}); d_max/L spans "
        f"{min(dm_all):.4f}-{max(dm_all):.4f} across both N")
    say(f"   -> NA3: {'n_L survives as the explanation' if abs(rho_pool) > NA3_RHO_BAR else 'n_L is EXCLUDED as the explanation'} "
        f"(bar |rho| {NA3_RHO_BAR})")
    say("   NOTE the within-block correlation is expected to stay strong in BOTH blocks")
    say("   whatever the answer, because n_L and y_S are collinear at fixed N by")
    say("   construction. Only the POOLED number can separate them, which is the whole")
    say("   reason the second N was run.")
    say()

    # ---------------------------------------------------------------- NA4
    say("=== NA4 -- the surviving candidate, within the N = 3000 block alone ===")
    for nm, v in (("y_S", [c[2] for c in COMPOSITIONS if c in s30]),
                  ("y_M", [c[1] for c in COMPOSITIONS if c in s30]),
                  ("y_L", [c[0] for c in COMPOSITIONS if c in s30]),
                  ("CV", [make_ternary(c).cv for c in COMPOSITIONS if c in s30]),
                  ("d32", [make_ternary(c).d32 for c in COMPOSITIONS if c in s30])):
        rho, p_ = spearmanr(v, [s30[c]["res"] for c in COMPOSITIONS if c in s30])
        say(f"   spearman(residual@3000, {nm:>4}) = {rho:+.4f}  (p {p_:.4f})")
    say()

    with open(os.path.join(OUT, "grade_ternary_N.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    print(f"-> {os.path.join(OUT, 'grade_ternary_N.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""predict_ternary.py -- the BLIND predictions of PREREG_TERNARY_2026-09-14.md.

Run BEFORE any ternary cell exists. It loads only the BINARY blocks (`phase12_mixr*`
and `phase12_extr*`) and the monodisperse `beta`, and writes the number each candidate
closure commits to for each registered ternary composition.

Nothing here can see a ternary measurement, because no ternary record exists when it is
run and the loaders it uses do not match the `tern` tag. That is the whole point: the
predictions go into the pre-registration, and `grade_ternary.py` compares them to the
data without recomputing them.

    python predict_ternary.py        -> out/predict_ternary.txt
"""

from __future__ import annotations

import os

import numpy as np

import grade_mixture_ext as gx
import mixture_closure as mc
import mixture_model as mm
from phase12_ternary import COMPOSITIONS, D_REL, MONO, make_ternary, y_string
from phase12_mixture import MIN_PHASE_COUNT

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
N_PROD = 1500
MU = 0.3

LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s)


def main() -> int:
    say("=== predict_ternary.py -- blind predictions for PREREG_TERNARY_2026-09-14.md ===")
    say()

    # --- the calibration inputs, all binary ------------------------------------------
    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)                 # ratios <= 3
    allpts = gx.ext_points(beta)                      # + ratios 4 and 5
    say(f"   beta = {beta:.6f}  (measured monodisperse phi, mu_p {MU}, p_hat "
        f"{mm.P_HAT:.0e}, N {mm.N_REQ}; n = {meas[('lognormal', 0.0)][2]} cells, "
        f"sem {meas[('lognormal', 0.0)][1]:.2e})")
    say(f"   PHI_SCALE = {SCALE:.1e} (declared in mixture_model.py, the protocol's own "
        f"reproducibility -- NOT the seed sem)")
    say(f"   binary calibration points: {len(pts23)} at ratios "
        f"{sorted({p[0] for p in pts23})} (mixr block)")
    say(f"                              {len(allpts)} at ratios "
        f"{sorted({p[0] for p in allpts})} (mixr + extr blocks)")
    say()

    cand, per, (ia, sa, ib, sb) = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    say(f"   J   (the SHIPPED two-exponent closure, fitted on the mixr binaries): "
        f"PA = {paJ:.4f}, PB = {pbJ:.4f}")
    tab, fa, fb = gx.tabulated_closure(beta, allpts)
    say("   TAB (the tabulated/nonparametric closure of grade_mixture_ext section 5),")
    say("        identifiable (a, b) at u = "
        + ", ".join(f"{u:.4f} -> ({v[0]:.4f}, {v[1]:.4f})" for u, v in sorted(tab.items()))
        + "; plus the proved limits (0,0) at u=0 and (1,1) at u=1, joined by PCHIP")
    say()

    # --- the three pair ratios a ternary at d = 1,2,4 actually asks about -------------
    say("=== the pairs this ladder exercises (u = d_fine/d_coarse) ===")
    say(f"   {'pair':>10}{'u':>8}{'a_J':>9}{'b_J':>9}{'a_TAB':>9}{'b_TAB':>9}"
        f"{'|da|':>8}{'|db|':>8}")
    pairs = [(D_REL[0], D_REL[1]), (D_REL[1], D_REL[2]), (D_REL[0], D_REL[2])]
    for df, dc in pairs:
        u = df / dc
        aJ = float(np.sqrt(1 - (1 - u) ** paJ))
        bJ = float(1 - (1 - u) ** pbJ)
        aT, bT = float(fa(u)), float(fb(u))
        say(f"   {f'{df:.0f}-{dc:.0f}':>10}{u:8.4f}{aJ:9.4f}{bJ:9.4f}{aT:9.4f}{bT:9.4f}"
            f"{abs(aJ - aT):8.4f}{abs(bJ - bT):8.4f}")
    say("   Both u values lie INSIDE the binary calibration range (u = 0.20 to 0.80), so")
    say("   the interaction VALUES are interpolated, not extrapolated. What is new is")
    say("   that three pairs act at once. This is the additivity question and nothing else.")
    say()

    # --- admissibility of each registered composition ---------------------------------
    say("=== admissibility (registered BEFORE the run, not applied after it) ===")
    say(f"   MIN_PHASE_COUNT = {MIN_PHASE_COUNT} grains, inherited verbatim from "
        f"phase12_mixture.py")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'n_L':>6}{'n_M':>6}{'n_S':>6}"
        f"{'d_L/d_S':>9}{'CV':>7}{'d32':>8}{'d_max/L':>9}  family string")
    refused = []
    for comp in list(COMPOSITIONS) + [(0.50, 0.30, 0.20)]:
        p = make_ternary(comp)
        counts = p.class_counts(N_PROD)
        # d_max/L at the DENSE state, the same diagnostic packing.py flags at 0.25.
        # L^3 = (pi/6) sum n_i d_i^3 / phi, with phi taken at beta (a lower bound on the
        # ternary phi, so this is the WORST case, i.e. the largest d_max/L).
        vol = float(np.pi / 6.0 * np.sum(counts * p.diameters ** 3))
        L = (vol / beta) ** (1.0 / 3.0)
        dmol = float(p.diameters[-1]) / L
        ok = int(counts.min()) >= MIN_PHASE_COUNT
        flag = "" if ok else "   <-- REFUSED, minority phase too small"
        if dmol > 0.25:
            flag += "   <-- FINITE SIZE RISK d_max/L > 0.25"
        if not ok:
            refused.append(comp)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}"
            f"{counts[2]:6d}{counts[1]:6d}{counts[0]:6d}"
            f"{p.size_ratio:9.2f}{p.cv:7.4f}{p.d32:8.4f}{dmol:9.4f}  {y_string(comp)}{flag}")
    say(f"   compositions registered: {len(COMPOSITIONS)}; refused: {len(refused)} "
        f"{refused}")
    say()

    # --- the IN-SAMPLE residual, which is where the bars come from --------------------
    say("=== in-sample residuals on the BINARY calibration set (the bars' origin) ===")
    say("   A ternary must be predicted AS WELL AS the binaries the closure was")
    say("   calibrated on. Anything else is a bar invented for this block.")
    insample = {}
    for name, fn, setname, pset in (
            ("J", lambda d, y: mm.packing_density(d, y, beta, paJ, pbJ),
             "its OWN calibration set, ratios <= 3", pts23),
            ("J", lambda d, y: mm.packing_density(d, y, beta, paJ, pbJ),
             "all binaries incl. r = 4, 5 (outside its validity)", allpts),
            ("TAB", lambda d, y: mc.packing_density_tabulated(d, y, beta, fa, fb),
             "all binaries (it is fitted per ratio)", allpts)):
        res = []
        for r, x, phi in pset:
            dd = np.array([r, 1.0])
            yy = np.array([1.0 - x, x])
            res.append((fn(dd, yy) - phi) / SCALE)
        res = np.asarray(res)
        insample[(name, setname)] = res
        say(f"   {name:>4} on {setname}:")
        say(f"        n = {res.size}, rms {float(np.sqrt(np.mean(res**2))):.2f}, "
            f"max |res| {float(np.max(np.abs(res))):.2f}, mean signed {float(res.mean()):+.2f}, "
            f"{int(np.sum(res > 0))} of {res.size} positive")
    say("   The ternary ladder's pairs are u = 0.50 and u = 0.25. u = 0.25 is r = 4, the")
    say("   ratio X1 found J UNDER-predicts by +7.5 rms, so J is the CONTROL here and TAB")
    say("   is the primary candidate: TAB is the only one calibrated where this block sits.")
    say()

    # --- a NAMED non-additive alternative, so the test has two arms -------------------
    say("=== the NON-ADDITIVE alternative: SEQUENTIAL (nested-binary) filling ===")
    say("   The additive model sums an interaction over every pair. The oldest")
    say("   alternative in the powder literature is SEQUENTIAL: the two coarsest classes")
    say("   pack first, and the finest class then fills the voids of that COMPOSITE as")
    say("   if it were one medium. It uses the SAME measured binary closure, adds NO free")
    say("   parameter, and reduces EXACTLY to the binary closure at K = 2 -- so the")
    say("   binary calibration cannot tell the two apart, which is precisely why a")
    say("   ternary is the cheapest test that can.")
    say("   It needs ONE convention: what diameter the composite presents to the fines.")
    say("   That convention is a knob, so BOTH defensible choices are computed and the")
    say("   test is declared not informative if they bracket the additive prediction.")

    def sequential(comp, dmode: str, closure: str) -> float:
        """(y_L, y_M, y_S) -> phi under the nested construction."""
        y_l, y_m, y_s = comp
        p = make_ternary(comp)
        d = p.diameters                      # ascending: [d_S, d_M, d_L]
        def dens(dd, yy):
            if closure == "J":
                return mm.packing_density(dd, yy, beta, paJ, pbJ)
            return mc.packing_density_tabulated(dd, yy, beta, fa, fb)
        # step 1: the coarse sub-pack {M, L}
        w = y_m + y_l
        if w <= 0:
            return beta
        ysub = np.array([y_m / w, y_l / w])
        dsub = np.array([d[1], d[2]])
        g_sub = dens(dsub, ysub)
        # its effective diameter, by the chosen convention
        nsub = ysub / dsub**3
        nsub = nsub / nsub.sum()
        m2 = float(np.sum(nsub * dsub**2)); m3 = float(np.sum(nsub * dsub**3))
        m4 = float(np.sum(nsub * dsub**4))
        d_eff = (m3 / m2) if dmode == "d32" else (m4 / m3)
        # step 2: the fines against that composite, with the composite's OWN density in
        # place of beta. The binary form of the model with beta -> g_sub for the coarse
        # branch; the fine branch still packs at beta, because the fines are monodisperse.
        u = d[0] / d_eff
        if closure == "J":
            a_ = float(np.sqrt(1 - (1 - u) ** paJ)); b_ = float(1 - (1 - u) ** pbJ)
        else:
            a_ = float(fa(u)); b_ = float(fb(u))
        x = y_s
        coarse = g_sub / (1.0 - (1.0 - a_) * x)
        fine = beta / (1.0 - (1.0 - beta + b_ * beta * (1.0 - 1.0 / beta)) * (1.0 - x))
        return float(min(coarse, fine))

    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'ADD(TAB)':>10}{'SEQ d32':>9}{'SEQ d43':>9}"
        f"{'(SEQd32-ADD)':>14}{'(SEQd43-ADD)':>14}")
    seq_sep = []
    for comp in COMPOSITIONS:
        p = make_ternary(comp)
        add = mc.packing_density_tabulated(p.diameters, p.volume_fractions, beta, fa, fb)
        s32 = sequential(comp, "d32", "TAB")
        s43 = sequential(comp, "d43", "TAB")
        seq_sep.append(max(abs(s32 - add), abs(s43 - add)) / SCALE)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{add:10.4f}{s32:9.4f}{s43:9.4f}"
            f"{(s32 - add) / SCALE:+14.1f}{(s43 - add) / SCALE:+14.1f}")
    say(f"   ADDITIVE vs SEQUENTIAL separation, in scale: min {min(seq_sep):.1f}, "
        f"median {float(np.median(seq_sep)):.1f}, max {max(seq_sep):.1f}")
    say()

    # --- THE PREDICTIONS ---------------------------------------------------------------
    say("=== BLIND PREDICTIONS ===")
    say(f"   {'y_L':>6}{'y_M':>6}{'y_S':>6}{'J':>9}{'TAB':>9}{'ideal':>9}"
        f"{'J-beta':>9}{'TAB-beta':>10}{'(J-TAB)/scale':>15}")
    pred: dict[tuple, dict[str, float]] = {}
    seps = []
    for comp in list(COMPOSITIONS) + [MONO]:
        p = make_ternary(comp)
        d, y = p.diameters, p.volume_fractions
        pj = mm.packing_density(d, y, beta, paJ, pbJ)
        pt = mc.packing_density_tabulated(d, y, beta, fa, fb)
        # The zero-parameter ideal-filling limit for THIS composition: every interaction
        # off. Not a candidate -- a ceiling, printed so the two candidates can be read
        # against the most a filling picture could give.
        pi_ = mc.packing_density_ab(d, y, beta, 0.0, 0.0)
        pred[comp] = {"J": pj, "TAB": pt, "ideal": pi_}
        if comp != MONO:
            seps.append(abs(pj - pt) / SCALE)
        say(f"   {comp[0]:6.2f}{comp[1]:6.2f}{comp[2]:6.2f}{pj:9.4f}{pt:9.4f}{pi_:9.4f}"
            f"{(pj - beta) / SCALE:+9.1f}{(pt - beta) / SCALE:+10.1f}"
            f"{(pj - pt) / SCALE:+15.1f}")
    say()

    # --- candidate separation, the thing that decides whether the test is informative ---
    say("=== CANDIDATE SEPARATION ===")
    rms_sep = float(np.sqrt(np.mean(np.square(seps))))
    say(f"   |J - TAB| over the six compositions, in units of PHI_SCALE = {SCALE:.1e}:")
    say(f"      min {min(seps):.2f}   median {float(np.median(seps)):.2f}   "
        f"max {max(seps):.2f}   rms {rms_sep:.2f}")
    say(f"   The n = 3 seed sd on a binary cell of this block is ~{1.0:.1f}-{2.0:.1f} "
        f"scale, so the standard error of a three-seed mean is ~0.6-1.2 scale.")
    say("   REGISTERED RULE (same form as X2 in PREREG_MIXTURE_EXT): the comparison of")
    say("   J against TAB is declared NOT INFORMATIVE unless the rms residual ratio")
    say("   loser/winner reaches 1.6. That is a statement about the OUTCOME, and it is")
    say("   separate from the statement about the DESIGN, which is the number above:")
    say(f"   the two candidates' predictions differ by {rms_sep:.2f} scale in rms. A")
    say("   separation below ~2 scale cannot be resolved by three seeds, so T1/T2 below")
    say("   are graded against beta and against each other's SIGN, not only their rms.")
    say()

    # --- what each hypothesis implies, written before the data ------------------------
    say("=== what each outcome means (written before any cell ran) ===")
    say("   ADDITIVE:      the ternary residual is the size the binary residual is, and")
    say("                  its SIGN is not systematically positive. The closure's")
    say("                  trouble with the continuous families is then a calibration")
    say("                  problem at small u, not an additivity failure.")
    say("   NON-ADDITIVE:  the ternary is UNDER-PREDICTED, in the same direction as the")
    say("                  continuous families (measured phi above predicted phi), by")
    say("                  more than the bar. Three pairs acting together fill better")
    say("                  than the sum of their separate effects.")
    say("   The continuous-family under-prediction the block is chasing: TAB misses the")
    say("   lognormal by -6.9 scale at CV 0.5 (measured ABOVE predicted), held-out chi2")
    say("   108.1 against J's 34.9.")
    say()

    with open(os.path.join(OUT, "predict_ternary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    say(f"-> {os.path.join(OUT, 'predict_ternary.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

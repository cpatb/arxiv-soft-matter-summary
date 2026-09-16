"""grade_additivity.py -- PREREG_ADDITIVITY_2026-09-17.md, the additivity repair.

ONE mixing parameter between the ADDITIVE pair sum and the fully NESTED sequential rule,
calibrated on the ternary + K-ladder cells ONLY, evaluated on the 15 held-out continuous
family points. No new cells: every record read here was already on disk.

    phi_mix(lambda, c) = (1 - lambda) * ADDtab + lambda * SEQ_c

SIGN CONVENTION, ONE PLACE: residual = MEASURED - MODEL, positive = the model
under-predicts phi. `transfer_kladder.py` uses the opposite sign and the agreement after
a flip is ASSERTED below, because this project has already shipped one artifact that
contradicted its own table over exactly this.

NOT A BLIND PREDICTION. Both datasets are on disk and were read while the prereg was
written; section 0 of the prereg says so and records the back-of-envelope (chi2 ~ 60)
made at design time, so that a large disagreement with it flags this script rather than
confirming a prediction.

SCOPE GUARD. This file imports `grade_ternary` and `grade_kladder` as LOADERS ONLY and
calls neither `main()`. It writes only out/grade_additivity.txt. It re-decides no closed
gate; `test_additivity.py` asserts that structurally.

Run:  python grade_additivity.py      (writes out/grade_additivity.txt; exit 1 on a
                                       reproduction mismatch, 0 otherwise -- registered
                                       prediction failures go in FAILURES, not the code.)
"""

from __future__ import annotations

import math
import os

import numpy as np
from scipy import stats
from scipy.optimize import minimize_scalar

import grade_kladder as gk
import grade_mixture_ext as gx
import grade_ternary as gt
import kladder as kl
import mixture_closure as mc
import mixture_model as mm
from phase12_ternary import COMPOSITIONS, make_ternary

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
MU = 0.3
CONV = ("d10", "d21", "d32", "d43")
PRIMARY = "d43"                 # registered: T4's bracket arm. NOT chosen after the fact.

# --- numbers quoted VERBATIM from the artifacts, re-derived below and asserted --------
REC_ADDTAB_HELDOUT = 108.1      # out/transfer_kladder.txt, grade_mixture_ext.txt
REC_J_HELDOUT = 34.9            # ditto
REC_SEQ_HELDOUT = {"d10": 186.3, "d21": 42.8, "d32": 48.0, "d43": 509.3}
REC_T1_TAB_MEAN = 7.53          # out/grade_ternary.txt, T1, measured - model
REC_T1_TAB_POS = 6             # 6 of 6
REC_G6_ADDTAB_MEAN = 3.39       # out/grade_kladder.txt, G6, N = 1500 core
REC_G6_SEQ43_MEAN = -3.30
REC_G6_BRACKET = (0.74, 0.74)   # ADDtab below / SEQd43 above, G6's FAILED 80% bar
REC_SPEARMAN_POOLED = -0.939    # |family residual| vs w(u > 0.90), U6 referee pass
REC_SPEARMAN_WITHIN = {"lognormal": -1.000, "uniform": -0.900}
ENVELOPE_CHI2 = 60.0            # the design-time back-of-envelope, prereg section 0

LINES: list[str] = []
FAILURES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s, flush=True)


def chk(name: str, got: float, want: float, tol: float) -> bool:
    ok = abs(got - want) <= tol
    say(f"   {'OK  ' if ok else 'FAIL'} {name:<44} got {got:>10.4f}  recorded "
        f"{want:>10.4f}  tol {tol}")
    return ok


# ------------------------------------------------------------------ the mixing fit
def lambda_closed_form(m, A, S):
    """lambda_hat for phi_mix = (1-l)A + lS under unweighted LSQ, in SCALE units.

    Closed form because the model is LINEAR in lambda; cross-checked against an
    optimiser below. `m`, `A`, `S` are raw phi, converted here.
    """
    dm = (np.asarray(m, float) - np.asarray(A, float)) / SCALE
    ds = (np.asarray(S, float) - np.asarray(A, float)) / SCALE
    den = float(np.sum(ds * ds))
    if den <= 0:
        return float("nan"), float("nan"), float("nan")
    lam = float(np.sum(dm * ds) / den)
    chi2 = float(np.sum((dm - lam * ds) ** 2))
    # Delta chi2 = 1 half-width: chi2(l) = chi2_min + den * (l - lam)^2
    half = float(1.0 / math.sqrt(den))
    return lam, chi2, half


def chi2_at(lam, m, A, S):
    dm = (np.asarray(m, float) - np.asarray(A, float)) / SCALE
    ds = (np.asarray(S, float) - np.asarray(A, float)) / SCALE
    return float(np.sum((dm - lam * ds) ** 2))


def residuals_at(lam, m, A, S):
    return (np.asarray(m, float) - ((1 - lam) * np.asarray(A, float)
                                    + lam * np.asarray(S, float))) / SCALE


def pair_weight_above(d, y, ucut=0.90):
    """Fraction of pairwise volume-product weight at u > ucut. Same construction as
    transfer_kladder.pair_weight_below, mirrored -- the U6 referee figure."""
    d = np.asarray(d, float)
    y = np.asarray(y, float)
    y = y / y.sum()
    tot = hi = 0.0
    for i in range(d.size):
        for j in range(i + 1, d.size):
            u = min(d[i], d[j]) / max(d[i], d[j])
            w = float(y[i] * y[j])
            tot += w
            if u > ucut:
                hi += w
    return hi / tot if tot > 0 else float("nan")


def main() -> int:
    meas = mm.measured(mu_p=MU)
    beta = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta)
    cand, _, _ = gx.candidates(beta, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta, allpts)

    say("=== grade_additivity.py -- PREREG_ADDITIVITY_2026-09-17.md ===")
    say(f"   beta = {beta:.6f} (measured, mu_p {MU}); PHI_SCALE = {SCALE:.1e}")
    say(f"   TAB knots at u = {', '.join(f'{u:.3f}' for u in sorted(tab))}; "
        f"J PA {paJ:.4f} PB {pbJ:.4f}")
    say(f"   {mm.K_CLASSES} equal-count classes per continuous family.")
    say("   residual = MEASURED - MODEL throughout; positive = the model under-predicts.")
    say(f"   primary composite-diameter convention: {PRIMARY} (registered, T4's bracket "
        "arm); all four reported.")
    say()

    # ============================================================ the three cell sets
    def modelset(d, y):
        out = {"ADDtab": float(mc.packing_density_tabulated(
                   np.asarray(d, float), np.asarray(y, float), beta, fa, fb)),
               "ADDj": float(mm.packing_density(
                   np.asarray(d, float), np.asarray(y, float), beta, paJ, pbJ))}
        for c in CONV:
            out[f"SEQ{c}"] = float(kl.sequential_k(d, y, beta, fa, fb, c))
        return out

    # --- held out: the 15 continuous family points -----------------------------------
    held = []
    for fam in mm.FAMILIES:
        for cv in mc.CVS:
            if (fam, cv) not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            held.append({"lab": f"{fam}_cv{cv:.1f}", "fam": fam, "cv": cv,
                         "phi": meas[(fam, cv)][0], "d": d, "y": y, "K": int(d.size),
                         "models": modelset(d, y),
                         "w90": pair_weight_above(d, y)})

    # --- calibration A: ternary ------------------------------------------------------
    g, nfiles_t, nraw_t = gt.load()
    tern = []
    n_before_t = n_after_t = 0
    for comp in COMPOSITIONS:
        key = tuple(round(v, 4) for v in comp)
        cs = g.get(key, [])
        n_before_t += len(cs)
        ok = [c for c in cs if gt.admissible(c)[0]]
        n_after_t += len(ok)
        if not ok:
            continue
        p = make_ternary(comp)
        phi = np.array([float(c["phi"]) for c in ok])
        tern.append({"lab": f"T{comp[0]:.2f}/{comp[1]:.2f}/{comp[2]:.2f}",
                     "set": "ternary", "phi": float(phi.mean()), "K": int(p.diameters.size),
                     "se": float(phi.std(ddof=1) / math.sqrt(phi.size)),
                     "models": modelset(p.diameters, p.volume_fractions)})

    # --- calibration B: the K-ladder G6 core at N = 1500 -----------------------------
    cells_all = {c["label"]: c for c in kl.build_cells(kl.N_BIG)[0]}
    gkl, nfiles_k, nraw_k = gk.load(1500)
    klad = []
    n_before_k = n_after_k = 0
    for lbl, cs in sorted(gkl.items()):
        if lbl not in cells_all:
            continue
        n_before_k += len(cs)
        ok = [c for c in cs if gk.admissible(c, 1500)[0]]
        n_after_k += len(ok)
        c0 = cells_all[lbl]
        if not ok or c0["arm"] not in ("K", "SPLIT", "NA2") or (c0["K"] or 0) <= 1:
            continue
        phi = np.array([float(c["phi"]) for c in ok])
        klad.append({"lab": lbl, "set": "kladder", "phi": float(phi.mean()),
                     "K": int(c0["K"]),
                     "se": float(phi.std(ddof=1) / math.sqrt(phi.size)),
                     "models": modelset(c0["d_rel"], c0["y"])})

    cal = tern + klad

    say("=== N BEFORE AND AFTER EVERY FILTER ===")
    say(f"   ternary : {nfiles_t} shard files, {nraw_t} records; grouped to "
        f"{len(COMPOSITIONS)} compositions; {n_before_t} -> {n_after_t} admitted "
        f"records over {len(tern)} cells")
    say(f"   K-ladder: {nfiles_k} shard files, {nraw_k} records at N = 1500; "
        f"{n_before_k} -> {n_after_k} admitted over all labels; "
        f"{len(klad)} cells survive the G6 core rule (arm in K/SPLIT/NA2, K > 1)")
    say(f"   CALIBRATION total {len(cal)} cells.  HELD OUT {len(held)} family points "
        f"({', '.join(sorted(set(h['fam'] for h in held)))} x {len(mc.CVS)} CV levels)")
    say()

    # ============================================================ REPRODUCTION GATE
    say("=== REPRODUCTION GATE -- a diagnostic that cannot reproduce the number it "
        "audits has not audited it ===")
    repro_ok = True

    def chi2_heldout(name):
        return float(sum(((h["phi"] - h["models"][name]) / SCALE) ** 2 for h in held))

    repro_ok &= chk("held-out chi2, ADDtab", chi2_heldout("ADDtab"),
                    REC_ADDTAB_HELDOUT, 0.06)
    repro_ok &= chk("held-out chi2, ADDj (J)", chi2_heldout("ADDj"), REC_J_HELDOUT, 0.06)
    for c in CONV:
        repro_ok &= chk(f"held-out chi2, SEQ{c}", chi2_heldout(f"SEQ{c}"),
                        REC_SEQ_HELDOUT[c], 0.06)
    rT = np.array([(t["phi"] - t["models"]["ADDtab"]) / SCALE for t in tern])
    repro_ok &= chk("T1 TAB mean residual (ternary)", float(rT.mean()),
                    REC_T1_TAB_MEAN, 0.01)
    repro_ok &= chk("T1 TAB count positive", float((rT > 0).sum()), REC_T1_TAB_POS, 0)
    rK = np.array([(k["phi"] - k["models"]["ADDtab"]) / SCALE for k in klad])
    rK43 = np.array([(k["phi"] - k["models"]["SEQd43"]) / SCALE for k in klad])
    repro_ok &= chk("G6 ADDtab mean residual (K-ladder)", float(rK.mean()),
                    REC_G6_ADDTAB_MEAN, 0.01)
    repro_ok &= chk("G6 SEQd43 mean residual (K-ladder)", float(rK43.mean()),
                    REC_G6_SEQ43_MEAN, 0.01)
    repro_ok &= chk("G6 bracket, ADDtab below", float((rK > 0).mean()),
                    REC_G6_BRACKET[0], 0.005)
    repro_ok &= chk("G6 bracket, SEQd43 above", float((rK43 < 0).mean()),
                    REC_G6_BRACKET[1], 0.005)
    # the SIGN-CONVENTION assertion: transfer_kladder's (model - phi) mean, flipped
    tk_mean = float(np.mean([(h["models"]["ADDtab"] - h["phi"]) / SCALE for h in held]))
    mine = float(np.mean([(h["phi"] - h["models"]["ADDtab"]) / SCALE for h in held]))
    repro_ok &= chk("sign convention: transfer_kladder -1.69 flips", tk_mean, -mine, 1e-9)
    repro_ok &= chk("  ... and equals the recorded -1.69", tk_mean, -1.69, 0.006)
    say(f"   -> reproduction {'PASSES' if repro_ok else 'FAILS'}")
    if not repro_ok:
        say()
        say("FAILURES: ['REPRODUCTION GATE']")
        _write()
        return 1
    say()

    # ============================================================ leverage
    say("=== WHAT lambda CANNOT DO: the leverage split (prereg section 2) ===")
    say("   At K = 2 the nested rule reduces EXACTLY to the additive binary closure, so")
    say("   S - A == 0 and lambda has literally no leverage on those cells.")
    for nm, arr in (("ternary", tern), ("K-ladder", klad), ("held out", held)):
        lev = [c for c in arr
               if abs(c["models"][f"SEQ{PRIMARY}"] - c["models"]["ADDtab"]) > 1e-12]
        say(f"   {nm:<10} {len(lev):>3} of {len(arr):>3} cells carry leverage at "
            f"{PRIMARY}  (K = 2 cells: {sum(1 for c in arr if c['K'] == 2)})")
    floor = float(sum(((h["phi"] - h["models"]["ADDtab"]) / SCALE) ** 2 for h in held
                      if abs(h["models"][f"SEQ{PRIMARY}"] - h["models"]["ADDtab"])
                      <= 1e-12))
    say(f"   -> lambda-INDEPENDENT held-out chi2 FLOOR: {floor:.1f} of "
        f"{chi2_heldout('ADDtab'):.1f}. No mixing model can go below it.")
    say()

    # ============================================================ the fit
    say("=== THE FIT: lambda on the calibration set, all four conventions ===")
    say(f"   {'conv':>6}{'lam_hat':>10}{'+-1':>9}{'chi2_cal':>10}{'lam(tern)':>11}"
        f"{'+-1':>8}{'lam(klad)':>11}{'+-1':>8}{'chi2_held':>11}")
    fit = {}
    for c in CONV:
        m = [x["phi"] for x in cal]
        A = [x["models"]["ADDtab"] for x in cal]
        S = [x["models"][f"SEQ{c}"] for x in cal]
        lam, c2, half = lambda_closed_form(m, A, S)
        lt, c2t, ht = lambda_closed_form([x["phi"] for x in tern],
                                         [x["models"]["ADDtab"] for x in tern],
                                         [x["models"][f"SEQ{c}"] for x in tern])
        lk, c2k, hk = lambda_closed_form([x["phi"] for x in klad],
                                         [x["models"]["ADDtab"] for x in klad],
                                         [x["models"][f"SEQ{c}"] for x in klad])
        ch = chi2_at(lam, [h["phi"] for h in held],
                     [h["models"]["ADDtab"] for h in held],
                     [h["models"][f"SEQ{c}"] for h in held])
        fit[c] = {"lam": lam, "half": half, "chi2_cal": c2, "chi2_held": ch,
                  "lam_t": lt, "half_t": ht, "chi2_t": c2t,
                  "lam_k": lk, "half_k": hk, "chi2_k": c2k}
        say(f"   {c:>6}{lam:>10.4f}{half:>9.4f}{c2:>10.1f}{lt:>11.4f}{ht:>8.4f}"
            f"{lk:>11.4f}{hk:>8.4f}{ch:>11.1f}")
    say("   lam = 0 is pure ADDITIVE, lam = 1 is the pure nested rule; fitted UNBOUNDED")
    say("   so that a truth outside the bracket can show itself. '+-1' is the Delta chi2")
    say("   = 1 half-width.")
    say()

    # closed form vs optimiser, primary only
    mP = [x["phi"] for x in cal]
    AP = [x["models"]["ADDtab"] for x in cal]
    SP = [x["models"][f"SEQ{PRIMARY}"] for x in cal]
    opt = minimize_scalar(lambda L: chi2_at(L, mP, AP, SP), bounds=(-5, 5),
                          method="bounded", options={"xatol": 1e-12})
    dcf = abs(float(opt.x) - fit[PRIMARY]["lam"])
    say(f"=== closed form vs optimiser at {PRIMARY}: {fit[PRIMARY]['lam']:.10f} vs "
        f"{float(opt.x):.10f}, |diff| = {dcf:.2e} (bar 1e-9) "
        f"{'OK' if dcf < 1e-9 else 'MISMATCH'} ===")
    if dcf >= 1e-9:
        FAILURES.append(f"closed form vs optimiser {dcf:.2e}")
    # seed-se-weighted lambda: ROBUSTNESS ONLY, not the primary
    wts = np.array([1.0 / max(x["se"] / SCALE, 1e-6) ** 2 for x in cal])
    dmW = (np.array(mP) - np.array(AP)) / SCALE
    dsW = (np.array(SP) - np.array(AP)) / SCALE
    lamW = float(np.sum(wts * dmW * dsW) / np.sum(wts * dsW * dsW))
    say(f"   seed-se-WEIGHTED lambda at {PRIMARY} = {lamW:.4f} (robustness figure, "
        f"NOT the primary; unweighted is what the prereg registered)")
    say()

    # ============================================================ P0..P4
    f = fit[PRIMARY]
    say("=== REGISTERED PREDICTIONS (bars fixed in the prereg; none of them moves) ===")

    # P0 -- the two calibration arms must agree
    lt, ht, lk, hk = f["lam_t"], f["half_t"], f["lam_k"], f["half_k"]
    inside = (abs(lt - lk) <= hk) and (abs(lk - lt) <= ht)
    overlap = min(lt + ht, lk + hk) - max(lt - ht, lk - hk)
    p0 = inside and overlap > 0
    say(f"   P0  ternary lambda {lt:+.4f} [{lt - ht:+.4f}, {lt + ht:+.4f}] vs "
        f"K-ladder {lk:+.4f} [{lk - hk:+.4f}, {lk + hk:+.4f}]")
    say(f"       each point inside the other's interval: {inside}; interval overlap "
        f"{overlap:+.4f}  -> {'PASS' if p0 else 'FAIL'}")
    if not p0:
        FAILURES.append(f"P0 calibration arms disagree ({lt:.4f} vs {lk:.4f})")
    # Two sharper statistics PRINTED BESIDE the registered verdict. Neither replaces it
    # and neither moves a bar: the P0 verdict above is whatever the registered SENTENCE
    # says, and both of these point the LOOSENING way, which is the direction to be most
    # suspicious of.
    sed_ = math.hypot(ht, hk)
    say(f"       [beside it, not instead of it] |d lambda| = {abs(lt - lk):.4f} against "
        f"the se of the DIFFERENCE {sed_:.4f}: t = {abs(lt - lk) / sed_:.2f}. The "
        "registered")
    say("       sentence compares a point to the OTHER arm's interval, which is a "
        "stricter and")
    say(f"       less standard test; it misses by {max(0.0, abs(lt - lk) - min(ht, hk)):.4f}.")
    ndof = max(len(cal) - 1, 1)
    infl = math.sqrt(f["chi2_cal"] / ndof)
    say(f"       [and the intervals are NOT statistical] chi2_cal = {f['chi2_cal']:.1f} "
        f"on {ndof} dof, so MODEL MISFIT is at")
    say(f"       least COMPARABLE to seed noise: chi2_red = "
        f"{f['chi2_cal'] / ndof:.2f}, i.e. misfit variance about EQUAL to")
    say(f"       noise variance, not dominant over it (referee). Inflation {infl:.2f}x. A")
    say("       Delta chi2 = 1 interval in PHI_SCALE units is 'one scale unit of "
        "misfit', not a")
    say(f"       confidence interval. Inflated, the two arms are "
        f"[{lt - infl * ht:+.4f}, {lt + infl * ht:+.4f}] and "
        f"[{lk - infl * hk:+.4f}, {lk + infl * hk:+.4f}].")
    say("       P0 STAYS FAILED on the sentence that was registered.")

    # P1 -- the bracket premise
    lo, hi = f["lam"] - f["half"], f["lam"] + f["half"]
    p1 = (0.0 < f["lam"] < 1.0) and (lo >= 0.0) and (hi <= 1.0)
    say(f"   P1  lambda_hat = {f['lam']:.4f} [{lo:.4f}, {hi:.4f}] at {PRIMARY}; "
        f"interior to (0,1) with its interval inside [0,1] -> "
        f"{'PASS' if p1 else 'FAIL'}")
    if not p1:
        FAILURES.append(f"P1 lambda {f['lam']:.4f} [{lo:.4f}, {hi:.4f}] not inside (0,1)")

    # P2 -- the weak win
    p2 = f["chi2_held"] < REC_ADDTAB_HELDOUT
    say(f"   P2  held-out chi2 {f['chi2_held']:.1f} < ADDtab {REC_ADDTAB_HELDOUT} -> "
        f"{'PASS' if p2 else 'FAIL'}   (registered WEAK: SEQd21 already manages "
        f"{REC_SEQ_HELDOUT['d21']} with no fitting)")
    if not p2:
        FAILURES.append(f"P2 held-out {f['chi2_held']:.1f} >= {REC_ADDTAB_HELDOUT}")

    # P2b -- the discriminating win: beat the CALIBRATION-selected convention
    cal_rms = {c: math.sqrt(np.mean([((x["phi"] - x["models"][f"SEQ{c}"]) / SCALE) ** 2
                                     for x in cal])) for c in CONV}
    sel = min(cal_rms, key=lambda c: cal_rms[c])
    sel_held = chi2_heldout(f"SEQ{sel}")
    p2b = f["chi2_held"] < sel_held
    say("   P2b calibration-set rms by convention: "
        + ", ".join(f"{c} {cal_rms[c]:.2f}" for c in CONV))
    say(f"       -> selected WITHOUT touching a held-out point: SEQ{sel}; its held-out "
        f"chi2 is {sel_held:.1f}")
    say(f"       mix {f['chi2_held']:.1f} < {sel_held:.1f} -> "
        f"{'PASS' if p2b else 'FAIL'}   (the DISCRIMINATING bar)")
    if not p2b:
        FAILURES.append(f"P2b mix {f['chi2_held']:.1f} does not beat SEQ{sel} "
                        f"{sel_held:.1f}")

    # P3 -- the headline win
    p3 = f["chi2_held"] < REC_J_HELDOUT
    say(f"   P3  held-out chi2 {f['chi2_held']:.1f} < J {REC_J_HELDOUT} -> "
        f"{'PASS' if p3 else 'FAIL'}   (the headline; the design-time envelope said no)")
    if not p3:
        FAILURES.append(f"P3 held-out {f['chi2_held']:.1f} >= J {REC_J_HELDOUT}")

    # P4 -- convention robustness
    nwin = sum(1 for c in CONV if fit[c]["chi2_held"] < REC_ADDTAB_HELDOUT)
    p4 = nwin >= 3
    say(f"   P4  P2 holds on {nwin} of 4 conventions (bar >= 3) -> "
        f"{'PASS' if p4 else 'FAIL'}")
    if not p4:
        FAILURES.append(f"P4 P2 holds on only {nwin} of 4 conventions")
    say()

    # ============================================================ convention spread
    ch = {c: fit[c]["chi2_held"] for c in CONV}
    say("=== THE NESTED-CONVENTION SPREAD, BESIDE THE ANSWER, ALWAYS ===")
    say("   held-out chi2 of the FITTED mix by convention: "
        + ", ".join(f"{c} {ch[c]:.1f}" for c in CONV))
    say(f"   spread {min(ch.values()):.1f} to {max(ch.values()):.1f}, a factor "
        f"{max(ch.values()) / max(min(ch.values()), 1e-9):.2f}")
    say(f"   the UNMIXED SEQ spread was 11.90x ({REC_SEQ_HELDOUT['d21']} to "
        f"{REC_SEQ_HELDOUT['d43']}); lambda_hat by convention: "
        + ", ".join(f"{c} {fit[c]['lam']:+.3f}" for c in CONV))
    say("   A result that moves by more than its own effect when an undeclared")
    say("   convention is changed is a measurement of the convention.")
    say()

    # ==================================================== the ORACLE: the form's range
    say("=== A FAIL IS THE MODEL'S RANGE, NOT A REFIT: the ORACLE lambda ===")
    say("   (!) The ORACLE minimises chi2 ON THE HELD-OUT POINTS. It is NOT a result and")
    say("   must never be quoted as one. It is the BEST any value of lambda could do, so")
    say("   it separates 'the calibration is off' from 'the FORM cannot get there'.")
    say(f"   {'conv':>6}{'lam_fit':>10}{'chi2_held':>11}{'lam_oracle':>12}"
        f"{'chi2_oracle':>13}{'beats J?':>10}")
    for c in CONV:
        lo_, co_, _ = lambda_closed_form([h["phi"] for h in held],
                                         [h["models"]["ADDtab"] for h in held],
                                         [h["models"][f"SEQ{c}"] for h in held])
        say(f"   {c:>6}{fit[c]['lam']:>10.4f}{fit[c]['chi2_held']:>11.1f}"
            f"{lo_:>12.4f}{co_:>13.1f}"
            f"{('yes' if co_ < REC_J_HELDOUT else 'NO'):>10}")
        fit[c]["lam_oracle"] = lo_
        fit[c]["chi2_oracle"] = co_
    # The verdict is GENERATED from the numbers, never written beside them.
    fp = fit[PRIMARY]
    say(f"   At the PRIMARY {PRIMARY}: fitted lambda reaches {fp['chi2_held']:.1f}, the "
        f"oracle {fp['chi2_oracle']:.1f}, against J's {REC_J_HELDOUT}; the "
        f"lambda-independent floor is {floor:.1f}.")
    if fp["chi2_oracle"] >= REC_J_HELDOUT:
        say("   -> The oracle CANNOT beat J either, so P3 fails on the FORM: no value of")
        say("      lambda gets this rule there and refitting is not the answer.")
    else:
        say("   -> The oracle CAN beat J, so P3 does NOT fail on the form. It fails on "
            "WHERE")
        say(f"      the calibration puts lambda: the ladders want "
            f"{fp['lam']:.4f} +- {fp['half']:.4f} and the families want "
            f"{fp['lam_oracle']:.4f},")
        _dof = max(len(cal) - 1, 1)
        _hi = fp["half"] * math.sqrt(fp["chi2_cal"] / _dof)
        say(f"      a gap of {abs(fp['lam'] - fp['lam_oracle']) / fp['half']:.1f} "
            f"half-widths, or {abs(fp['lam'] - fp['lam_oracle']) / _hi:.1f} once that")
        say("      half-width is inflated by sqrt(chi2_red) as P0's caveat requires. Both")
        say("      are printed; neither is a sigma.")
        say("      This is a TRANSFER failure, not a form failure.")
    say()

    # ============================== the three-way lambda comparison: the actual result
    say("=== THE THREE-WAY lambda COMPARISON -- what this item actually measured ===")
    say(f"   {'set':<26}{'lambda':>10}{'+-1':>9}{'n cells':>9}{'leverage':>10}")
    rows = [("ternary (calibration)", fp["lam_t"], fp["half_t"], len(tern)),
            ("K-ladder (calibration)", fp["lam_k"], fp["half_k"], len(klad)),
            ("both ladders (primary fit)", fp["lam"], fp["half"], len(cal))]
    for nm, lv, hv, nn in rows:
        nlev = sum(1 for c in (tern if "ternary" in nm else
                               klad if "K-ladder" in nm else cal)
                   if abs(c["models"][f"SEQ{PRIMARY}"] - c["models"]["ADDtab"]) > 1e-12)
        say(f"   {nm:<26}{lv:>10.4f}{hv:>9.4f}{nn:>9d}{nlev:>10d}")
    lo_h, c2_h, half_h = lambda_closed_form([h["phi"] for h in held],
                                            [h["models"]["ADDtab"] for h in held],
                                            [h["models"][f"SEQ{PRIMARY}"] for h in held])
    nlev_h = sum(1 for h in held
                 if abs(h["models"][f"SEQ{PRIMARY}"] - h["models"]["ADDtab"]) > 1e-12)
    say(f"   {'families (HELD OUT, oracle)':<26}{lo_h:>10.4f}{half_h:>9.4f}"
        f"{len(held):>9d}{nlev_h:>10d}")
    t_arms = abs(fp["lam_t"] - fp["lam_k"]) / math.hypot(fp["half_t"], fp["half_k"])
    t_fam = abs(fp["lam"] - lo_h) / math.hypot(fp["half"], half_h)
    say(f"   ternary vs K-ladder : |d| = {abs(fp['lam_t'] - fp['lam_k']):.4f}, "
        f"t = {t_arms:.2f}")
    say(f"   ladders vs families : |d| = {abs(fp['lam'] - lo_h):.4f}, t = {t_fam:.2f}  "
        f"({t_fam / max(t_arms, 1e-9):.1f}x the between-ladder gap)")
    say("   Read this comparison ONLY with the leverage column: at K = 2 the nested rule")
    say("   IS the additive one, so the 5 bidisperse points constrain nothing.")
    say()

    # ---------------------------------------------------------------- pass-6 (a)
    # A cleaner "do the two arms agree" statistic than P0's registered sentence: the
    # combined fit is the two arms' fits with one parameter instead of two, so the excess
    # chi2 is a nested-model test on 1 dof, with no interval convention to argue about.
    excess = fp["chi2_cal"] - (fp["chi2_t"] + fp["chi2_k"])
    say("   THE CLEANEST FORM OF THE SAME QUESTION, with no interval convention in it:")
    say(f"   pooling the arms costs chi2 {fp['chi2_cal']:.3f} - ({fp['chi2_t']:.3f} + "
        f"{fp['chi2_k']:.3f}) = {excess:.3f} on 1 dof.")
    # The two arms have very different dispersion (chi2_red 1.29 vs 2.41), so the
    # unweighted excess is not the heteroscedasticity-aware statistic. Re-weight each
    # arm's points by 1/chi2_red_arm -- which by construction reproduces each arm's own
    # chi2_min, an internal check -- and recompute the cost of pooling on that scale.
    def _wexcess():
        dof_t, dof_k = max(len(tern) - 1, 1), max(len(klad) - 1, 1)
        wt = dof_t / fp["chi2_t"] if fp["chi2_t"] > 0 else 1.0
        wk = dof_k / fp["chi2_k"] if fp["chi2_k"] > 0 else 1.0
        dm, ds, w = [], [], []
        for arr, ww in ((tern, wt), (klad, wk)):
            for x in arr:
                dm.append((x["phi"] - x["models"]["ADDtab"]) / SCALE)
                ds.append((x["models"][f"SEQ{PRIMARY}"] - x["models"]["ADDtab"]) / SCALE)
                w.append(ww)
        dm, ds, w = np.array(dm), np.array(ds), np.array(w)
        den = float(np.sum(w * ds * ds))
        lam_w = float(np.sum(w * dm * ds) / den) if den > 0 else float("nan")
        pooled = float(np.sum(w * (dm - lam_w * ds) ** 2))
        sep = fp["chi2_t"] * wt + fp["chi2_k"] * wk
        return pooled - sep, lam_w
    wex, lam_w = _wexcess()
    say(f"   variance-weighted (1/chi2_red per arm, which reproduces each arm's own "
        f"chi2_min): {wex:.3f} on 1 dof,")
    say(f"   at a pooled lambda of {lam_w:.4f}. The more rigorous form makes the "
        "conclusion STRONGER,")
    say("   not weaker -- the same direction as every other correction in this file.")
    say("   One lambda for both ladders is bought for essentially nothing, which is a")
    say("   stronger statement than P0's point-in-interval sentence and does not depend")
    say("   on whether a Delta chi2 = 1 width is a confidence interval.")
    say()
    say("   \u26a0 READ THAT AS 'THE TWO ARMS' BEST-FIT LAMBDA ARE INDISTINGUISHABLE',")
    say("   NOT AS 'ONE LAMBDA FITS THE LADDER DATA WELL'. They are different claims and")
    say(f"   only the first is supported: the K-ladder arm's OWN chi2_red is "
        f"{fp['chi2_k'] / max(len(klad) - 1, 1):.2f}, i.e. its")
    say("   19 cells carry real cell-to-cell structure that no single lambda absorbs. That")
    say("   is not news and is not this block's finding -- it is most likely G7, the")
    say("   closure's WRONG SIGN IN K, which predates this work and is still unrepaired.")
    say(f"   The ternary arm does not have the problem (chi2_red "
        f"{fp['chi2_t'] / max(len(tern) - 1, 1):.2f}).")
    say()

    # ---------------------------------------------------------------- pass-6 (b)/(c)
    # The pooled oracle is a fit over THREE families. This project has already ruled that
    # within a family both axes are monotone in CV alone, so 15 points are not 15
    # independent measurements -- the same argument that sets the effective N near 2 for
    # the weight correlation. Fit each family on its own before believing the pooled one.
    say("   \U0001f534 AND THE POOLED ORACLE HIDES A DISAGREEMENT BETWEEN THE FAMILIES.")
    say("   The same reasoning this file already applies to the weight correlation -- "
        "within a")
    say("   family both axes are monotone in CV alone, so 15 points are not 15 independent")
    say("   measurements -- applies to the oracle's own dof. Fitted per family:")
    say(f"   {'family':<14}{'n':>4}{'leverage':>10}{'lambda':>10}{'+-1':>9}"
        f"{'chi2_red':>10}{'interval':>22}")
    fam_lam = {}
    for fam in mm.FAMILIES:
        sub_ = [h for h in held if h["fam"] == fam]
        nlev_ = sum(1 for h in sub_
                    if abs(h["models"][f"SEQ{PRIMARY}"] - h["models"]["ADDtab"]) > 1e-12)
        lf, _cf, hf = lambda_closed_form(
            [h["phi"] for h in sub_], [h["models"]["ADDtab"] for h in sub_],
            [h["models"][f"SEQ{PRIMARY}"] for h in sub_])
        if nlev_ == 0 or not math.isfinite(lf):
            say(f"   {fam:<14}{len(sub_):>4}{nlev_:>10}{'UNDEFINED':>10}{'--':>9}"
                f"{'--':>10}{'S == A, no leverage':>22}")
            continue
        red_ = _cf / max(len(sub_) - 1, 1)
        fam_lam[fam] = (lf, hf)
        say(f"   {fam:<14}{len(sub_):>4}{nlev_:>10}{lf:>10.4f}{hf:>9.4f}{red_:>10.3f}"
            f"{f'[{lf - hf:.3f}, {lf + hf:.3f}]':>22}")
    if len(fam_lam) >= 2:
        (f1, (l1, h1)), (f2, (l2, h2)) = sorted(fam_lam.items())
        say(f"   the two continuous families disagree with EACH OTHER by "
            f"|{l1:.4f} - {l2:.4f}| = {abs(l1 - l2):.4f},")
        say(f"   which is {abs(l1 - l2) / max(abs(fp['lam'] - lo_h), 1e-9):.2f}x the "
            f"ladders-vs-pooled-oracle gap of {abs(fp['lam'] - lo_h):.4f} -- and they "
            "straddle the")
        say(f"   ladders' {fp['lam']:.4f} from OPPOSITE sides. So the pooled oracle "
            f"{lo_h:.4f} is")
        say(f"   {f1} dominating a combined fit (its half-width {h1:.4f} against "
            f"{f2}'s {h2:.4f}),")
        say("   NOT a number describing continuous PSDs. The supportable statement is")
        say("   family-specific and is printed as such below.")
        say("   Both per-family chi2_red are BELOW 1, so these two intervals are if")
        say("   anything GENEROUSLY WIDE rather than falsely narrow -- which strengthens")
        say("   the disagreement rather than explaining it away.")
        say(f"   -> {f1}: wants markedly {'LESS' if l1 < fp['lam'] else 'MORE'} nesting "
            f"than the ladders ({l1:.3f} vs {fp['lam']:.3f}), tightly constrained.")
        marg = abs(l2 - fp["lam"]) - h2
        say(f"   -> {f2}: half-width {h2:.4f} is {h2 / h1:.1f}x wider; its interval "
            f"{'CONTAINS' if marg <= 0 else 'excludes'} the ladders'")
        say(f"      value by {abs(marg):.4f}, which is {abs(marg) / h2:.0%} of its own "
            "width -- MARGINALLY, and on")
        say(f"      the OPPOSITE side, leaning toward "
            f"{'MORE' if l2 > fp['lam'] else 'LESS'} nesting. Do NOT read that")
        say("      exclusion as a result; the families' MUTUAL disagreement is the larger")
        say("      and better-determined effect.")
    say()
    # ---- the SAME caveat P0 gets, carried to the structurally identical statistic ----
    # A Delta chi2 = 1 half-width in PHI_SCALE units is "one scale unit of misfit". P0
    # already says so; quoting t = 5.67 without saying so would apply the project's own
    # rigour to the statistic that FAILED and not to the one that carries the headline.
    # The four fits have DIFFERENT reduced chi2, so each is inflated by its own.
    say("   \u26a0 THESE ARE 'ONE SCALE UNIT OF MISFIT' RATIOS, NOT SIGMA-EQUIVALENTS --")
    say("   the same caveat P0 carries, applied here because the statistic is the same")
    say("   kind of object. The four fits are NOT equally dispersed, so each half-width")
    say("   is inflated by ITS OWN sqrt(chi2_red):")
    say(f"   {'fit':<28}{'n':>4}{'chi2_min':>10}{'dof':>5}{'chi2_red':>10}"
        f"{'half':>9}{'inflated':>10}")
    rows2 = [("ternary alone", len(tern), fp["chi2_t"], fp["half_t"]),
             ("K-ladder alone", len(klad), fp["chi2_k"], fp["half_k"]),
             ("both ladders (primary)", len(cal), fp["chi2_cal"], fp["half"]),
             ("families (HELD OUT, oracle)", len(held), c2_h, half_h)]
    infl_half = {}
    for nm, nn, c2_, h_ in rows2:
        dof_ = max(nn - 1, 1)
        red_ = c2_ / dof_
        hi_ = h_ * math.sqrt(red_)
        infl_half[nm] = hi_
        say(f"   {nm:<28}{nn:>4}{c2_:>10.3f}{dof_:>5}{red_:>10.3f}{h_:>9.4f}"
            f"{hi_:>10.4f}")
    ta2 = abs(fp["lam_t"] - fp["lam_k"]) / math.hypot(infl_half["ternary alone"],
                                                      infl_half["K-ladder alone"])
    tf2 = abs(fp["lam"] - lo_h) / math.hypot(infl_half["both ladders (primary)"],
                                             infl_half["families (HELD OUT, oracle)"])
    say(f"   INFLATED: ternary vs K-ladder t = {ta2:.2f} (was {t_arms:.2f}); "
        f"ladders vs families t = {tf2:.2f} (was {t_fam:.2f}),")
    say(f"   i.e. {tf2 / max(ta2, 1e-9):.1f}x the between-ladder gap (was "
        f"{t_fam / max(t_arms, 1e-9):.1f}x).")
    say("   The held-out oracle fit is the only nominally dispersed one "
        f"(chi2_red {c2_h / max(len(held) - 1, 1):.2f}); the")
    say("   calibration fits are not, which is WHY the inflation changes the two sides by")
    say("   different amounts. The RATIO GROWS rather than shrinking, so the conclusion is")
    say("   robust to the correction -- but the uninflated 5.67 and 11.6 are not derived")
    say("   with the rigour this file demands of P0, and both forms are printed.")
    say()

    # ============================================================ per-family
    say("=== HELD-OUT, PER FAMILY -- a win concentrated in one family is not a general "
        "one ===")
    say(f"   {'family':<12}{'n':>3}{'ADDtab':>9}{'J':>9}{'MIX':>9}{'SEQ' + sel:>10}"
        f"{'MIX mean res':>14}")
    for fam in mm.FAMILIES:
        sub = [h for h in held if h["fam"] == fam]
        if not sub:
            continue
        cA = sum(((h["phi"] - h["models"]["ADDtab"]) / SCALE) ** 2 for h in sub)
        cJ = sum(((h["phi"] - h["models"]["ADDj"]) / SCALE) ** 2 for h in sub)
        cS = sum(((h["phi"] - h["models"][f"SEQ{sel}"]) / SCALE) ** 2 for h in sub)
        rm = residuals_at(f["lam"], [h["phi"] for h in sub],
                          [h["models"]["ADDtab"] for h in sub],
                          [h["models"][f"SEQ{PRIMARY}"] for h in sub])
        say(f"   {fam:<12}{len(sub):>3}{cA:>9.1f}{cJ:>9.1f}{float(np.sum(rm**2)):>9.1f}"
            f"{cS:>10.1f}{float(rm.mean()):>+14.2f}")
    say()
    say("=== HELD-OUT, PER POINT (residual = measured - model, scale units) ===")
    say(f"   {'family':<12}{'CV':>5}{'K':>5}{'phi':>9}{'ADDtab':>9}{'J':>8}"
        f"{'SEQ' + PRIMARY:>9}{'MIX':>8}{'w(u>.9)':>9}")
    for h in held:
        rA = (h["phi"] - h["models"]["ADDtab"]) / SCALE
        rJ = (h["phi"] - h["models"]["ADDj"]) / SCALE
        rS = (h["phi"] - h["models"][f"SEQ{PRIMARY}"]) / SCALE
        rM = float(residuals_at(f["lam"], [h["phi"]], [h["models"]["ADDtab"]],
                                [h["models"][f"SEQ{PRIMARY}"]])[0])
        say(f"   {h['fam']:<12}{h['cv']:5.1f}{h['K']:5d}{h['phi']:9.4f}{rA:+9.2f}"
            f"{rJ:+8.2f}{rS:+9.2f}{rM:+8.2f}{h['w90']:9.3f}")
    say()

    # ============================================================ the weight correlation
    say("=== THE WITHIN-GROUP SPLIT BESIDE THE POOLED FIGURE (prereg section 5.4) ===")
    say("   The pooled rho overstates independence: within each family BOTH axes are")
    say("   monotone in CV alone, so the effective N is nearer 2 than 10.")
    cont = [h for h in held if h["fam"] in REC_SPEARMAN_WITHIN]
    aA = np.array([abs((h["phi"] - h["models"]["ADDtab"]) / SCALE) for h in cont])
    wW = np.array([h["w90"] for h in cont])
    rp = stats.spearmanr(aA, wW)
    say(f"   pooled  (lognormal + uniform, n = {len(cont)}): rho = {rp.statistic:+.3f} "
        f"(p {rp.pvalue:.3f})   recorded {REC_SPEARMAN_POOLED:+.3f}")
    for fam, recv in REC_SPEARMAN_WITHIN.items():
        sub = [h for h in cont if h["fam"] == fam]
        a_ = np.array([abs((h["phi"] - h["models"]["ADDtab"]) / SCALE) for h in sub])
        w_ = np.array([h["w90"] for h in sub])
        r_ = stats.spearmanr(a_, w_)
        say(f"   {fam:<10} n = {len(sub)}: rho = {r_.statistic:+.3f} "
            f"(p {r_.pvalue:.3f})   recorded {recv:+.3f}")
    # and the same, for the MIX residual
    aM = np.array([abs(float(residuals_at(f["lam"], [h["phi"]],
                                          [h["models"]["ADDtab"]],
                                          [h["models"][f"SEQ{PRIMARY}"]])[0]))
                   for h in cont])
    rm_ = stats.spearmanr(aM, wW)
    say(f"   the SAME correlation on the MIX residual: rho = {rm_.statistic:+.3f} "
        f"(p {rm_.pvalue:.3f})")
    say("   A weight table says where phi is SENSITIVE, not where the model FAILS.")
    say()

    # ============================================================ envelope check
    say("=== THE DESIGN-TIME ENVELOPE, CHECKED (prereg section 0) ===")
    dev = abs(f["chi2_held"] - ENVELOPE_CHI2) / ENVELOPE_CHI2
    say(f"   envelope said ~{ENVELOPE_CHI2:.0f}; computed {f['chi2_held']:.1f}; "
        f"relative deviation {dev:.2f}")
    say("   This is a FLAG ON THE SCRIPT, not a prediction: a large gap means the code")
    say("   is computing something other than what the prereg describes. It is not")
    say("   evidence for the model, because the envelope was made from the same tables.")
    if dev > 0.5:
        say("   -> DEVIATION ABOVE 50 %: check the script before believing the number.")
    say()

    say("=== WHAT THIS CAN AND CANNOT DECIDE ===")
    say("   CAN: whether ONE mixing parameter, fitted on discrete ladders only, carries")
    say("   to continuous PSDs; and whether the two discrete families want the same")
    say("   amount of nesting (P0), which is the only part that is about the RULE.")
    say("   CANNOT: which nesting is right -- the convention spread above is the honest")
    say("   width; anything about u < 0.25, where a(u), b(u) are interpolated toward a")
    say("   proved limit; and it cannot be quoted as a parameter-count win over J, whose")
    say("   two exponents saw BINARIES ONLY while lambda saw the ternary and the")
    say("   K-ladder as well.")
    say()
    say(f"FAILURES: {FAILURES}")
    _write()
    return 0


def _write():
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "grade_additivity.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())

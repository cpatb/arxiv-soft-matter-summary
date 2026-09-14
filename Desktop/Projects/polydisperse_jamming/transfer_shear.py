"""transfer_shear.py -- does the filling closure transfer to the SHEAR protocol?

THE QUESTION
------------
The pairwise-filling closure is calibrated entirely on ISOTROPIC compression at
p_hat = 1e-3, N = 1500, mu_p = 0.3. The project's HEADLINE -- the jamming shift
Delta phi_c(CV) -- is measured on a completely different protocol: constant-pressure
simple shear, rate-matched, at p* ~ 3e-4 and N = 1000, where the critical state is read
as the steady-state phi rather than as an isotropic compression endpoint.

The closure has exactly one protocol-dependent input, `beta`: the MONODISPERSE packing
density of whatever protocol it is being asked about. Its shape -- the two interaction
functions a(u) and b(u) -- is claimed to be geometric, and the mu_p = 0 block already
found the composition effect friction-independent. So the transfer test is:

    read beta on the SHEAR protocol, change NOTHING else, and ask whether the closure
    predicts the shear protocol's own phi_c(CV).

That is a cross-protocol prediction whose calibration set contains no shear point at all.

(!) THIS IS NOT A PRE-REGISTERED GATE. Both datasets were on disk and have been read
before this script was written. It is a transfer CHECK, and it is reported as one. What
it CAN establish is whether the closure's shape survives a change of protocol at all;
what it CANNOT do is be quoted as a blind prediction.

THE YARDSTICK, DECLARED BEFORE ANY NUMBER IS COMPUTED
-----------------------------------------------------
The isotropic PHI_SCALE = 2.1e-3 is the ISOTROPIC protocol's reproducibility and it is
the wrong scale here. The shear protocol's own scale is printed beside every residual:
the pooled sd over the six cells (two starting phi x three seeds) that each phi_c is
averaged over, which is 1.02e-3, 1.31e-3 and 4.34e-3 at CV 0.0, 0.2 and 0.5. Residuals
are reported in BOTH, and the SHEAR sd is the one any verdict is read against, because a
model cannot be asked to beat a measurement's own spread.

THREE CAVEATS THAT ARE PART OF THE RESULT, NOT FOOTNOTES
--------------------------------------------------------
  * The CV rungs of the rate-matched shear campaign are 0.0, 0.2, 0.5 -- NOT 0.3. The
    Status card's next-action line says "CV 0/0.3/0.5"; `out/constp_ratematch.txt` says
    0.2, and that file is the artifact. Evaluating the closure at 0.3 and grading it
    against a number measured at 0.2 would be the error this project has recorded twice.
  * N = 1000 on the shear protocol against N = 1500 for beta_iso and the whole binary
    calibration. beta is re-read on the shear protocol so the N difference is absorbed
    into it at CV 0; it is NOT absorbed at CV 0.2 and 0.5, and that is a live systematic.
  * The shear phi_c is a STEADY STATE under shear and beta_iso is a compression endpoint.
    They are different physical quantities that happen to share a symbol. The transfer
    claim is about the SHAPE of the CV dependence, not about the two betas agreeing --
    and they do not (0.5911 vs 0.5947, 1.7 isotropic scale apart).

Run:  python transfer_shear.py      (writes out/transfer_shear.txt)
"""

from __future__ import annotations

import os

import numpy as np

import grade_mixture_ext as gx
import mixture_closure as mc
import mixture_model as mm

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
MU = 0.3

# --- the shear protocol's measured critical state, from out/constp_ratematch.txt ------
# Quoted with the sd over the six cells each is pooled from, because a phi_c without its
# own spread cannot be the denominator of anything.
SHEAR_PHI_C = {0.0: (0.591127, 0.001020, 6),
               0.2: (0.601989, 0.001313, 6),
               0.5: (0.644447, 0.004336, 6)}
SHEAR_N = 1000
SHEAR_SOURCE = "out/constp_ratematch.txt, grade_constp_rm.py (R1-R5), N = 1000"

LINES: list[str] = []


def say(s: str = "") -> None:
    LINES.append(s)
    print(s)


def main() -> int:
    say("=== transfer_shear.py -- the filling closure on the SHEAR protocol ===")
    say(f"    shear source: {SHEAR_SOURCE}")
    say("    NOT a pre-registered gate: both datasets were on disk before this ran.")
    say()

    meas = mm.measured(mu_p=MU)
    beta_iso = meas[("lognormal", 0.0)][0]
    pts23 = mc.binary_points(mu_p=MU)
    allpts = gx.ext_points(beta_iso)
    cand, _, _ = gx.candidates(beta_iso, pts23)
    paJ, pbJ = cand["J"](0)
    tab, fa, fb = gx.tabulated_closure(beta_iso, allpts)

    beta_shear = SHEAR_PHI_C[0.0][0]
    say("=== the one input that changes ===")
    say(f"   beta_iso   = {beta_iso:.6f}   (isotropic compression, p_hat 1e-3, N 1500, "
        f"mu_p {MU}, n {meas[('lognormal', 0.0)][2]})")
    say(f"   beta_shear = {beta_shear:.6f}   (constant-pressure rate-matched shear, "
        f"N {SHEAR_N}, n {SHEAR_PHI_C[0.0][2]}, sd {SHEAR_PHI_C[0.0][1]:.2e})")
    say(f"   difference = {(beta_shear - beta_iso) / mm.PHI_SCALE:+.1f} isotropic scale "
        f"({beta_shear - beta_iso:+.6f}). The two betas are NOT the same number, which is")
    say("   why the closure is re-based rather than transplanted.")
    say(f"   Everything else is held FIXED: J = (PA {paJ:.4f}, PB {pbJ:.4f}); TAB at "
        f"u = {', '.join(f'{u:.3f}' for u in sorted(tab))}.")
    say("   Nothing is refitted. No bar is moved. No shear point enters any fit.")
    say()

    # --- the prediction -----------------------------------------------------------------
    say("=== absolute phi_c: measured vs predicted ===")
    say(f"   {'CV':>5}{'phi_c meas':>12}{'sd':>10}{'J':>10}{'TAB':>10}"
        f"{'J-meas':>9}{'TAB-meas':>10}{'  (in shear sd)':>16}")
    rows = []
    for cv in sorted(SHEAR_PHI_C):
        phi_m, sd, n = SHEAR_PHI_C[cv]
        d, y = mm.psd_classes("lognormal", cv)
        pj = mm.packing_density(d, y, beta_shear, paJ, pbJ)
        pt = mc.packing_density_tabulated(d, y, beta_shear, fa, fb)
        rows.append((cv, phi_m, sd, pj, pt))
        say(f"   {cv:5.1f}{phi_m:12.6f}{sd:10.2e}{pj:10.6f}{pt:10.6f}"
            f"{pj - phi_m:+9.4f}{pt - phi_m:+10.4f}"
            f"   J {(pj - phi_m) / sd:+.1f}, TAB {(pt - phi_m) / sd:+.1f}")
    say("   (CV 0 is beta_shear by construction and is not a test -- it is the input.)")
    say()

    # --- the quantity the headline actually is ------------------------------------------
    say("=== THE HEADLINE QUANTITY: Delta phi_c(CV) ===")
    say("   This is what the project quotes, and it is the right comparison: it removes")
    say("   beta from both sides, so the N = 1000 vs 1500 difference and the")
    say("   compression-endpoint vs steady-state difference cancel to first order.")
    say(f"   {'CV':>5}{'d_meas':>10}{'sd':>9}{'d_J':>10}{'d_TAB':>10}"
        f"{'J err':>9}{'TAB err':>9}{'   in shear sd':>16}")
    phi0 = SHEAR_PHI_C[0.0][0]
    sd0 = SHEAR_PHI_C[0.0][1]
    for cv, phi_m, sd, pj, pt in rows:
        if cv == 0.0:
            continue
        dm = phi_m - phi0
        # sd of a DIFFERENCE of two independent pooled means, not the per-point sd
        sem = float(np.sqrt((sd / np.sqrt(SHEAR_PHI_C[cv][2])) ** 2
                            + (sd0 / np.sqrt(SHEAR_PHI_C[0.0][2])) ** 2))
        dj = pj - beta_shear
        dt = pt - beta_shear
        say(f"   {cv:5.1f}{dm:+10.4f}{sem:9.2e}{dj:+10.4f}{dt:+10.4f}"
            f"{dj - dm:+9.4f}{dt - dm:+9.4f}"
            f"   J {(dj - dm) / sem:+.1f}, TAB {(dt - dm) / sem:+.1f}")
    say("   sd column is the standard error of the DIFFERENCE (both pooled means), not")
    say("   the per-cell sd -- lessons.md 2026-09-02.")
    say()

    # --- the same closure on its OWN protocol, as the control ---------------------------
    # --- the control, at the RIGHT pressure ----------------------------------------------
    say("=== CONTROL at the MATCHED PRESSURE (the correct comparator) ===")
    say("   🔴 An earlier version of this file used only the p_hat = 1e-3 isotropic")
    say("   ladder as its control. That is the pressure the CLOSURE is calibrated at, but")
    say("   NOT the pressure the shear protocol runs at (p* ~ 3e-4) -- so it compared a")
    say("   shear result against an isotropic result taken at a different pressure, when a")
    say("   pressure-matched isotropic ladder already exists in this campaign and is used")
    say("   by grade_constp_rm.py's own R4. Found by the physics-referee pass, 2026-09-14.")
    say("   Both controls are printed below; the MATCHED one is the comparator.")
    try:
        from grade_constp_prod import phi_iso_at_p           # noqa: PLC0415
        iso_p = phi_iso_at_p(3.0e-4)
    except Exception as exc:                                  # noqa: BLE001
        iso_p = {}
        say(f"   *** the matched ladder could not be loaded: {exc!r}")
    if iso_p and 0.0 in iso_p:
        beta_isop = float(iso_p[0.0])
        say("   phi_iso(p* = 3e-4), N = 1500, mu_p 0.3: "
            + ", ".join(f"CV {k} -> {v:.6f}" for k, v in sorted(iso_p.items())))
        say(f"   {'CV':>5}{'d_meas':>10}{'d_TAB':>10}{'% captured':>12}")
        for cv in sorted(k for k in iso_p if k != 0.0):
            d, y = mm.psd_classes("lognormal", cv)
            dt = mc.packing_density_tabulated(d, y, beta_isop, fa, fb) - beta_isop
            dm = iso_p[cv] - beta_isop
            say(f"   {cv:5.1f}{dm:+10.4f}{dt:+10.4f}{100 * dt / dm:11.0f} %")
        say("   Against the shear protocol's 87 % / 75 %, the MATCHED isotropic control")
        say("   is what the protocol comparison has to be read against -- not the")
        say("   1e-3 numbers printed below, which are at the closure's calibration")
        say("   pressure and are kept only because the closure was fitted there.")
    say()

    say("=== the UNMATCHED p_hat = 1e-3 ladder, kept for continuity, NOT the comparator ===")
    say("   Same code, same exponents, beta_iso instead of beta_shear. If the shear")
    say("   residual is much worse than this one, the protocol is the cause; if it is")
    say("   the same size, the closure's known family error is, and the shear transfer")
    say("   has added nothing.")
    say(f"   {'CV':>5}{'phi meas':>11}{'J':>10}{'TAB':>10}{'J res':>9}{'TAB res':>9}"
        f"   (isotropic scale 2.1e-3)")
    iso_res = {"J": [], "TAB": []}
    for cv in (0.0, 0.2, 0.5):
        if ("lognormal", cv) not in meas:
            say(f"   {cv:5.1f}   no isotropic lognormal record at this CV")
            continue
        phi_m = meas[("lognormal", cv)][0]
        d, y = mm.psd_classes("lognormal", cv)
        pj = mm.packing_density(d, y, beta_iso, paJ, pbJ)
        pt = mc.packing_density_tabulated(d, y, beta_iso, fa, fb)
        iso_res["J"].append((pj - phi_m) / mm.PHI_SCALE)
        iso_res["TAB"].append((pt - phi_m) / mm.PHI_SCALE)
        say(f"   {cv:5.1f}{phi_m:11.6f}{pj:10.6f}{pt:10.6f}"
            f"{(pj - phi_m) / mm.PHI_SCALE:+9.1f}{(pt - phi_m) / mm.PHI_SCALE:+9.1f}")
    say()

    # --- the comparison that decides what this block says -------------------------------
    say("=== READING ===")
    say("   The two protocols are compared on Delta phi_c, in each protocol's OWN units:")
    for cv, phi_m, sd, pj, pt in rows:
        if cv == 0.0:
            continue
        dm = phi_m - phi0
        dt = pt - beta_shear
        say(f"   CV {cv}: shear measured {dm:+.4f}, TAB predicts {dt:+.4f} "
            f"-- the model gets {100 * dt / dm:.0f} % of the measured shift")
    for cv in (0.2, 0.5):
        if ("lognormal", cv) not in meas or ("lognormal", 0.0) not in meas:
            continue
        dm_i = meas[("lognormal", cv)][0] - beta_iso
        d, y = mm.psd_classes("lognormal", cv)
        dt_i = mc.packing_density_tabulated(d, y, beta_iso, fa, fb) - beta_iso
        say(f"   CV {cv}: ISOTROPIC measured {dm_i:+.4f}, TAB predicts {dt_i:+.4f} "
            f"-- {100 * dt_i / dm_i:.0f} % of the measured shift")
    say()
    say("   A model that captures a similar FRACTION of the shift on both protocols has")
    say("   a shape that transfers and a magnitude that is wrong on both; a model that")
    say("   captures very different fractions has a protocol dependence in its shape.")
    say("   Which of those it is, is read off the two lines above -- not asserted here.")
    say()

    # --- HOW MUCH OF THIS TEST IS THE RE-BASING? ----------------------------------------
    say("=== HOW SENSITIVE IS THE PREDICTION TO beta AT ALL? ===")
    say("   The whole transfer turns on re-basing the closure from beta_iso to")
    say("   beta_shear. If the PREDICTED Delta phi barely moves when beta does, then the")
    say("   re-basing is not what is being tested and saying 'the closure transfers'")
    say("   would credit it for something it did not do. Measured, not assumed:")
    say(f"   {'CV':>5}{'d_pred(beta_iso)':>19}{'d_pred(beta_shear)':>20}{'change':>11}"
        f"{'  as a fraction of the shear measurement se':>44}")
    dbeta = beta_shear - beta_iso
    slopes = []
    for cv in (0.2, 0.5):
        d, y = mm.psd_classes("lognormal", cv)
        di = mc.packing_density_tabulated(d, y, beta_iso, fa, fb) - beta_iso
        ds = mc.packing_density_tabulated(d, y, beta_shear, fa, fb) - beta_shear
        sd = SHEAR_PHI_C[cv][1]
        sem = float(np.sqrt((sd / np.sqrt(SHEAR_PHI_C[cv][2])) ** 2
                            + (sd0 / np.sqrt(SHEAR_PHI_C[0.0][2])) ** 2))
        slopes.append((ds - di) / dbeta)
        say(f"   {cv:5.1f}{di:+19.6f}{ds:+20.6f}{ds - di:+11.6f}"
            f"{abs(ds - di) / sem:43.2f}")
    say(f"   d(Delta_pred)/d(beta) = {slopes[0]:+.3f} and {slopes[1]:+.3f}: the predicted")
    say("   shift is nearly INDEPENDENT of beta, so changing beta by the whole")
    say(f"   {abs(dbeta):.4f} between protocols moves the prediction by ~1e-4, which is")
    say("   a few per cent of the measurement's own standard error.")
    say()
    say("   CONSEQUENCE, and it weakens the reading above: this comparison CANNOT FAIL")
    say("   through the re-basing, because the re-basing does essentially nothing. What")
    say("   it actually compares is ONE fixed model prediction against TWO different")
    say("   measurements. So the honest statement is NOT 'the closure transfers across")
    say("   protocols' -- it is:")
    say("     (1) the closure's predicted Delta phi is protocol-BLIND to 1e-4;")
    say("     (2) the two protocols' MEASURED Delta phi differ from each other "
        "(isotropic")
    say("         +0.0144 vs shear +0.0109 at CV 0.2), and the single prediction sits")
    say("         between them;")
    say("     (3) it under-predicts BOTH, by 13-34 %, which is the same direction and a")
    say("         similar size to the ternary's non-additivity residual.")
    say("   A test of whether the SHAPE is protocol-dependent would have to move")
    say("   something the prediction is actually sensitive to, and beta is not it.")
    say()

    with open(os.path.join(OUT, "transfer_shear.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    print(f"-> {os.path.join(OUT, 'transfer_shear.txt')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

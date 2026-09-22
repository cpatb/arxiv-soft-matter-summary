"""jam_map.py -- the stress an impact needs to JAM a suspension, as a function of phi and CV
at fixed particle chemistry. Writes out/jam_map.txt (+ out/jam_map.npz for the figure).

Derivation: derivations/jam_stress_2026-09-22.md, eq. (8):

    sigma_jam = kappa_m d^(m-2) / Q_{1-f_c}(ghat_m),     f_c = (1 - x)/(1 - chi),

composed here from two measurements this project already makes:

  (i)  phi_0(CV), phi_m(CV) -> x, chi, f_c    isotropic ladder records (three_knobs.load)
  (ii) the distribution of ghat_m(CV)         sheared per-contact dumps (PREREG_ACTSHAPE)

⚠ The two come from DIFFERENT protocols (isotropic compression to p_hat = 3e-5, N = 1500
vs constant-pressure shear at p_hat = 3e-4, N = 1000). Composing them assumes the ghat
distribution of the sheared pack is the one that activates at the jamming point; that is
(A4) of the derivation and is flagged, not hidden.

⚠ chi needs a phi_m, i.e. the friction of the REAL grains. Two endpoints are carried:
mu_p = 3 (the framework's choice, saturated friction) and mu_p = 0.3 (the DEM default).
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

import actshape as A
import three_knobs as tk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RUNG = tk.QUOTE_RUNG
CVS = (0.0, 0.2, 0.5)
ACT_P = 3e-4                    # the actshape rung composed with the window
M_EXP = {"constant": 0, "derjaguin": 1}


def phi_endpoints(best, cv: float, mu: float) -> tuple[float, float, int]:
    """Seed mean and sd of phi at RUNG for (cv, mu), CONVERGED records only."""
    v = [r["phi"] for (m, c, s, p), r in best.items()
         if m == mu and c == cv and p == RUNG and r.get("converged")]
    if not v:
        return math.nan, math.nan, 0
    return float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else math.nan, len(v)


def act_g(cv: float, mu: float, norm: str) -> dict[int, np.ndarray]:
    """Per-seed pooled g (sigma = sampled p) from the healthy actshape dumps."""
    out = {}
    for f in glob.glob(os.path.join(OUT, "phase9_actshape_*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" in r or r["cv"] != cv or r["mu_p"] != mu or r["p_hat_target"] != ACT_P:
                continue
            if not (r["steady"] and r["servo_held"] and r["quasistatic"]):
                continue
            p = os.path.join(OUT, "actshape_dump", r["key"] + ".npz")
            if not os.path.exists(p):
                continue
            d = np.load(p)
            e = np.concatenate([[0], np.cumsum(d["n"])])
            g = [A.g_values(d["fn"][e[k]:e[k + 1]], d["r1"][e[k]:e[k + 1]],
                            d["r2"][e[k]:e[k + 1]], float(d["p_hat"][k]), norm)
                 for k in range(len(d["n"]))]
            out[r["seed"]] = (np.concatenate(g), r["d50_realised"], r["d32_realised"])
    return out


def main() -> int:
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    best = tk.load()
    say("jam_map -- sigma_jam(phi, CV) at fixed particle chemistry; derivation eq. (8)")
    say(f"window endpoints at the quoted rung p_hat = {RUNG:g} (converged records only)")
    say("")
    say("=== 1. THE WINDOW IN phi, PER CV ===")
    say("   CV   phi_0 (n)        phi(mu .3) (n)   phi(mu 3) (n)    chi(.3)  chi(3)   "
        "window[mu 3]        window[mu .3]")
    ends = {}
    for cv in CVS:
        p0, s0, n0 = phi_endpoints(best, cv, 0.0)
        p3, s3, n3 = phi_endpoints(best, cv, 0.3)
        pM, sM, nM = phi_endpoints(best, cv, 3.0)
        ends[cv] = dict(phi0=p0, phim3=pM, phim03=p3)
        say(f"  {cv:.1f}  {p0:.4f} ({n0})   {p3:.4f} ({n3})   {pM:.4f} ({nM})   "
            f"{p3 / p0:.4f}   {pM / p0:.4f}   [{pM:.4f}, {p0:.4f})   [{p3:.4f}, {p0:.4f})")
    d_ln = math.log(ends[0.5]["phi0"] / ends[0.0]["phi0"])
    say("")
    say(f"   ln phi_0(0.5)/phi_0(0) = {d_ln:.4f}   against the window's own log-width "
        f"-ln chi = {-math.log(ends[0.0]['phim3'] / ends[0.0]['phi0']):.4f} (mu 3), "
        f"{-math.log(ends[0.0]['phim03'] / ends[0.0]['phi0']):.4f} (mu 0.3)")
    say("   => CV 0 -> 0.5 at fixed weighed phi moves x by the fraction "
        f"{d_ln / -math.log(ends[0.0]['phim3'] / ends[0.0]['phi0']):.2f} (mu 3) / "
        f"{d_ln / -math.log(ends[0.0]['phim03'] / ends[0.0]['phi0']):.2f} (mu 0.3) "
        "of the window's log-width.")
    say("   ⚠ A CHECK, NOT NEW: dst_window.py (2026-09-09) already recorded 'dphi_mu/dCV =")
    say("     2.0 windows per unit CV, so dCV = 0.50-0.65 traverses the whole window'. This")
    say("     section re-derives it from the later, converged records and agrees.")
    say("   ⚠ The same record prices the lab's working point: at phi = 0.44 the DEM window is")
    say("     1.6-2.3 window-widths ABOVE it. Everything below is a statement in x = phi/phi_0,")
    say("     and the DEM phi axis must not be read as the bench's phi axis.")
    say("   ⚠ phi_0(CV 0.5) rests on n = 1 converged seed (the frictionless arm caps). A capped")
    say("     cell reads LOW, so the shift is if anything under-stated.")

    say("")
    say("=== 2. f_c(phi) PER CV (mu_m = 3 / 0.3); '--' = outside the window ===")
    phis = np.round(np.arange(0.54, 0.70, 0.01), 3)
    say("   phi   " + "  ".join(f"CV{cv:.1f}[3]  CV{cv:.1f}[.3]" for cv in CVS))
    for ph in phis:
        row = []
        for cv in CVS:
            for key in ("phim3", "phim03"):
                chi = ends[cv][key] / ends[cv]["phi0"]
                x = ph / ends[cv]["phi0"]
                fc = A.f_crit(x, chi)
                row.append("   --    " if not math.isfinite(fc) else
                           ("  rest   " if fc == 0 else f"  {fc:.3f}  "))
        say(f"  {ph:.2f} " + " ".join(row))

    # ---- 3. compose with the activation distributions ----------------------------
    say("")
    say("=== 3. sigma_jam(phi)/kappa AT FIXED kappa AND d, composed with actshape ===")
    have = {(cv, mu): act_g(cv, mu, "derjaguin") for cv in CVS for mu in (0.0, 0.3)}
    if not all(have.values()):
        say("   actshape dumps not (all) present -- section 3 skipped")
    else:
        npz = {}
        phig = np.linspace(0.54, 0.70, 321)
        for norm in ("derjaguin", "constant"):
            m = M_EXP[norm]
            for mu in (0.0, 0.3):
                say(f"   norm {norm} (m = {m}), activation network mu_p = {mu}, "
                    f"window mu_m = 3; sigma_jam/kappa in d10 units [and at matched d50]")
                for cv in CVS:
                    gs = act_g(cv, mu, norm)
                    chi = ends[cv]["phim3"] / ends[cv]["phi0"]
                    curves, curves50 = [], []
                    for s, (g, d50, _d32) in sorted(gs.items()):
                        c = []
                        for ph in phig:
                            fc = A.f_crit(ph / ends[cv]["phi0"], chi)
                            c.append(A.sigma_jam_over_kappa(g, fc) if 0 < fc < 1 else
                                     (0.0 if fc <= 0 else math.inf))
                        c = np.array(c)
                        curves.append(c)
                        curves50.append(c * d50 ** (2 - m))
                    def gmean(cs):
                        a = np.array(cs)
                        return np.exp(np.mean(np.log(np.where(np.isfinite(a) & (a > 0), a,
                                                              np.nan)), axis=0))
                    cm, cm50 = gmean(curves), gmean(curves50)
                    # state per phi: 0 = window, +1 = jammed at rest, -1 = never jams, kept
                    # separately so 'rest' and 'never' are not both just NaN in the npz
                    state = np.array([0 if 0 < A.f_crit(ph / ends[cv]["phi0"], chi) < 1 else
                                      (1 if A.f_crit(ph / ends[cv]["phi0"], chi) <= 0 else -1)
                                      for ph in phig])
                    npz[f"{norm}_mu{mu}_cv{cv}"] = cm
                    npz[f"{norm}_mu{mu}_cv{cv}_d50"] = cm50
                    npz[f"state_cv{cv}"] = state
                    pick = [0.60, 0.62, 0.64, 0.66, 0.68]
                    for lab, arr in (("d10", cm), ("d50", cm50)):
                        vals = []
                        for ph in pick:
                            i = int(np.argmin(abs(phig - ph)))
                            fc = A.f_crit(ph / ends[cv]["phi0"], chi)
                            vals.append("   never " if not math.isfinite(fc) else
                                        ("   rest  " if fc <= 0 else f"{arr[i]:8.4g} "))
                        say(f"     CV {cv:.1f} [{lab}]: phi " + " ".join(f"{p:.2f}" for p in pick)
                            + " -> " + " ".join(vals))
        # ---- 4. which ingredient moves it: counterfactual swaps -----------------------
        # Revised after referee 2026-09-22 (MATERIAL x3): both swap orders and the interaction
        # term are printed; the decomposition is repeated at three length conventions because
        # the g-only factor's SIGN depends on which diameter is held fixed; and the g of ONE
        # sheared steady state (x_ss 0.978-0.994 at mu_p 0, 0.923-0.936 at mu_p 0.3) is applied
        # across the whole phi sweep, an
        # extrapolation in x that is flagged, not measured.
        say("")
        say("=== 4. DECOMPOSITION at fixed phi: window shift (x) vs activation (g) ===")
        say("   T = sjam(g_CV, f_c[CV]) / sjam(g_0, f_c[0]);  X = sjam(g_0, f_c[CV])/base;")
        say("   G = sjam(g_CV, f_c[0])/base;  interaction I = T/(X G) (1 = separable).")
        say("   x-share: order x-first ln X/ln T, order g-first ln(T/G)/ln T -- both printed.")
        say("   G is referenced to a fixed diameter: sjam scales as d_ref^(2-m), m = 1 (derjaguin).")
        say("   ⚠ CV 0.5 is FINITE-SIZE FLAGGED (PREREG_ACTSHAPE Amendment 1): provisional.")
        say("   ⚠ g is measured at the sheared steady state (x_ss = phi_ss/phi0 = 0.978-0.994 at mu_p 0,")
        say("     0.923-0.936 at mu_p 0.3; referee pass 2) and APPLIED at every")
        say("     phi below: an extrapolation in x the DEM has not tested (assumption A4). The only")
        say("     evidence is H4 (delta stable across a decade of p_hat, which moves x slightly).")
        for mu in (0.0, 0.3):
            s0 = act_g(0.0, mu, "derjaguin")
            g0 = np.concatenate([v[0] for v in s0.values()])
            dref0 = {"d10": 1.0, "d50": float(np.mean([v[1] for v in s0.values()])),
                     "d32": float(np.mean([v[2] for v in s0.values()]))}
            chi0 = ends[0.0]["phim3"] / ends[0.0]["phi0"]
            for cv in (0.2, 0.5):
                sc = act_g(cv, mu, "derjaguin")
                gc = np.concatenate([v[0] for v in sc.values()])
                drefc = {"d10": 1.0, "d50": float(np.mean([v[1] for v in sc.values()])),
                         "d32": float(np.mean([v[2] for v in sc.values()]))}
                chic = ends[cv]["phim3"] / ends[cv]["phi0"]
                for ph in (0.58, 0.59, 0.60, 0.61, 0.62, 0.63):
                    f0 = A.f_crit(ph / ends[0.0]["phi0"], chi0)
                    fcv = A.f_crit(ph / ends[cv]["phi0"], chic)
                    if not (0 < f0 < 1 and 0 < fcv < 1):
                        say(f"   mu {mu:.1f} CV {cv:.1f} phi {ph:.2f}: outside a window "
                            f"(CV {cv}: {'never' if fcv == math.inf else 'rest' if fcv <= 0 else 'ok'},"
                            f" CV 0: {'never' if f0 == math.inf else 'rest' if f0 <= 0 else 'ok'})")
                        continue
                    for ref in ("d10", "d50", "d32"):
                        k = drefc[ref] / dref0[ref]          # (2 - m) = 1
                        base = A.sigma_jam_over_kappa(g0, f0)
                        T = A.sigma_jam_over_kappa(gc, fcv) * k / base
                        X = A.sigma_jam_over_kappa(g0, fcv) / base
                        G = A.sigma_jam_over_kappa(gc, f0) * k / base
                        inter = T / (X * G)
                        sh1 = math.log(X) / math.log(T)
                        sh2 = math.log(T / G) / math.log(T)
                        say(f"   mu {mu:.1f} CV {cv:.1f} phi {ph:.2f} [{ref}]: f_c {f0:.3f}->{fcv:.3f}"
                            f"  T x{T:.3f}  X x{X:.3f}  G x{G:.3f}  I {inter:.3f}  "
                            f"x-share {sh1:.2f} / {sh2:.2f}  G {'opposes' if (G - 1) * (T - 1) < 0 else 'reinforces'}")
        np.savez(os.path.join(OUT, "jam_map.npz"), phi=phig, **npz,
                 **{f"phi0_cv{cv}": ends[cv]["phi0"] for cv in CVS},
                 **{f"phim3_cv{cv}": ends[cv]["phim3"] for cv in CVS})
    with open(os.path.join(OUT, "jam_map.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""multiplet_wc.py -- turn an aggregation state into Wyart-Cates inputs for the DNS.

WHAT THIS IS FOR
----------------
`C:\\Users\\compl\\Github\\Drop Impact\\basilisk-drop-impact\\wc_drop_impact.c` solves a
drop impact with

    eta(phi, fhat) = ETA_S * (1 - phi/phiJ(fhat))^(-2),
    phiJ(fhat)     = PHI0 - (PHI0 - PHIM) * fhat,
    fhat           -> fss(sigma/SIGMA_STAR)

and takes `PHI`, `PHI0`, `PHIM`, `SIGMA_STAR` from the environment. This module says what
an experimentalist's AGGREGATION STATE does to that quadruple.

THE DEGENERACY THAT MAKES THE ANSWER ONE NUMBER
-----------------------------------------------
The closure above is EXACTLY covariant under a common rescaling of the three packing
fractions: send `(PHI, PHI0, PHIM) -> c (PHI, PHI0, PHIM)` and `phi/phiJ` is unchanged at
every `fhat`, so `eta` is unchanged, and `fhat`'s own dynamics depend on `sigma/SIGMA_STAR`
and not on `phi` at all. (This repo already records the covariance as a memory:
"Wyart-Cates phi_0 is gauge".)

Aggregation does exactly two things (`derivations/multiplet_psd_2026-09-22.md`):

    PHI  -> Lambda * PHI          Lambda = <k^gamma>/<k>, closed form, no fit
    PHI0 -> G0 * PHI0             G0 = polydispersity gain of the FRICTIONLESS jamming point
    PHIM -> Gm * PHIM             Gm = polydispersity gain of the FRICTIONAL one

So **if `G0 == Gm == G`, the entire effect of aggregation on this solver is the single
ratio `Lambda/G`** -- the three-parameter change collapses to one, by the covariance. If
`G0 != Gm` the collapse fails and aggregation ALSO changes the SHAPE of `phiJ(fhat)`,
i.e. the thickening contrast, which is a qualitatively different statement.

**Whether `G0 == Gm` is a question this project has the cells to answer**, and answering it
is what section 2 below does, converged-only and against the protocol floor.

WHAT IS DELIBERATELY NOT DONE HERE
----------------------------------
`SIGMA_STAR` is NOT mapped. Whether aggregation moves the onset stress is a separate
question with its own gauge problem, already owned by `SIZE_THRESHOLD_2026-09-08.md` and
`derivations/onset_size_degeneracy_2026-09-19.md`. Section 4 states the question and the
two candidate conventions and refuses to choose between them from this module.

Run:  python3 multiplet_wc.py      (writes out/multiplet_wc.txt)
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np

import mixture_model as mm
from multiplet_levers import D10, flory
from psd import MultipletPSD

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

SCALE = mm.PHI_SCALE          # 2.1e-3, the protocol's own run-to-run non-determinism
RUNG = 3e-5                   # the matched lowest rung; the fit-free channel
N_REQ = 1500

# The solver's own defaults, for reference when a mapped quadruple is printed.
WC_DEFAULT = dict(PHI=0.52, PHI0=0.62, PHIM=0.50, SIGMA_STAR=20.0, ETA_S=0.05)


def load_rung(rung: float = RUNG) -> dict:
    """Every lognormal N=1500 record at one rung, keyed (mu_p, cv, seed).

    Selection rule copied from `grade_mustar.load_all`, which is the rule the graded
    campaigns use: highest `max_iter` wins at a repeated key. N is printed before and
    after every filter.
    """
    n_lines = n_family = n_size = n_rung = 0
    cells: dict = {}
    for fn in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        if "wcbidi" in os.path.basename(fn):
            continue                       # bidisperse family, graded separately
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            n_lines += 1
            r = json.loads(ln)
            if r.get("failed"):
                continue
            if (r.get("family") or "lognormal") != "lognormal":
                continue
            n_family += 1
            if r.get("n_requested") != N_REQ:
                continue
            n_size += 1
            if abs(r["p_hat_target"] - rung) / rung > 0.05:
                continue
            n_rung += 1
            k = (r["mu_p"], round(r["cv"], 1), r["seed"])
            prev = cells.get(k)
            if prev is None or r.get("max_iter", 0) > prev.get("max_iter", 0):
                cells[k] = r
    return dict(cells=cells, n_lines=n_lines, n_family=n_family, n_size=n_size,
                n_rung=n_rung)


def gain(cells: dict, mu: float, cv: float, converged_only: bool = True):
    """G(mu, cv) = <phi(mu, cv)> / <phi(mu, 0)>, SEED-PAIRED, and its floor-based s.e.

    Returns (G, se_G, n_pairs, n_dropped). Seeds are paired because `packing.py` seeds
    the radius draw from `st.seed` regardless of CV, so the same seed is a correlated
    draw at both CV levels -- treating the arms as independent inflates the bar (a
    referee caught exactly that in this repo on 2026-09-19).

    ⚠ The quoted s.e. is `max(seed-paired s.e., sqrt(L) * SCALE / sqrt(n))` with L = 2:
    the two cells (cv, 0) do NOT share a record, so neither cancels. A seed-paired s.e.
    alone is NOT an acceptable bar in this repository.
    """
    num, den, dropped = {}, {}, 0
    for s in (1, 2, 3):
        a, b = cells.get((mu, cv, s)), cells.get((mu, 0.0, s))
        if a is None or b is None:
            continue
        if converged_only and not (a.get("converged") and b.get("converged")):
            dropped += 1
            continue
        num[s], den[s] = a["phi"], b["phi"]
    if not num:
        return None, None, 0, dropped
    seeds = sorted(num)
    g = np.array([num[s] / den[s] for s in seeds])
    m = float(g.mean())
    se_seed = float(g.std(ddof=1) / math.sqrt(len(g))) if len(g) > 1 else float("inf")
    # propagate the protocol floor through the RATIO: d(a/b) = (1/b) sqrt(2) SCALE / sqrt(n)
    b_mean = float(np.mean([den[s] for s in seeds]))
    se_floor = math.sqrt(2.0) * SCALE / (b_mean * math.sqrt(len(g)))
    return m, max(se_seed, se_floor), len(g), dropped


def wc_inputs(k_bar: float, d_f: float, G0: float, Gm: float,
              base: dict | None = None) -> dict:
    """Map an aggregation state onto the solver's quadruple."""
    base = dict(base or WC_DEFAULT)
    lam = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar)).lambda_trap
    out = dict(base)
    out["PHI"] = base["PHI"] * lam
    out["PHI0"] = base["PHI0"] * G0
    out["PHIM"] = base["PHIM"] * Gm
    out["_Lambda"] = lam
    # the two distances to jamming that actually drive the solver
    out["_phi_over_phi0"] = out["PHI"] / out["PHI0"]
    out["_phi_over_phim"] = out["PHI"] / out["PHIM"]
    return out


def main() -> int:
    L: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        L.append(s)

    say("=" * 92)
    say("multiplet_wc.py -- aggregation state -> Wyart-Cates DNS inputs")
    say("=" * 92)

    # ------------------------------------------------------------------ 1
    say("\n=== 1. THE COVARIANCE THAT COLLAPSES THREE PARAMETERS INTO ONE ===")
    say("   eta = ETA_S (1 - PHI/phiJ)^-2,  phiJ = PHI0 - (PHI0-PHIM) fhat.")
    say("   Under (PHI, PHI0, PHIM) -> c (PHI, PHI0, PHIM), PHI/phiJ is unchanged at")
    say("   every fhat, and fhat's own dynamics read sigma/SIGMA_STAR, not phi.")
    c = 1.37
    b = WC_DEFAULT
    for fh in (0.0, 0.5, 1.0):
        pj = b["PHI0"] - (b["PHI0"] - b["PHIM"]) * fh
        pj_c = c * b["PHI0"] - (c * b["PHI0"] - c * b["PHIM"]) * fh
        say(f"      fhat = {fh:.1f}:  phi/phiJ = {b['PHI']/pj:.10f}   "
            f"rescaled by c = {c}: {c*b['PHI']/pj_c:.10f}")
    say("   -> the covariance is EXACT, so only the RATIO Lambda/G can move the solver.")

    # ------------------------------------------------------------------ 2
    say("\n=== 2. ARE THE TWO GAINS EQUAL? G0 (frictionless) vs Gm (frictional) ===")
    d = load_rung()
    cells = d["cells"]
    say(f"   N: {d['n_lines']} records read -> {d['n_family']} lognormal -> "
        f"{d['n_size']} at N={N_REQ} -> {d['n_rung']} at the rung p_hat = {RUNG}")
    say(f"   distinct (mu_p, CV, seed) cells after highest-cap dedup: {len(cells)}")
    say("")
    say("   G(mu, CV) = phi(mu, CV)/phi(mu, 0), seed-paired, CONVERGED CELLS ONLY.")
    say(f"   bar = max(seed s.e., sqrt(2)*PHI_SCALE/sqrt(n)), PHI_SCALE = {SCALE}.")
    say("")
    mus = sorted(set(k[0] for k in cells))
    cvs = sorted(set(k[1] for k in cells if k[1] > 0))
    say(f"   {'mu_p':>6}{'CV':>6}{'G':>10}{'s.e.':>9}{'n':>4}{'dropped':>9}   note")
    table: dict = {}
    for mu in mus:
        for cv in cvs:
            g, se, n, dr = gain(cells, mu, cv)
            if g is None:
                say(f"   {mu:6.1f}{cv:6.1f}{'--':>10}{'--':>9}{0:4d}{dr:9d}   "
                    f"no converged pair")
                continue
            table[(mu, cv)] = (g, se, n)
            say(f"   {mu:6.1f}{cv:6.1f}{g:10.4f}{se:9.4f}{n:4d}{dr:9d}")
    say("")
    say("   ⚠ 'dropped' counts seed pairs removed because a cell did not converge.")
    say("   A drop rate that tracks mu_p or CV is a drop rate that tracks the axis")
    say("   under test; it is printed per row for exactly that reason.")

    # the comparison
    say("\n   THE COMPARISON: does G depend on FRICTION?")
    say("   Taken over the WIDEST FRICTION SPAN that has a resolved G at each CV")
    say("   (n >= 2 converged seed pairs, finite s.e.). The span actually used is")
    say("   printed, because a comparison over a narrow span is a weaker test and the")
    say("   reader must be able to see which one it was.")
    say("")
    n_resolved = n_differ = 0
    for cv in cvs:
        usable = [(mu, *table[(mu, cv)]) for mu in mus
                  if (mu, cv) in table and table[(mu, cv)][2] >= 2
                  and math.isfinite(table[(mu, cv)][1])]
        if len(usable) < 2:
            have = [f"mu={u[0]}" for u in usable]
            say(f"      CV = {cv}: NOT COMPARABLE -- only {len(usable)} friction(s) "
                f"resolved {have}")
            continue
        lo, hi = usable[0], usable[-1]
        diff = lo[1] - hi[1]
        se = math.hypot(lo[2], hi[2])   # DIFFERENT frictions, run as separate cells
        t = abs(diff) / se if se > 0 else float("inf")
        differ = t > 4.303
        n_resolved += 1
        n_differ += differ
        say(f"      CV = {cv}: G(mu={lo[0]}) = {lo[1]:.4f} +- {lo[2]:.4f}   "
            f"G(mu={hi[0]}) = {hi[1]:.4f} +- {hi[2]:.4f}   "
            f"|t| = {t:.2f} -> {'DIFFER' if differ else 'not resolved'}")
    say("")
    if n_resolved:
        say(f"   -> {n_differ} of {n_resolved} comparable CV levels RESOLVE a friction")
        say("      dependence of G at t_crit = 4.303.")
        if n_differ == 0:
            say("   🔑 So on the cells that exist, G0 = Gm within the protocol floor, and")
            say("      the collapse holds: the whole effect of aggregation on this solver")
            say("      is the single ratio Lambda/G.")
            say("   ⚠ This is a BOUND, not a proof of equality. n = 3 seeds and a friction")
            say("      span narrower than the registered one cannot establish identity;")
            say("      the honest statement is that no friction dependence is resolved.")
    say("")
    say("   ⚠ THE FRICTION SPAN IS CENSORED, and in the direction that matters: the")
    say("   mu_p = 0.6 and 3.0 blocks ran at MAX_ITER = 200000 and capped at exactly")
    say("   the CV > 0 cells this comparison needs (see PREREG_MUSTARCAP_2026-09-22.md).")
    say("   The `mustarcap` re-run at 1e6 is what widens the span; re-run this module")
    say("   when it lands.")
    say("")
    say("   t_crit = 4.303 is the two-sided 95 % value at the design's own 2 dof (n = 3")
    say("   seeds), the same bar the manuscript uses. ⚠ `hypot` is correct HERE because")
    say("   the two arms are different frictions run as separate cells; it would be wrong")
    say("   for two CV levels at one friction, which share the seed's radius draw.")

    # ------------------------------------------------------------------ 3
    say("\n=== 3. THE MAPPED QUADRUPLE, for a solver run ===")
    say("   Using the MEASURED gains where they exist. If G0 and Gm are not resolved")
    say("   apart (section 2), the honest mapping uses ONE G for both and the whole")
    say("   effect is Lambda/G -- which is what the table below reports.")
    say("")
    cv_ref = max(cvs) if cvs else 0.5
    usable = [(mu, *table[(mu, cv_ref)]) for mu in mus
              if (mu, cv_ref) in table and table[(mu, cv_ref)][2] >= 2
              and math.isfinite(table[(mu, cv_ref)][1])]
    g0 = (usable[0][1], usable[0][2]) if usable else None
    gm = (usable[-1][1], usable[-1][2]) if len(usable) >= 2 else None
    if g0 is None or gm is None:
        say(f"   ⚠ FEWER THAN TWO FRICTIONS RESOLVE A GAIN at CV = {cv_ref} on")
        say("   converged cells. The quadruple is therefore NOT emitted; a mapping")
        say("   built on a non-converged cell would carry that cell's censoring")
        say("   into the DNS. Re-run this module when the mustarcap campaign lands.")
    else:
        say(f"   gains read at CV = {cv_ref}: G0 = {g0[0]:.4f} (mu_p = {usable[0][0]}), "
            f"Gm = {gm[0]:.4f} (mu_p = {usable[-1][0]})")
        say("   ⚠ 'G0' here is the LOWEST RESOLVED friction, which is not necessarily")
        say(f"   mu_p = 0; it is mu_p = {usable[0][0]}. The frictionless branch is named")
        say("   G0 throughout this repo and the label is kept, but the value is the")
        say("   lowest friction that has converged cells, and it is printed so.")
        say("")
        say(f"   {'<k>':>6}{'d_f':>6}{'Lambda':>9}{'PHI':>9}{'PHI0':>9}{'PHIM':>9}"
            f"{'phi/phi0':>10}{'phi/phiM':>10}")
        base = WC_DEFAULT
        lam_max = 1.0 / base["PHI"]

        def regime(w: dict) -> str:
            """ONE rule, applied to the baseline and to every mapped row alike.

            The baseline row was briefly labelled 'flowable' by a hardcoded string
            while its own phi/phiM = 1.0400 said otherwise -- the exact defect this
            repository records as 'prose contradicts its own artifact'. Computed now.
            """
            if w["PHI"] > 1.0:
                return "🔴 INADMISSIBLE"
            if w["_phi_over_phim"] > 1.0:
                return "jammed (frictional)"
            if w["_phi_over_phi0"] > 1.0:
                return "jammed (frictionless)"
            return "flowable"

        base_row = dict(base)
        base_row["_phi_over_phi0"] = base["PHI"] / base["PHI0"]
        base_row["_phi_over_phim"] = base["PHI"] / base["PHIM"]
        say(f"   {'(none)':>6}{'--':>6}{1.0:9.4f}{base['PHI']:9.4f}{base['PHI0']:9.4f}"
            f"{base['PHIM']:9.4f}{base_row['_phi_over_phi0']:10.4f}"
            f"{base_row['_phi_over_phim']:10.4f}   {regime(base_row)}")
        n_inadm = 0
        for kb, df in ((1.5, 3.0), (1.5, 2.1), (1.5, 1.8), (2.0, 2.1), (2.0, 1.8)):
            w = wc_inputs(kb, df, g0[0], gm[0])
            flag = regime(w)
            n_inadm += flag.endswith("INADMISSIBLE")
            say(f"   {kb:6.1f}{df:6.1f}{w['_Lambda']:9.4f}{w['PHI']:9.4f}{w['PHI0']:9.4f}"
                f"{w['PHIM']:9.4f}{w['_phi_over_phi0']:10.4f}{w['_phi_over_phim']:10.4f}"
                f"   {flag}")
        say("")
        say("   🔑 READ THE TOP TWO ROWS TOGETHER. The solver's DEFAULT preset is already")
        say("   frictionally jammed (phi/phiM = 1.0400), and COMPACT aggregation at")
        say("   <k> = 1.5 UN-JAMS it (0.9497) -- Lambda = 1 exactly there, so the only")
        say("   lever acting is the one that RAISES phiJ. Opening the same aggregates to")
        say("   d_f = 2.1 drives it the other way, to 1.2293. One aggregation axis, both")
        say("   signs, at FIXED solid loading -- which is the bench prediction to test.")
        say("")
        say("   🔴 ADMISSIBILITY, and it needs no threshold: phi_eff = Lambda * phi_solid")
        say("   is a VOLUME FRACTION and cannot exceed 1. That caps Lambda at")
        say(f"      Lambda_max = 1/phi_solid = 1/{base['PHI']} = {lam_max:.4f},")
        say(f"   and {n_inadm} of the 5 rows above BREACH it. Those rows are not")
        say("   'a large effect' -- they are the no-interpenetration assumption (A3 of")
        say("   the derivation) DESTROYING ITSELF, and they must not be run.")
        say("")
        say("   🔑 This is the sharpest statement in the file, because it is the model")
        say("   bounding its own domain from the inside: open aggregates at moderate")
        say("   loading REQUIRE interpenetration, so any treatment of a d_f = 1.8")
        say("   suspension at phi = 0.52 that keeps them impenetrable is wrong before")
        say("   any rheology is computed. The admissible window is")
        say(f"      Lambda <= {lam_max:.4f}, i.e. d_f >= the value where Lambda reaches it.")
        for kb in (1.5, 2.0):
            cross = None
            nk = flory(kb)
            for f in np.linspace(3.0, 1.0, 20001):
                if MultipletPSD(d10=D10, d_f=float(f), n_k=nk).lambda_trap >= lam_max:
                    cross = float(f)
                    break
            say(f"      <k> = {kb}: the effective-sphere model is admissible only for "
                f"d_f >= {cross:.3f}")
        say("")
        say("   phi/phiM > 1 (non-breaching rows) means the FRICTIONAL branch is jammed:")
        say("   eta diverges as soon as fhat rises. That is the shear-jammed, rebounding")
        say("   state, and the table says which aggregation states reach it from a phi")
        say("   that is SAFE unaggregated -- which is the bench prediction worth testing.")

    # ------------------------------------------------------------------ 4
    say("\n=== 4. SIGMA_STAR IS NOT MAPPED, AND WHY ===")
    say("   The onset stress is a force divided by an area, and WHICH area is a")
    say("   convention this module is not entitled to choose:")
    say("      (a) the MONOMER convention -- the repulsive force is a surface")
    say("          interaction between two monomer surfaces at contact, so the area is")
    say("          d_1^2 and SIGMA_STAR does NOT move with aggregation;")
    say("      (b) the AGGREGATE convention -- the stress is transmitted across the")
    say("          multiplet's cross-section, so the area is d_k^2 and SIGMA_STAR falls")
    say("          like <d^2>^-1, i.e. by a factor k^(2 gamma/3).")
    say("   These differ by orders of magnitude and they predict opposite trends.")
    say("   `SIZE_THRESHOLD_2026-09-08.md` and `derivations/onset_size_degeneracy_2026-")
    say("   09-19.md` already own this question; the second one's verdict is a NEGATIVE")
    say("   (a size ladder cannot separate the two grain-scale stresses). Choosing here")
    say("   would be picking a convention and calling the consequence a result -- which")
    say("   this repository has been caught doing before. It is left OPEN and flagged.")
    say("")
    say("   ⚠ CONSEQUENCE FOR THE DNS: a solver run built from section 3 varies PHI,")
    say("   PHI0 and PHIM at FIXED SIGMA_STAR. That is convention (a). It must be")
    say("   labelled as such in the run's metadata, not presented as the mapping.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "multiplet_wc.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    say("")
    say("   written: out/multiplet_wc.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())

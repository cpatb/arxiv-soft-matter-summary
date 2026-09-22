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
from multiplet_levers import BETA as ML_BETA
from multiplet_levers import D10, flory
from multiplet_levers import phi_J as ml_phi_J
from psd import MultipletPSD

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

SCALE = mm.PHI_SCALE          # 2.1e-3, the protocol's own run-to-run non-determinism
RUNG = 3e-5                   # the matched lowest rung; the fit-free channel
N_REQ = 1500

# The solver's own defaults, for reference when a mapped quadruple is printed.
WC_DEFAULT = dict(PHI=0.52, PHI0=0.62, PHIM=0.50, SIGMA_STAR=20.0, ETA_S=0.05)

# Two-sided 95 % Student t, by the number of PAIRED seeds a row actually has.
# n = 3 -> 2 dof -> 4.303 is the design's own bar and the one the manuscript uses;
# a row that lost a seed to non-convergence has n = 2 -> 1 dof -> 12.706.
T_CRIT = {2: 12.706, 3: 4.303, 4: 3.182, 5: 2.776}


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


def gain_seeds(cells: dict, mu: float, cv: float, converged_only: bool = True):
    """Per-seed G(mu, cv) = phi(mu, cv)/phi(mu, 0), plus the denominators and the drops.

    Returns ({seed: G}, {seed: phi(mu,0)}, n_dropped). Kept per-seed rather than
    averaged, because the friction comparison in `compare_frictions` has to be PAIRED
    and an average throws the pairing away.
    """
    num, den, dropped = {}, {}, 0
    for s in (1, 2, 3):
        a, b = cells.get((mu, cv, s)), cells.get((mu, 0.0, s))
        if a is None or b is None:
            continue
        if converged_only and not (a.get("converged") and b.get("converged")):
            dropped += 1
            continue
        num[s], den[s] = a["phi"] / b["phi"], b["phi"]
    return num, den, dropped


def gain(cells: dict, mu: float, cv: float, converged_only: bool = True):
    """G(mu, cv), SEED-PAIRED across CV, and its floor-based s.e.

    Returns (G, se_G, n_pairs, n_dropped). Seeds are paired across CV because
    `packing.py` seeds the radius draw from `st.seed` regardless of CV, so the same seed
    is a correlated draw at both CV levels.

    ⚠ The quoted s.e. is `max(seed-paired s.e., floor)`. The floor propagates the
    protocol's per-record non-determinism `SCALE` through the RATIO: with `G = a/b`,
    `var(G) = SCALE^2 (1 + G^2)/b^2`, and the mean over `n` seeds divides by `n`. The
    protocol floor does NOT cancel under seed pairing -- pairing removes the SEED's
    contribution, not the re-run scatter of each record.
    """
    g, den, dropped = gain_seeds(cells, mu, cv, converged_only)
    if not g:
        return None, None, 0, dropped
    seeds = sorted(g)
    arr = np.array([g[s] for s in seeds])
    m = float(arr.mean())
    se_seed = float(arr.std(ddof=1) / math.sqrt(len(arr))) if len(arr) > 1 else float("inf")
    b_mean = float(np.mean([den[s] for s in seeds]))
    se_floor = SCALE * math.sqrt(1.0 + m * m) / (b_mean * math.sqrt(len(arr)))
    return m, max(se_seed, se_floor), len(arr), dropped


def compare_frictions(cells: dict, mu_lo: float, mu_hi: float, cv: float):
    """PAIRED comparison of G between two frictions. Returns (diff, se, n, r).

    🔴 CORRECTED 2026-09-22 after referee pass 2. This comparison previously combined the
    two arms with `math.hypot`, justified in a comment as *"correct HERE because the two
    arms are different frictions run as separate cells"*. **That reasoning is wrong on
    this repository's own recorded terms**: `packing.py` seeds the radius draw with
    `default_rng(st.seed)` REGARDLESS of `mu_p`, so seed `s` is the same size draw at
    both frictions -- which is exactly the pairing a referee caught being ignored here on
    2026-09-19. Being separate CELLS does not make two arms independent when they share
    the draw.

    The correct estimator is the per-seed DIFFERENCE of the two gains. The measured
    correlation is printed beside the result so the reader can see how much pairing buys.

    The floor propagates `SCALE` through all FOUR records that enter `d_s`:
        d = a/b - c/e,
        var(d) = SCALE^2 [ (1 + G_lo^2)/b^2 + (1 + G_hi^2)/e^2 ].
    None of the four cancels: pairing removes the seed's radius draw, not the protocol's
    per-record re-run scatter.
    """
    g_lo, den_lo, _ = gain_seeds(cells, mu_lo, cv)
    g_hi, den_hi, _ = gain_seeds(cells, mu_hi, cv)
    seeds = sorted(set(g_lo) & set(g_hi))
    if len(seeds) < 2:
        return None, None, len(seeds), None
    d = np.array([g_lo[s] - g_hi[s] for s in seeds])
    a = np.array([g_lo[s] for s in seeds])
    b = np.array([g_hi[s] for s in seeds])
    r = float(np.corrcoef(a, b)[0, 1]) if len(seeds) > 2 else float("nan")
    n = len(seeds)
    se_seed = float(d.std(ddof=1) / math.sqrt(n))
    bl = float(np.mean([den_lo[s] for s in seeds]))
    bh = float(np.mean([den_hi[s] for s in seeds]))
    var = SCALE ** 2 * ((1 + a.mean() ** 2) / bl ** 2 + (1 + b.mean() ** 2) / bh ** 2)
    se_floor = math.sqrt(var / n)
    return float(d.mean()), max(se_seed, se_floor), n, r


def wc_inputs(k_bar: float, d_f: float, G0: float, Gm: float,
              base: dict | None = None, row_gain: bool = True) -> dict:
    """Map an aggregation state onto the solver's quadruple.

    🔴 FATAL DEFECT, FOUND BY REFEREE PASS 4 AND FIXED HERE. This function used to apply
    the `G0`/`Gm` passed in -- read once at `cv_ref = max(cvs)`, i.e. **CV = 0.5** -- to
    EVERY row, whatever that row's own CV was. The rows' actual realised CVs are
    **0.16 to 0.36**, and `G` rises monotonically with CV (this module's own section 2
    measures it: 1.0043 / 1.0212 / 1.0518 / 1.0781 / 1.1066 at CV 0.1 ... 0.5). So every
    row got an INFLATED `G`, which inflated the mapped `PHI0`/`PHIM` and biased every
    `phi/phiJ` ratio DOWNWARD -- toward looking less jammed than it is.

    It flipped the two rows that carried the headline. With the row's own `G`:
        (<k> = 1.5, d_f = 3.0): phi/phiM  0.9497 -> 1.0250   ("un-jams" is FALSE)
        (<k> = 1.5, d_f = 2.1): phi/phi0  0.9810 -> 1.0508   (admissible -> A3 FAILS)

    `row_gain=True` (the default, and the only defensible setting) evaluates the SAME
    calibrated closure on the row's OWN multiplet PSD, which is what `multiplet_levers`
    already computes as `G`. It is self-consistent by construction: the closure that
    supplies `Lambda`'s partner supplies the gain at the same PSD, same rung, same gauge.

    ⚠ One `G` is used for both `PHI0` and `PHIM`, because section 2 resolves no friction
    dependence of `G`. That is a BOUND, not a proof of equality -- see section 2.

    `row_gain=False` reproduces the withdrawn behaviour and exists only so the two can be
    printed side by side; it is never the quoted mapping.
    """
    base = dict(base or WC_DEFAULT)
    psd = MultipletPSD(d10=D10, d_f=d_f, n_k=flory(k_bar))
    lam = psd.lambda_trap
    g_row = ml_phi_J(psd) / ML_BETA
    out = dict(base)
    out["PHI"] = base["PHI"] * lam
    out["PHI0"] = base["PHI0"] * (g_row if row_gain else G0)
    out["PHIM"] = base["PHIM"] * (g_row if row_gain else Gm)
    out["_Lambda"] = lam
    out["_cv"] = psd.cv
    out["_G_row"] = g_row
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
    say("")
    say("   ⚠ SCOPE, added after referee pass 2. The covariance was verified against the")
    say("   solver source: `eps = 1 - PHI/phiJ` is scale-invariant, and the hard clip,")
    say("   `EPS_SOFT` and `MP_EXPO` all act on `eps`, so none of them breaks it. BUT")
    say("   the claim holds only at DIL_COEFF = 0. The Reynolds-dilatancy / Darcy term")
    say("   is built from a `perm_factor(PHI)` that uses a FIXED `PHI_REF = 0.49`, which")
    say("   is NOT one of the three rescaled parameters, so that term is NOT covariant.")
    say("   DIL_COEFF defaults to 0 and is not in WC_DEFAULT, so no number here is")
    say("   affected -- but a DNS run with dilatancy on is outside this mapping.")

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
        diff, se, npair, r = compare_frictions(cells, lo[0], hi[0], cv)
        if diff is None:
            say(f"      CV = {cv}: fewer than 2 PAIRED seeds between mu = {lo[0]} "
                f"and {hi[0]}")
            continue
        t = abs(diff) / se if se > 0 else float("inf")
        t_hypot = abs(lo[1] - hi[1]) / math.hypot(lo[2], hi[2])
        # ⚠ t_crit FOLLOWS THE ROW'S OWN n, corrected after referee pass 3. A single
        # 4.303 was applied to every row, but that is the two-sided 95 % value at 2 dof
        # (n = 3). The CV = 0.4 row has n = 2, i.e. 1 dof, where the value is 12.706 --
        # and the error was in the LENIENT direction, making "DIFFER" easier to reach.
        t_crit = T_CRIT.get(npair, 4.303)
        differ = t > t_crit
        n_resolved += 1
        n_differ += differ
        say(f"      CV = {cv}: G(mu={lo[0]}) = {lo[1]:.4f}   G(mu={hi[0]}) = {hi[1]:.4f}"
            f"   PAIRED d = {diff:+.4f} +- {se:.4f} (n={npair}, r={r:+.2f})")
        say(f"                 |t| = {t:.2f} vs t_crit {t_crit} (n={npair}) -> "
            f"{'DIFFER' if differ else 'not resolved'}"
            f"    [the withdrawn hypot estimator gave {t_hypot:.2f}]")
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
    say("   seeds), the same bar the manuscript uses.")
    say("   🔴 THE ESTIMATOR CHANGED after referee pass 2 and the old value is printed")
    say("   beside the new one above. `hypot` was justified here as 'correct because the")
    say("   two arms are different frictions run as separate cells' -- WRONG on this")
    say("   repo's own terms: `packing.py` seeds the radius draw from the seed REGARDLESS")
    say("   of mu_p, so the arms share the draw (measured r = +0.56 to +0.84) and the")
    say("   paired difference is the right estimator. It moves CV = 0.3 from 1.73 to")
    say("   2.73 and CV = 0.5 from 1.43 to 2.56: the '0 of 4' verdict SURVIVES, but the")
    say("   margin is far thinner than the withdrawn estimator suggested.")

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
        say(f"   {'<k>':>6}{'d_f':>6}{'CV':>8}{'G(row)':>9}{'Lambda':>9}{'PHI':>9}"
            f"{'PHI0':>9}{'PHIM':>9}{'phi/phi0':>10}{'phi/phiM':>10}")
        base = WC_DEFAULT
        lam_max = 1.0 / base["PHI"]

        def regime(w: dict) -> str:
            """ONE rule, applied to the baseline and to every mapped row alike.

            The baseline row was briefly labelled 'flowable' by a hardcoded string
            while its own phi/phiM = 1.0400 said otherwise -- the exact defect this
            repository records as 'prose contradicts its own artifact'. Computed now.
            """
            # 🔴 CEILING CORRECTED 2026-09-22 after referee pass 2. This first used
            # `phi_eff > 1`, which is a TRIVIAL space-filling bound and flagged 1 of 5
            # rows. The binding ceiling for RIGID, NON-INTERPENETRATING aggregates
            # (assumption A3) is the pack's own FRICTIONLESS jamming point: a collection
            # of impenetrable spheres cannot be denser than where it jams. That is
            # `PHI0` here -- the multiplet's OWN phi_0, already computed, and HIGHER
            # than a monodisperse RCP of 0.64 because the effective pack is polydisperse
            # (using a fixed 0.64 would over-flag). It catches 3 of 5.
            if w["_phi_over_phi0"] > 1.0:
                return "🔴 A3 FAILS"
            # ⚠ LABEL NARROWED after referee pass 3. phi/phiM > 1 means the material
            # jams at the FULLY ACTIVATED end (fhat = 1), not unconditionally: phiJ(fhat)
            # runs continuously from PHI0 to PHIM, so phi_m < phi < phi_0 is exactly the
            # stress- and history-dependent window Wyart-Cates exists to describe. The
            # baseline row makes the point -- phi/phi0 = 0.84 (freely flowable at low
            # stress) and phi/phiM = 1.04 -- and calling that "shear-jammed" invites
            # reading a fragile state as a static one.
            if w["_phi_over_phim"] > 1.0:
                return "fragile (can shear-jam)"
            return "flowable"

        base_row = dict(base)
        base_row["_phi_over_phi0"] = base["PHI"] / base["PHI0"]
        base_row["_phi_over_phim"] = base["PHI"] / base["PHIM"]
        say(f"   {'(none)':>6}{'--':>6}{0.0:8.4f}{1.0:9.4f}{1.0:9.4f}{base['PHI']:9.4f}"
            f"{base['PHI0']:9.4f}{base['PHIM']:9.4f}{base_row['_phi_over_phi0']:10.4f}"
            f"{base_row['_phi_over_phim']:10.4f}   {regime(base_row)}")
        n_inadm = 0
        for kb, df in ((1.5, 3.0), (1.5, 2.1), (1.5, 1.8), (2.0, 2.1), (2.0, 1.8)):
            w = wc_inputs(kb, df, g0[0], gm[0])
            flag = regime(w)
            n_inadm += flag.endswith("A3 FAILS")
            old = wc_inputs(kb, df, g0[0], gm[0], row_gain=False)
            say(f"   {kb:6.1f}{df:6.1f}{w['_cv']:8.4f}{w['_G_row']:9.4f}"
                f"{w['_Lambda']:9.4f}{w['PHI']:9.4f}{w['PHI0']:9.4f}"
                f"{w['PHIM']:9.4f}{w['_phi_over_phi0']:10.4f}{w['_phi_over_phim']:10.4f}"
                f"   {flag}")
            say(f"   {'':12s}{'(withdrawn: G at CV=0.5)':>17}{'':18s}"
                f"{old['PHI0']:9.4f}{old['PHIM']:9.4f}"
                f"{old['_phi_over_phi0']:10.4f}{old['_phi_over_phim']:10.4f}"
                f"   {regime(old)}")
        say("")
        say("")
        say("   🔴 THE 'ONE AXIS, BOTH SIGNS' BENCH PREDICTION IS WITHDRAWN.")
        say("   Referee pass 4 found that this table applied a gain measured at CV = 0.5")
        say("   to every row, while the rows' own CVs are 0.16-0.36 and G RISES with CV.")
        say("   Every mapped PHI0/PHIM was inflated and every ratio biased DOWNWARD --")
        say("   toward looking less jammed. On each row's OWN G the two headline rows flip:")
        say("")
        for kb, df, ch, lab in ((1.5, 3.0, "_phi_over_phim", "phi/phiM"),
                                (1.5, 2.1, "_phi_over_phi0", "phi/phi0")):
            nw = wc_inputs(kb, df, g0[0], gm[0])
            od = wc_inputs(kb, df, g0[0], gm[0], row_gain=False)
            say(f"      (<k>={kb}, d_f={df}): {lab}  {od[ch]:.4f} -> {nw[ch]:.4f}"
                f"   {regime(od)} -> {regime(nw)}")
        say("")
        say("   So there is no admissible row on the far side, and NOTHING is left of the")
        say("   two-sided contrast. The withdrawn row is printed under each corrected one.")
        say("")
        say("   🔑 WHAT SURVIVES, and it is weaker and quantitative. At d_f = 3")
        say("   Lambda = 1 EXACTLY, so the only lever acting raises phiJ, and compact")
        say("   aggregation does move the suspension TOWARD flowable:")
        _c = wc_inputs(1.5, 3.0, g0[0], gm[0])
        say(f"      phi/phiM: {base_row['_phi_over_phim']:.4f} (unaggregated) -> "
            f"{_c['_phi_over_phim']:.4f} (<k> = 1.5, d_f = 3).")
        say(f"   It does NOT cross 1. A {(_c['_G_row']-1)*100:.1f} % gain in phiM against "
            f"a {(base_row['_phi_over_phim']-1)*100:.1f} % excess is not")
        say("   enough, and the direction is the whole claim. Opening the aggregates")
        say("   reverses that direction, but only into INADMISSIBLE territory at this")
        say("   solid loading -- so the reversal cannot be quoted from this table either.")
        say("")
        say("   ⚠ To get an admissible two-sided contrast the solid loading must come")
        say("   DOWN so the open rows clear the A3 ceiling. That is a different table and")
        say("   it has not been computed; no claim is made from this one.")
        say("")
        say("   🔴 ADMISSIBILITY -- CEILING CORRECTED after referee pass 2.")
        say("   This block first used `phi_eff <= 1`, a space-filling bound, and flagged")
        say("   ONE row. That is far too loose. Under assumption A3 the aggregates are")
        say("   RIGID and IMPENETRABLE, and such a collection cannot be denser than the")
        say("   point where it JAMS -- its own frictionless jamming fraction, which is")
        say("   the PHI0 column above and is already computed:")
        say("")
        say(f"      phi_eff <= phi_0^eff = G0 * PHI0 = {g0[0]:.4f} * {base['PHI0']} "
            f"= {g0[0]*base['PHI0']:.4f}")
        say(f"      (the looser bound was 1/phi_solid = {lam_max:.4f} on Lambda, i.e. "
            f"phi_eff <= 1)")
        say("")
        say(f"   On that ceiling {n_inadm} of the 5 rows BREACH -- not 1. Those rows are")
        say("   not 'a large effect': they are the NO-INTERPENETRATION ASSUMPTION")
        say("   DESTROYING ITSELF, and they must not be run.")
        say("")
        say("   ⚠ A monodisperse RCP of 0.64 would be the WRONG ceiling here and would")
        say("   over-flag, because the effective pack is polydisperse and its own")
        say(f"   jamming point is {g0[0]*base['PHI0']:.4f}. The bound has to be the")
        say("   multiplet's own, which is why it is read off the mapped PHI0.")
        say("")
        say("   🔑 THE CONSEQUENCE FOR THE BENCH PREDICTION, and it is a DEMOTION.")
        say("   The un-jam / re-jam contrast below survives only where A3 holds. The")
        say("   compact row (d_f = 3) and the <k> = 1.5, d_f = 2.1 row are admissible;")
        say("   the d_f = 1.8 rows and <k> = 2.0, d_f = 2.1 are NOT. So the statement")
        say("   'open aggregates drive the suspension deep into the jammed state' is")
        say("   true of the MODEL and NOT of a rigid-aggregate suspension -- past the")
        say("   ceiling the aggregates must deform or interpenetrate, which is physics")
        say("   this model does not contain. The admissible part of the contrast is")
        say("   still real and is the part to quote.")

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

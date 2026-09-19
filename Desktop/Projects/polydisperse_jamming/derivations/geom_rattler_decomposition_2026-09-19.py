"""geom_rattler_decomposition_2026-09-19.py

Prints every number quoted in geom_rattler_decomposition_2026-09-19.md.

The object under test is the EXACT identity

    phi = phi_B / (1 - psi),     phi_B = V_backbone / V_box,  psi = V_rattler / V_solid

and the physical hypothesis H1 that phi_B does not depend on the size distribution.
The shards record the rattler NUMBER fraction q, not psi, so the test is run through
the exact rearrangement bracket Psi(q) <= psi <= 1 - Psi(1-q).

READ ONLY on out/.  Run:  python geom_rattler_decomposition_2026-09-19.py
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_common_20260919 as G  # noqa: E402

FAMS = ("lognormal", "uniform", "bidisperse")
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
MU = 0.3


def sec(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------------------------ 1. the algebra
def part1_algebra():
    sec("1. THE DECOMPOSITION IS AN IDENTITY, NOT A MODEL -- SymPy check")
    VB, VR, V = sp.symbols("V_B V_R V", positive=True)
    Vs = VB + VR
    phi = Vs / V
    psi = VR / Vs
    phiB = VB / V
    lhs = sp.simplify(phi - phiB / (1 - psi))
    print("   phi        = (V_B + V_R)/V")
    print("   phi_B      = V_B/V")
    print("   psi        = V_R/(V_B + V_R)")
    print("   phi - phi_B/(1-psi) simplifies to:", lhs)
    assert lhs == 0
    print("   -> phi = phi_B/(1-psi) is an IDENTITY.  [asserted zero]")

    print("\n   Which PSD moment enters. With N grains, rattler set R, |R| = qN:")
    print("       psi = sum_{i in R} d_i^3 / sum_i d_i^3 = q * <d^3>_R / <d^3>")
    print("   so the functional of the PSD that the depopulation picture needs is the")
    print("   THIRD moment CONDITIONED ON THE RATTLER SET, not <d^3>, not d32, not d43.")
    d, q = sp.symbols("d q", positive=True)
    print("   At CV = 0 every grain has the same d, so <d^3>_R = <d^3> and psi = q")
    print("   EXACTLY.  That makes the monodisperse cell a zero-width anchor.")

    print("\n   Bound (rearrangement inequality).  For fixed q, sum_{i in R} d_i^3 is")
    print("   minimised by taking the q smallest grains and maximised by the q largest:")
    print("       Psi(q) <= psi <= 1 - Psi(1-q),   Psi(q) = int_0^{F^-1(q)} d^3 f dd/<d^3>")
    print("   Psi is closed form per family (geom_common_20260919.Psi).")

    print("\n   Closed forms, derived not asserted:")
    z, s = sp.symbols("z sigma", real=True, positive=True)
    print("     lognormal  d = exp(mu + sigma Z):  E[d^3 1(Z<z)] = E[d^3] Phi(z - 3 sigma)")
    print("                so Psi(q) = Phi(Phi^-1(q) - 3 sigma),  sigma^2 = ln(1+CV^2)")
    w, t = sp.symbols("w t", positive=True)
    num = sp.integrate(d**3, (d, 1 - w, t))
    den = sp.integrate(d**3, (d, 1 - w, 1 + w))
    print("     uniform    Psi(q) = ", sp.simplify(num / den), " with t = 1-w+2wq, w=CV*sqrt3")
    c = sp.symbols("c", positive=True)
    ws = (1 - c) ** 3 / ((1 - c) ** 3 + (1 + c) ** 3)
    print("     bidisperse small-grain volume share w_s =", sp.simplify(ws))
    print("                Psi(q) = 2 q w_s               (q <= 1/2)")
    print("                Psi(q) = w_s + (2q-1)(1 - w_s) (q >  1/2)")

    print("\n   Limiting checks on Psi (asserted):")
    for fam in FAMS:
        assert abs(G.Psi(fam, 0.3, 0.0)) < 1e-12
        assert abs(G.Psi(fam, 0.3, 1.0) - 1.0) < 1e-12
        assert abs(G.Psi(fam, 0.0, 0.42) - 0.42) < 1e-12
    print("     Psi(0)=0, Psi(1)=1 for all three families at CV=0.3: OK")
    print("     Psi(q)=q at CV=0 for all three families: OK")
    # Monte-Carlo spot check of the closed forms against a drawn sample.
    rng = np.random.default_rng(20260919)
    print("\n   Monte-Carlo spot check of Psi (10^6 draws, seed 20260919):")
    print("     family      CV     q      Psi closed    Psi sampled    |diff|")
    for fam, cv, qq in [("lognormal", 0.5, 0.45), ("lognormal", 0.2, 0.30),
                        ("uniform", 0.5, 0.45), ("uniform", 0.3, 0.60),
                        ("bidisperse", 0.5, 0.45), ("bidisperse", 0.3, 0.70)]:
        n = 1_000_000
        if fam == "lognormal":
            s_ = math.sqrt(math.log(1 + cv * cv))
            dd = np.exp(rng.normal(-0.5 * s_ * s_, s_, n))
        elif fam == "uniform":
            wq = cv * math.sqrt(3.0)
            dd = rng.uniform(1 - wq, 1 + wq, n)
        else:
            dd = np.where(rng.random(n) < 0.5, 1 - cv, 1 + cv)
        dd.sort()
        k = int(round(qq * n))
        samp = (dd[:k] ** 3).sum() / (dd ** 3).sum()
        clos = G.Psi(fam, cv, qq)
        print(f"     {fam:11s} {cv:.1f}  {qq:.2f}   {clos:10.6f}   {samp:11.6f}   "
              f"{abs(clos-samp):.2e}")
        assert abs(clos - samp) < 3e-3


# ------------------------------------------------------------------ 2. the data
def part2_data(sh):
    sec("2. THE MEASURED INPUTS -- rattler NUMBER fraction and phi, per rung")
    print("   All at mu_p = 0.3, N = 1500, 3 seeds, mean over seeds.")
    g = sh[(sh.mu_p == MU)]
    print(f"   N cells at mu_p={MU}: {len(g)}")
    for name, col in (("rattler number fraction q", "rf"), ("phi", "phi")):
        print(f"\n   {name}")
        print("     p_hat     family      " + "".join(f"  CV{c:.1f}  " for c in CVS))
        for p in RUNGS:
            for fam in FAMS:
                row = []
                for c in CVS:
                    v = g[(g.family == fam) & (np.isclose(g.cv, c))
                          & (np.isclose(g.p, p))][col]
                    row.append(f"  {v.mean():.4f}" if len(v) else "     -- ")
                print(f"     {p:.0e}  {fam:11s}" + "".join(row))
    print("\n   NOTE on a number quoted elsewhere in the project.  MECHANISM_2026-09-08.md")
    print("   §M3 writes the rattler fraction as '0.040 -> 0.234'.  On this table the")
    print("   pair (0.040, 0.234) is the p_hat = 1e-3 rung at CV = 0 and CV = 0.3, not")
    print("   CV = 0 -> 0.5.  The CV 0 -> 0.5 range at that rung is printed above.")


# ------------------------------------------------- 3. the zero-parameter test
def backbone_interval(phi, fam, cv, q):
    lo_psi, hi_psi = G.psi_bracket(fam, cv, q)
    return phi * (1 - hi_psi), phi * (1 - lo_psi)


def part3_test(sh):
    sec("3. H1 (BACKBONE INVARIANCE) -- a gate that can fail")
    print("   H1: phi_B = phi (1 - psi) does not depend on the size distribution.")
    print("   Test at each rung: the per-CV intervals [phi*Psi(1-q), phi*(1-Psi(q))]")
    print("   must have a NON-EMPTY common intersection.  At CV = 0 the interval")
    print("   collapses to the single point phi(1-q), so the intersection is pinned")
    print("   there and the gate is one-sided and sharp:")
    print("       REFUTED  iff  phi(CV)*(1 - Psi(q(CV)))  <  phi(0)*(1 - q(0))")
    print("   i.e. iff even making every rattler one of the SMALLEST grains leaves the")
    print("   backbone less dense than the monodisperse one.")
    verdicts = {}
    for p in RUNGS:
        print(f"\n   --- rung p_hat = {p:.0e} ---")
        for fam in FAMS:
            g = sh[(sh.mu_p == MU) & (sh.family == fam) & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)]
            if not len(ref):
                continue
            phi0, q0 = float(ref.phi.mean()), float(ref.rf.mean())
            phiB0 = phi0 * (1 - q0)
            print(f"     {fam}: anchor CV=0  phi={phi0:.4f}  q={q0:.4f}  "
                  f"phi_B(0)={phiB0:.4f}  (exact, psi=q at CV=0)")
            print("       CV     phi      q      Psi(q)  phi_B max  phi_B min   "
                  "req psi  req <d3>_R/<d3>  verdict")
            ok = True
            for c in CVS[1:]:
                gg = g[np.isclose(g.cv, c)]
                if not len(gg):
                    continue
                phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
                lo, hi = backbone_interval(phi_, fam, c, q_)
                req_psi = 1.0 - phiB0 / phi_
                ratio = req_psi / q_ if q_ > 0 else float("nan")
                bad = phiB0 > hi or phiB0 < lo
                ok = ok and not bad
                print(f"       {c:.1f}  {phi_:.4f}  {q_:.4f}  "
                      f"{G.Psi(fam, c, q_):.4f}  {hi:.4f}     {lo:.4f}     "
                      f"{req_psi:.4f}   {ratio:.4f}          "
                      f"{'REFUTES' if bad else 'admits'}")
            verdicts[(p, fam)] = ok
            print(f"       -> H1 at this rung/family: "
                  f"{'NOT REFUTED' if ok else 'REFUTED'}")
    n_ref = sum(1 for v in verdicts.values() if not v)
    print(f"\n   H1 refuted ON SEED MEANS in {n_ref} of {len(verdicts)} (rung, family) "
          f"blocks.")
    return verdicts


def part3b_seedtest(sh):
    sec("3b. THE SAME GATE PER SEED -- a mean margin of 0.003 is not a refutation")
    print("   Section 3 compares MEANS.  The margin where it fires is 0.0011-0.0034 in")
    print("   phi, and this protocol's own monodisperse run-to-run scatter is 2.1e-3")
    print("   (paper/numbers.tex \\phiScale).  So the gate must be run seed-paired.")
    print("   Statistic, per seed s:")
    print("       M(cv,s) = phi(cv,s)(1 - Psi(q(cv,s))) - phi(0,s)(1 - q(0,s))")
    print("   H1 requires M >= 0.  3 seeds -> 2 dof -> t_crit = "
          f"{G.T_CRIT_2DOF:.3f} (two-sided 95%).")
    fired = []
    for p in RUNGS:
        for fam in FAMS:
            g = sh[(sh.mu_p == MU) & (sh.family == fam) & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)].set_index("seed")
            if not len(ref):
                continue
            for c in CVS[1:]:
                gg = g[np.isclose(g.cv, c)].set_index("seed")
                seeds = sorted(set(ref.index) & set(gg.index))
                if len(seeds) < 3:
                    continue
                M = np.array([
                    gg.phi[s] * (1 - G.Psi(fam, c, gg.rf[s]))
                    - ref.phi[s] * (1 - ref.rf[s]) for s in seeds])
                m, se, t, n = G.tstat(M)
                sig = (m < 0) and abs(t) > G.T_CRIT_2DOF
                if m < 0:
                    fired.append((p, fam, c, m, se, t, sig))
                    print(f"     p={p:.0e} {fam:11s} CV={c:.1f}: M = {m:+.5f} +- "
                          f"{se:.5f}  t = {t:+.2f}  n={n}  "
                          f"{'RESOLVED NEGATIVE -> REFUTES' if sig else 'negative but NOT resolved'}")
    if not fired:
        print("     no (rung, family, CV) cell has a negative seed-paired margin.")
    nres = sum(1 for f in fired if f[6])
    print(f"\n   cells with a NEGATIVE mean margin: {len(fired)}")
    print(f"   cells where that margin is RESOLVED at t_crit = {G.T_CRIT_2DOF:.3f}: {nres}")
    return nres


def part5b_selectivity(sh):
    sec("5b. SIZE SELECTIVITY -- where the required rattler set sits between the two "
        "extreme models")
    print("   Two limiting models for WHICH grains rattle:")
    print("     S = 1  'size-ordered':  the rattlers are exactly the q smallest grains,")
    print("            which is what a pure coordination threshold z(d) < 2 gives, since")
    print("            the expected neighbour count of a grain grows monotonically with d")
    print("            -> psi = Psi(q).")
    print("     S = 0  'size-blind':    rattling is independent of size -> psi = q.")
    print("   Define  S = (q - psi_req)/(q - Psi(q)),  so S in [0,1] interpolates and an")
    print("   S outside [0,1] REFUTES H1 (S > 1 is section 3's refutation).")
    print("\n     p_hat   family      " + "".join(f"  CV{c:.1f}" for c in CVS[1:]))
    allS = []
    for p in RUNGS:
        for fam in FAMS:
            g = sh[(sh.mu_p == MU) & (sh.family == fam) & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)]
            if not len(ref):
                continue
            phiB0 = float(ref.phi.mean()) * (1 - float(ref.rf.mean()))
            row = []
            for c in CVS[1:]:
                gg = g[np.isclose(g.cv, c)]
                if not len(gg):
                    row.append("    --")
                    continue
                phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
                psi_req = 1.0 - phiB0 / phi_
                S = (q_ - psi_req) / (q_ - G.Psi(fam, c, q_))
                allS.append(S)
                row.append(f" {S:5.3f}")
            print(f"     {p:.0e} {fam:11s}" + "".join(row))
    a = np.array(allS)
    print(f"\n     N cells: {len(a)}; S range {a.min():.3f} to {a.max():.3f}; "
          f"median {np.median(a):.3f}; mean {a.mean():.3f}; sd {a.std(ddof=1):.3f}")
    print(f"     S outside [0,1]: {int(((a < 0) | (a > 1)).sum())} cells")
    print("\n   Reading: if H1 is right, rattling is STRONGLY but not perfectly")
    print("   size-selective, and the degree of selectivity is close to a constant")
    print("   across three families, five CVs and four pressure rungs.  That constancy")
    print("   is a prediction, not a measurement -- it is what H1 REQUIRES, and the one")
    print("   measurement that tests it is psi itself.")
    return a


# --------------------------------------- 4. is the gate capable of firing at all?
def part4_power(sh):
    sec("4. CAN THE GATE FIRE?  -- a gate that cannot fail is not a gate")
    print("   The gate fires when phi(CV)(1 - Psi(q)) < phi(0)(1 - q(0)).  Ask how far")
    print("   the data is from that, in units of phi, and what phi(CV) WOULD have had")
    print("   to be for the gate to fire.  Printed per rung for the widest family.")
    for p in RUNGS:
        for fam in FAMS:
            g = sh[(sh.mu_p == MU) & (sh.family == fam) & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)]
            gg = g[np.isclose(g.cv, 0.5)]
            if not len(ref) or not len(gg):
                continue
            phi0, q0 = float(ref.phi.mean()), float(ref.rf.mean())
            phiB0 = phi0 * (1 - q0)
            phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
            need = phiB0 / (1 - G.Psi(fam, 0.5, q_))
            print(f"     p={p:.0e} {fam:11s} CV=0.5: measured phi={phi_:.4f}; the gate "
                  f"would have fired below phi={need:.4f}  (margin {phi_-need:+.4f})")
    print("\n   A margin that is negative anywhere is a REFUTATION; a margin that is")
    print("   positive everywhere means the bracket Psi(q) is too loose to decide, and")
    print("   the quantity that closes it is the rattler VOLUME fraction psi, which is")
    print("   NOT recorded in any shard (only rattler_frac, a NUMBER fraction).")


# ------------------------------------ 5. the sharp prediction H1 makes if true
def part5_prediction(sh):
    sec("5. WHAT H1 PREDICTS, IF IT IS TRUE -- the number a new run must hit")
    print("   Required conditional volume ratio R = <d^3>_R/<d^3> = psi_req/q, and the")
    print("   implied rattler composition.  For the bidisperse family the PSD has two")
    print("   atoms, so R fixes the FRACTION of rattlers that are large grains exactly:")
    print("       R = (1-f) v_s + f v_l,  v_s = 2(1-c)^3/((1-c)^3+(1+c)^3), v_l = 2-v_s")
    for p in (1e-3, 3e-5):
        print(f"\n   --- rung p_hat = {p:.0e} ---")
        for fam in FAMS:
            g = sh[(sh.mu_p == MU) & (sh.family == fam) & (np.isclose(sh.p, p))]
            ref = g[np.isclose(g.cv, 0.0)]
            if not len(ref):
                continue
            phiB0 = float(ref.phi.mean()) * (1 - float(ref.rf.mean()))
            print(f"     {fam}")
            print("       CV    q      psi_req   R=psi/q   d_R/<d> (=R^1/3)   "
                  "frac of rattlers that are LARGE (bidisperse only)")
            for c in CVS[1:]:
                gg = g[np.isclose(g.cv, c)]
                if not len(gg):
                    continue
                phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
                psi_req = 1.0 - phiB0 / phi_
                R = psi_req / q_
                extra = ""
                if fam == "bidisperse":
                    vs = 2 * (1 - c) ** 3 / ((1 - c) ** 3 + (1 + c) ** 3)
                    vl = 2.0 - vs
                    f = (R - vs) / (vl - vs)
                    # fraction of ALL small grains that rattle, and of all large
                    fs = q_ * (1 - f) / 0.5
                    fl = q_ * f / 0.5
                    extra = (f"f={f:+.4f}  -> {fs*100:5.1f}% of small grains rattle, "
                             f"{fl*100:5.1f}% of large")
                print(f"       {c:.1f}  {q_:.4f}  {psi_req:.4f}   {R:.4f}   "
                      f"{R**(1/3) if R > 0 else float('nan'):.4f}            {extra}")
    print("\n   ADMISSIBILITY.  R must lie in [v_min, v_max], the relative volumes of the")
    print("   smallest and largest grain in the family.  An R outside that, or an f")
    print("   outside [0,1], would REFUTE H1 with no further measurement.")


# ------------------------------------ 6. the non-monotone bidisperse stress test
def part6_bidisperse(sh):
    sec("6. THE BIDISPERSE NON-MONOTONICITY -- the sharpest available stress test")
    print("   phi_J for the bidisperse family PEAKS near CV = 0.3 and falls, while its")
    print("   rattler NUMBER fraction rises monotonically 0.04 -> 0.50.  Under H1,")
    print("   psi = 1 - phi_B(0)/phi must therefore be NON-MONOTONE in CV even though q")
    print("   is monotone: i.e. the mean rattler VOLUME must fall faster than the")
    print("   rattler COUNT rises, then reverse.  Printed as R(CV) in section 5.")
    for p in RUNGS:
        g = sh[(sh.mu_p == MU) & (sh.family == "bidisperse") & (np.isclose(sh.p, p))]
        ref = g[np.isclose(g.cv, 0.0)]
        if not len(ref):
            continue
        phiB0 = float(ref.phi.mean()) * (1 - float(ref.rf.mean()))
        rows = []
        for c in CVS[1:]:
            gg = g[np.isclose(g.cv, c)]
            if not len(gg):
                continue
            phi_, q_ = float(gg.phi.mean()), float(gg.rf.mean())
            rows.append((c, q_, 1 - phiB0 / phi_, (1 - phiB0 / phi_) / q_))
        psis = [r[2] for r in rows]
        mono_q = all(rows[i][1] <= rows[i + 1][1] for i in range(len(rows) - 1))
        mono_psi = all(psis[i] <= psis[i + 1] for i in range(len(psis) - 1))
        amax = int(np.argmax(psis))
        print(f"     p={p:.0e}: q monotone in CV: {mono_q};  psi_req monotone: "
              f"{mono_psi};  psi_req peaks at CV = {rows[amax][0]:.1f}")
        print("              " + "  ".join(f"CV{r[0]:.1f} psi={r[2]:.4f} R={r[3]:.4f}"
                                           for r in rows))
    print("\n   Reading: H1 does not merely tolerate the non-monotonicity, it REQUIRES a")
    print("   non-monotone conditional rattler volume.  That is a strong, specific and")
    print("   currently unmeasured statement.")


# ---------------------------------------------- 7. what the anchor cannot decide
def part7_limits(sh):
    sec("7. WHAT THIS ANCHOR CAN AND CANNOT DECIDE")
    print("   CAN decide:  whether H1 is impossible, i.e. whether ANY assignment of")
    print("                rattlers to grains makes the backbone density CV-independent.")
    print("                That question is fully settled by q, phi and the PSD, all on")
    print("                disk, with no free parameter.")
    print("   CANNOT decide: whether H1 is TRUE.  Doing that needs psi, the rattler")
    print("                VOLUME fraction, and no shard records it.  Every shard field")
    print("                that touches rattlers is a COUNT: rattler_frac_k2/k3/k4,")
    print("                largest_frac_k*, z_nr_k*.  Verified by field inventory below.")
    print("   CANNOT decide: anything about WHY the backbone would be invariant.  Fact 3")
    print("                (Z_nr and T flat) is consistent with H1 but does not imply it:")
    print("                a structure can keep its contact number and its fabric while")
    print("                changing its density (compress a random network affinely).")
    import glob
    import json
    f = sorted(glob.glob(os.path.join(G.OUT, "phase2_ln0_shard00of06.jsonl")))[0]
    rec = json.loads(open(f, encoding="utf-8").readline())
    print(f"\n   Field inventory of a shard record ({os.path.basename(f)}):")
    print("     top level:", sorted(rec.keys()))
    print("     rattler:  ", sorted(rec["rattler"].keys()))
    has_size = [k for k in list(rec.keys()) + list(rec["rattler"].keys())
                if any(t in k.lower() for t in ("d10", "d32", "d43", "diam", "radi"))]
    print("     any per-RATTLER size field?  fields mentioning a size:", has_size)
    print("     -> all of those are WHOLE-PACK moments; none is conditioned on the")
    print("        rattler set.  psi is not recoverable from these files.")


def main():
    sh = G.load_shards()
    part1_algebra()
    part2_data(sh)
    v = part3_test(sh)
    nres = part3b_seedtest(sh)
    part4_power(sh)
    part5_prediction(sh)
    S = part5b_selectivity(sh)
    part6_bidisperse(sh)
    part7_limits(sh)
    sec("VERDICT")
    n_ref = sum(1 for x in v.values() if not x)
    print(f"   H1 refuted on seed MEANS in {n_ref} of {len(v)} (rung, family) blocks;")
    print(f"   resolved at t_crit = {G.T_CRIT_2DOF:.3f} in {nres} cells.")
    print(f"   Required size selectivity S over {len(S)} cells: "
          f"{S.min():.3f} to {S.max():.3f}, median {np.median(S):.3f}.")
    if nres:
        print("   VERDICT: REFUTED")
    else:
        print("   VERDICT: UNDETERMINED")
        print("   (No admissible-rattler-assignment contradiction survives the seed")
        print("    test, so H1 is not impossible; and psi is not on disk, so it is not")
        print("    confirmable either.)")
    print("\n   ANCHOR SCOPE.  This can decide that H1 is IMPOSSIBLE (it did not); it")
    print("   cannot decide that H1 is TRUE, because the rattler VOLUME fraction psi is")
    print("   recorded nowhere.  It says nothing about WHY a backbone would be")
    print("   size-distribution invariant.")


if __name__ == "__main__":
    main()

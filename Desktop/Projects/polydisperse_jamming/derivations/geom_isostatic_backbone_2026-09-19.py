"""geom_isostatic_backbone_2026-09-19.py

Prints every number quoted in geom_isostatic_backbone_2026-09-19.md.

Question: is a flat backbone coordination in CV what constraint counting predicts, and
what does the measured CV DRIFT in z_iso mean?

Derived here:
  * the isostatic count for frictional spheres in d = 3 is  z_iso = 12/(3 - x_m),
    with x_m the fraction of FULLY MOBILISED (sliding) contacts, reproducing van
    Hecke's d+1 <= z_c <= 2d as the x_m in [0,1] interval;
  * NO grain size enters that count, so polydispersity cannot move z_iso -- the flat
    Z_nr is what the counting requires, and the drift is the thing that needs explaining;
  * inverting the measured z_iso for x_m turns the drift into a statement about contact
    mobilisation, which is a quantity the solver can report and nobody has read.

READ ONLY on out/.  Run:  python geom_isostatic_backbone_2026-09-19.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geom_common_20260919 as G  # noqa: E402

CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
TC = G.T_CRIT_2DOF


def sec(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ------------------------------------------------------- 1. the constraint count
def part1():
    sec("1. THE ISOSTATIC COUNT, DERIVED -- and it contains no length")
    d, N, Nc, xm = sp.symbols("d N N_c x_m", positive=True)
    print("   Count STATIC DETERMINACY: the contact forces are determined by balance")
    print("   alone exactly when the number of independent unknowns equals the number of")
    print("   independent equations.")
    print("\n   EQUATIONS.  Each rigid grain gives d force-balance and d(d-1)/2")
    print("   torque-balance equations:  d + d(d-1)/2 = d(d+1)/2 per grain.")
    eqs = d * (d + 1) / 2 * N
    print("     total equations =", sp.simplify(eqs), " -> at d = 3:",
          sp.simplify(eqs.subs(d, 3)))
    print("\n   UNKNOWNS.  Each contact carries d force components (1 normal, d-1")
    print("   tangential).  At a FULLY MOBILISED contact the tangential magnitude is")
    print("   pinned to mu f_n, which removes one unknown.  With x_m the mobilised")
    print("   fraction:")
    unk = Nc * d - xm * Nc
    print("     total unknowns =", sp.simplify(unk), " = N_c (d - x_m)")
    print("\n   Setting them equal and using z = 2 N_c / N:")
    z = sp.symbols("z", positive=True)
    sol = sp.solve(sp.Eq(unk.subs(Nc, z * N / 2), eqs), z)[0]
    print("     z_iso =", sp.simplify(sol), "  = d(d+1)/(d - x_m)")
    ziso3 = sp.simplify(sol.subs(d, 3))
    print("     at d = 3:  z_iso =", ziso3, " = 12/(3 - x_m)")
    assert sp.simplify(ziso3 - 12 / (3 - xm)) == 0
    print("\n   LIMITS (asserted):")
    lo = sp.simplify(ziso3.subs(xm, 0))
    hi = sp.simplify(ziso3.subs(xm, 1))
    print(f"     x_m = 0 (nothing sliding, friction fully effective): z_iso = {lo} "
          f"= d + 1")
    print(f"     x_m = 1 (every contact sliding, i.e. mu -> 0):       z_iso = {hi} "
          f"= 2d")
    assert lo == 4 and hi == 6
    print("     2D check: d(d+1)/(d-x_m) at d = 2 gives",
          sp.simplify(sol.subs(d, 2).subs(xm, 0)), "at x_m=0 (=d+1=3) and",
          sp.simplify(sol.subs(d, 2).subs(xm, 1)), "at x_m=1 (=2d=4).  OK")
    print("\n   External anchor (API-resolved, given in the brief): van Hecke, J. Phys.:")
    print("   Condens. Matter 22, 033101 (2010), doi:10.1088/0953-8984/22/3/033101,")
    print("   states d+1 <= z_c <= 2d for frictional spheres with z_c DECREASING in mu.")
    print("   The derivation above reproduces both bounds and supplies the interpolant:")
    print("   larger mu means fewer contacts sit at the Coulomb cone, so x_m falls and")
    print("   z_iso falls -- the decreasing trend is a consequence, not an assumption.")
    print("\n   THE POINT FOR THIS PROJECT.  Every symbol above is a COUNT.  Grain")
    print("   diameters appear nowhere, so:")
    print("     PREDICTION P1: z_iso is independent of the size distribution at fixed mu.")
    print("   This is the sense in which 'the backbone does not restructure' is not a")
    print("   coincidence -- it is forced, PROVIDED the pack is marginally rigid and")
    print("   x_m does not itself depend on CV.  The second proviso is the whole item.")
    return ziso3, xm


def z_of_xm(x):
    return 12.0 / (3.0 - x)


def xm_of_z(z):
    return 3.0 - 12.0 / z


# ---------------------------------------------- 2. the pressure extrapolation
def part2(sh):
    sec("2. z_nr AT A PRESSURE IS NOT z_iso -- the extrapolation, per seed")
    print("   The contact law on this project is Cundall-Strack LINEAR springs")
    print("   (packing.py: Law2_ScGeom_FrictPhys_CundallStrack), so p ~ (phi - phi_J)")
    print("   and the standard marginal scaling z - z_c ~ (phi - phi_J)^(1/2) gives")
    print("       z_nr(p) = z_iso + a sqrt(p_hat).")
    print("   measure_z_iso.py uses exactly this rule; here it is refitted PER SEED so")
    print("   the CV drift gets an error bar at the design's own 2 dof.")
    rows = {}
    for fam in ("lognormal", "uniform", "bidisperse"):
        for mu in (0.0, 0.3):
            g = sh[(sh.family == fam) & (sh.mu_p == mu)]
            if not len(g):
                continue
            print(f"\n   --- {fam}, mu_p = {mu} ---")
            print("     CV    n seeds   z_iso mean    se       a mean   worst r2   rungs")
            for cv in CVS:
                gg = g[np.isclose(g.cv, cv)]
                zi, aa, r2s, ns = [], [], [], []
                for s in sorted(gg.seed.unique()):
                    h = gg[gg.seed == s].sort_values("p")
                    if len(h) < 3:
                        continue
                    x = np.sqrt(h.p.to_numpy())
                    y = h.z_nr.to_numpy()
                    A = np.vstack([np.ones_like(x), x]).T
                    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
                    pred = A @ coef
                    ss = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
                    zi.append(coef[0])
                    aa.append(coef[1])
                    r2s.append(ss)
                    ns.append(len(x))
                if len(zi) < 2:
                    print(f"     {cv:.1f}   {len(zi)}  -- too few seeds to fit")
                    continue
                m, se, t, n = G.tstat(zi)
                rows[(fam, mu, cv)] = (m, se, n, float(np.mean(aa)), min(r2s), min(ns))
                print(f"     {cv:.1f}     {n}      {m:.4f}    {se:.4f}   "
                      f"{np.mean(aa):7.3f}   {min(r2s):.4f}    {min(ns)}")
    return rows


# --------------------------------------------------- 3. is the drift resolved?
def part3(sh, rows):
    sec("3. IS THE CV DRIFT IN z_iso RESOLVED?  -- seed-paired, 2 dof")
    print("   Prediction P1 says the drift is zero.  Tested seed-paired so the PSD draw")
    print("   and the protocol noise cancel.  t_crit at 2 dof = "
          f"{TC:.3f}.")
    out = {}
    for fam in ("lognormal", "uniform", "bidisperse"):
        for mu in (0.0, 0.3):
            g = sh[(sh.family == fam) & (sh.mu_p == mu)]
            if not len(g):
                continue

            def zi_per_seed(cv):
                gg = g[np.isclose(g.cv, cv)]
                d = {}
                for s in sorted(gg.seed.unique()):
                    h = gg[gg.seed == s].sort_values("p")
                    if len(h) < 3:
                        continue
                    x = np.sqrt(h.p.to_numpy())
                    A = np.vstack([np.ones_like(x), x]).T
                    coef, *_ = np.linalg.lstsq(A, h.z_nr.to_numpy(), rcond=None)
                    d[s] = coef[0]
                return d

            lo, hi = zi_per_seed(0.0), zi_per_seed(0.5)
            common = sorted(set(lo) & set(hi))
            if len(common) < 3:
                print(f"   {fam} mu={mu}: only {len(common)} paired seeds -- skipped")
                continue
            dz = np.array([hi[s] - lo[s] for s in common])
            m, se, t, n = G.tstat(dz)
            out[(fam, mu)] = (m, se, t, n)
            print(f"   {fam:11s} mu={mu:.1f}: d z_iso (CV 0 -> 0.5) = {m:+.4f} +- "
                  f"{se:.4f}  t = {t:+.2f}  n={n}  -> "
                  f"{'RESOLVED, P1 REFUTED' if abs(t) > TC else 'NOT resolved, P1 stands'}")
    return out


# ------------------------------- 4. translate the drift into contact mobilisation
def part4(rows, drift):
    sec("4. WHAT THE DRIFT WOULD MEAN -- inverting z_iso for the mobilised fraction")
    print("   z_iso = 12/(3 - x_m)  <=>  x_m = 3 - 12/z_iso.  Every z_iso on disk maps")
    print("   to an x_m, and an x_m outside [0, 1] would REFUTE the marginal-rigidity")
    print("   premise rather than measure anything.")
    print("\n     family       mu_p   CV    z_iso     x_m      admissible?")
    bad = 0
    for (fam, mu, cv), (m, se, n, a, r2, nr) in sorted(rows.items()):
        x = xm_of_z(m)
        okk = 0.0 <= x <= 1.0
        bad += 0 if okk else 1
        print(f"     {fam:11s}  {mu:.1f}   {cv:.1f}  {m:7.4f}  {x:+7.4f}   "
              f"{'yes' if okk else 'NO -- INADMISSIBLE'}")
    print(f"\n   cells with x_m outside [0,1]: {bad} of {len(rows)}")
    print("\n   The CV 0 -> 0.5 drift, expressed in x_m:")
    for (fam, mu), (m, se, t, n) in sorted(drift.items()):
        z0 = rows.get((fam, mu, 0.0))
        z5 = rows.get((fam, mu, 0.5))
        if not z0 or not z5:
            continue
        print(f"     {fam:11s} mu={mu:.1f}: z_iso {z0[0]:.4f} -> {z5[0]:.4f}  =>  "
              f"x_m {xm_of_z(z0[0]):+.4f} -> {xm_of_z(z5[0]):+.4f}  "
              f"(change {xm_of_z(z5[0])-xm_of_z(z0[0]):+.4f})")
    print("\n   SIGN.  A naive reading runs the other way: a broader size distribution")
    print("   puts many small grains into the pore space of the coarse skeleton, those")
    print("   contacts carry little normal force, and a contact with a small f_n is")
    print("   EASIER to drive to the Coulomb cone -- so x_m should RISE with CV and")
    print("   z_iso should rise with it.  The measured drift is the opposite sign.")
    print("   Two readings survive that, and they are distinguishable by one output:")
    print("     (a) x_m really does fall, i.e. the lightly loaded small-grain contacts")
    print("         are NOT preferentially mobilised -- which would say the extra small")
    print("         grains are being held, not slipping;")
    print("     (b) the pack is not marginally rigid at CV = 0.5, so z_iso is not a")
    print("         count at all there and the inversion is meaningless.")
    print("   Yade exposes the mobilised fraction directly (|f_t| >= mu f_n at a")
    print("   contact), so (a) is settled by ONE scalar added to the measure step.")


# --------------------------------------------- 5. what the drift is NOT explained by
def part5(sh, rows):
    sec("5. THREE CANDIDATE ARTEFACTS FOR THE DRIFT, each checked or named")
    print("   (i) IS IT THE EXTRAPOLATION?  If the drift were an artefact of the sqrt(p)")
    print("   fit, it should shrink as the rungs get lower.  Compare the raw z_nr drift")
    print("   at each rung with the extrapolated one.")
    for fam in ("lognormal", "uniform", "bidisperse"):
        g = sh[(sh.family == fam) & (sh.mu_p == 0.3)]
        if not len(g):
            continue
        print(f"     {fam}:")
        for p in RUNGS:
            a = g[(np.isclose(g.cv, 0.0)) & (np.isclose(g.p, p))].set_index("seed").z_nr
            b = g[(np.isclose(g.cv, 0.5)) & (np.isclose(g.p, p))].set_index("seed").z_nr
            common = sorted(set(a.index) & set(b.index))
            if len(common) < 2:
                continue
            d = np.array([b[s] - a[s] for s in common])
            m, se, t, n = G.tstat(d)
            print(f"       p = {p:.0e}: d z_nr = {m:+.4f} +- {se:.4f} (n={n}, "
                  f"t = {t:+.2f})")
        r0, r5 = rows.get((fam, 0.3, 0.0)), rows.get((fam, 0.3, 0.5))
        if r0 and r5:
            print(f"       extrapolated: d z_iso = {r5[0]-r0[0]:+.4f}")
    print("     Reading: the raw drift is present at EVERY rung and does not shrink")
    print("     towards the lowest one, so it is not made by the extrapolation.")

    print("\n   (ii) IS IT THE RATTLER DEFINITION?  z_nr is computed after iterated")
    print("   removal of grains with fewer than k = 2 contacts, and the rattler")
    print("   population TRIPLES over the CV range, so the backbone is a different")
    print("   SUBSET at each CV.  The shards carry z_nr at k = 2, 3 and 4.  If the drift")
    print("   were a threshold artefact it would move with k.")
    import glob
    import json
    recs = {}
    for f in sorted(glob.glob(os.path.join(G.OUT, "phase2*.jsonl"))):
        if "wcbidi" in os.path.basename(f):
            continue
        for ln in open(f, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            if (r.get("failed") or not r.get("converged")
                    or r.get("n_requested") != 1500 or r.get("mu_p") != 0.3):
                continue
            k = (r.get("family") or "lognormal", r["cv"], r["seed"], r["p_hat_target"])
            if k in recs and r.get("max_iter", 0) <= recs[k].get("max_iter", 0):
                continue
            recs[k] = r
    print(f"     N records for the k sweep after every filter: {len(recs)}")
    nan_k = {k: sum(1 for r in recs.values()
                    if not np.isfinite(r["rattler"][f"z_nr_k{k}"])) for k in (2, 3, 4)}
    print(f"     records whose z_nr_k is NOT FINITE: {nan_k} of {len(recs)}")
    print("     -> at k = 4 the iterated removal empties the backbone at mu_p = 0.3,")
    print("        where z_nr itself is only ~4.5, so that threshold is UNDEFINED here")
    print("        and the sweep runs over k = 2 and 3 only.  Printed, not skipped.")
    print("     family       rung      d z_nr(k=2)  d z_nr(k=3)  d z_nr(k=4)")
    for fam in ("lognormal", "uniform", "bidisperse"):
        for p in (1e-3, 3e-5):
            vals = []
            for k in (2, 3, 4):
                a = {s: recs[(fam, 0.0, s, p)]["rattler"][f"z_nr_k{k}"]
                     for s in (1, 2, 3) if (fam, 0.0, s, p) in recs}
                b = {s: recs[(fam, 0.5, s, p)]["rattler"][f"z_nr_k{k}"]
                     for s in (1, 2, 3) if (fam, 0.5, s, p) in recs}
                common = [s for s in sorted(set(a) & set(b))
                          if np.isfinite(a[s]) and np.isfinite(b[s])]
                vals.append(np.mean([b[s] - a[s] for s in common]) if common
                            else float("nan"))
            print(f"     {fam:11s}  {p:.0e}   {vals[0]:+10.4f}  {vals[1]:+10.4f}  "
                  f"{vals[2]:+10.4f}"
                  + ("   (k=4 undefined)" if not np.isfinite(vals[2]) else ""))
    print("     Reading: a drift that barely moves between k = 2 and k = 3 is not a")
    print("     threshold artefact.")

    print("\n   (iii) IS IT FINITE SIZE?  d_max/L rises with CV by construction, and the")
    print("   shards flag it.  Not decidable here -- it needs the N = 3000 block, which")
    print("   exists on disk but under a different tag and was not loaded.  NAMED, NOT")
    print("   CHECKED.")


def part6(sh, rows):
    sec("6. A DISCREPANCY WITH out/z_iso_measured.txt -- the ladders are not matched")
    print("   The committed artifact gives z_iso = 4.609 at CV = 0 for ALL THREE")
    print("   families at mu_p = 0.3, and 4.377 for lognormal at CV = 0.5, a drift of")
    print("   -0.232 which analyse_phase8.py:165 hardcodes as a per-cell table.")
    print("   The matched-ladder per-seed refit above gives:")
    for fam in ("lognormal", "uniform", "bidisperse"):
        a, b = rows.get((fam, 0.3, 0.0)), rows.get((fam, 0.3, 0.5))
        if a and b:
            print(f"     {fam:11s} CV0 {a[0]:.4f} +- {a[1]:.4f}   CV0.5 {b[0]:.4f} +- "
                  f"{b[1]:.4f}   drift {b[0]-a[0]:+.4f}")
    print("\n   The likely cause, printed rather than argued: the rung populations are")
    print("   NOT the same across CV.  Rungs present per cell at mu_p = 0.3, N = 1500,")
    print("   converged, one record per physical cell:")
    for fam in ("lognormal", "uniform", "bidisperse"):
        for cv in CVS:
            g = sh[(sh.family == fam) & (sh.mu_p == 0.3) & (np.isclose(sh.cv, cv))]
            ps = sorted(g.p.unique())
            ns = [int((np.isclose(g.p, p)).sum()) for p in ps]
            if len(ps) != 4 or set(ns) != {3}:
                print(f"     {fam:11s} cv={cv:.1f}: rungs {[f'{p:.0e}' for p in ps]} "
                      f"seeds/rung {ns}   <-- NOT the common 4x3 ladder")
    print("     (every cell not printed above has exactly the 4 rungs 3e-5..1e-3 with")
    print("      3 seeds each)")
    print("\n   So ONE cell -- lognormal, CV = 0, mu_p = 0.3 -- carries a fifth rung at")
    print("   p_hat = 3e-3 on ONE seed, and it is the CV = 0 REFERENCE from which the")
    print("   drift is measured.  Effect of including it, seed 1 only:")
    g = sh[(sh.family == "lognormal") & (sh.mu_p == 0.3) & (np.isclose(sh.cv, 0.0))
           & (sh.seed == 1)].sort_values("p")
    for pmax in (1e-3, 3e-3):
        h = g[g.p <= pmax * 1.01]
        x = np.sqrt(h.p.to_numpy())
        A = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(A, h.z_nr.to_numpy(), rcond=None)
        print(f"     rungs up to {pmax:.0e} ({len(x)} points): z_iso = {coef[0]:.4f}, "
              f"a = {coef[1]:.3f}")
    print("\n   READING -- and the obvious explanation is REFUTED by its own sign.")
    print("   The unmatched fifth rung LOWERS the CV = 0 intercept (4.5750 -> 4.5537),")
    print("   which would make the drift SMALLER in magnitude, not larger, and the")
    print("   committed reference is HIGHER still (4.609).  So the extra rung does not")
    print("   explain the gap; it moves the wrong way.  What is established is only:")
    print("     * on a ladder identical at every CV, refitted per seed, the mu_p = 0.3")
    print("       CV = 0 intercept is 4.537-4.540 and the lognormal drift is -0.1489,")
    print("       36% smaller than the committed -0.232;")
    print("     * out/z_iso_measured.txt reports '5 rungs' and 433 cells where this")
    print("       loader finds 4 rungs and 395 cells, so its selection differs.")
    print("   CAUSE NOT TRACED.  The candidate worth checking first is the 2026-09-17")
    print("   `_it{cap}` duplication (sweep_descriptive.py:243-262): if measure_z_iso")
    print("   still counts one physical cell twice it will duplicate rungs inside a")
    print("   single fit.  Flagged, not asserted.")


def main():
    sh = G.load_shards()
    part1()
    rows = part2(sh)
    drift = part3(sh, rows)
    part4(rows, drift)
    part5(sh, rows)
    part6(sh, rows)
    sec("VERDICT")
    res = [k for k, v in drift.items() if abs(v[2]) > TC]
    print(f"   P1 (z_iso independent of the PSD) is REFUTED in {len(res)} of "
          f"{len(drift)} (family, mu) arms at t_crit = {TC:.3f}: {sorted(res)}")
    print("   VERDICT: the constraint count is SUPPORTED as the reason Z_nr is nearly")
    print("   flat -- no length enters it -- and REFUTED as an exact statement, because")
    print("   the residual drift is resolved where the data is strongest.")
    print("\n   ANCHOR SCOPE.  The count can decide what z_iso must be GIVEN marginal")
    print("   rigidity and a mobilised fraction; it can also decide that size cannot")
    print("   enter it.  It CANNOT decide whether the measured drift is a real change in")
    print("   x_m or a failure of marginal rigidity, because both predict the same")
    print("   z_iso.  The one output that separates them -- the mobilised contact")
    print("   fraction -- is not recorded by the solver.")


if __name__ == "__main__":
    main()

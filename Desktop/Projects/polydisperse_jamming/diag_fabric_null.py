"""diag_fabric_null.py -- size the fabric campaign BEFORE it is registered or run.

Three things, all with zero DEM cells, written so PREREG_FABRIC_2026-09-18.md can quote
printed numbers instead of asserted ones.

1. THE FINITE-SAMPLE FLOOR HAS A CLOSED FORM, and it is the reason the raw channel is
   unusable. For `n` INDEPENDENT unit vectors uniform on the sphere,

       F_ij = (1/n) sum_a n_i^a n_j^a ,   dev = F - I/3 ,   A^2 = 1.5 * sum_ij dev_ij^2

   has  E[A^2] = 1/n  EXACTLY. Derivation, every step:

     E[n_i n_j] = delta_ij / 3           (isotropy)
     E[n_i^4]   = 1/5 ,  E[n_i^2 n_j^2] = 1/15  (i != j)
       -- from <n_i^2> = 1/3 and the standard moments of a uniform direction.
     Var(F_ii)      = (1/5 - 1/9)/n = (4/45)/n
     Cov(F_ii,F_jj) = (1/15 - 1/9)/n = -(2/45)/n            (i != j)
     Var(F_ij)      = (1/15 - 0)/n  = (1/15)/n              (i != j)
     E[sum_i dev_ii^2]     = 3 * (4/45)/n = (12/45)/n = (4/15)/n
     E[sum_{i!=j} dev_ij^2] = 6 * (1/15)/n = (6/15)/n
     E[sum_ij dev_ij^2]    = (10/15)/n = (2/3)/n
     E[A^2] = 1.5 * (2/3)/n = 1/n .                                        [QED]

   Consistency check built into the algebra: Var(sum_i F_ii) = 3*(4/45)/n +
   6*(-2/45)/n = 0, as it must be, because sum_i F_ii = 1 identically for unit normals.

   So the PRE-REGISTERED observable is

       T = n_contacts * A^2 ,     E[T] = 1 under the independent-normal null,

   which is free of the contact-count confound by construction. T > 1 is orientational
   CORRELATION in the contact network (an effective sample size n_eff = n / T).

2. THE CONFOUND, measured on committed shards rather than argued. `n_contacts` falls
   with CV because polydispersity depopulates the network, so the FLOOR alone drifts
   upward with CV and a raw `aniso` vs CV slope reads that drift as physics.

3. THE DESIGN. Resolving power of the slope of T against CV at S seeds, against the
   pre-registered effect size. Printed BEFORE the prereg is written, because this
   project's own rule is that a bar whose resolving power is unknown decides nothing.

Writes out/diag_fabric_null.txt. Reads committed shards. Runs no cell.
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np

from measure import fabric_tensor

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

MC_SEED = 20260918
MC_REPS = 4000
NULL_NS = (500, 1000, 2000, 3000, 4000, 6000)

# The production grid this is sizing.
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
MU = 0.3
N_BODIES = 1500
P_PRIMARY = 1e-3
P_SECOND = 1e-4

# Pre-registered effect size, in T units, stated as a CHANGE ACROSS THE CV RANGE:
# T = 1 is the value an orientationally structureless network would return, so
# DELTA_T_TARGET = 0.5 is "half of the whole isotropic floor's worth of extra
# correlation, appearing over CV 0 -> 0.5". It is a declared convention, not a
# measurement, and the verdict grammar reports the ACHIEVED resolving power beside it.
DELTA_T_TARGET = 0.5


def mc_null(n: int, reps: int = MC_REPS, seed: int = MC_SEED):
    """Monte-Carlo the null for `n` independent uniform normals: (mean A, sd A, mean T)."""
    rng = np.random.default_rng(seed + n)
    v = rng.normal(size=(reps, n, 3))
    v /= np.linalg.norm(v, axis=2, keepdims=True)
    F = np.einsum("rni,rnj->rij", v, v) / n
    dev = F - np.eye(3) / 3.0
    a2 = 1.5 * np.einsum("rij,rij->r", dev, dev)
    a = np.sqrt(a2)
    # The vectorised form above must agree with the SHIPPED `fabric_tensor`, which is
    # what produced every number in the campaign -- a null computed by a private copy of
    # the formula certifies the copy (lessons.md 2026-09-05). Checked on a subsample.
    worst = max(abs(fabric_tensor(v[r])["aniso"] - a[r]) for r in range(min(50, reps)))
    assert worst < 1e-12, f"MC null disagrees with measure.fabric_tensor by {worst:g}"
    return float(a.mean()), float(a.std(ddof=1)), float((n * a2).mean()), float((n * a2).std(ddof=1))


def load_shards():
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT, "phase2_*.jsonl"))):
        for ln in open(f, encoding="utf-8"):
            try:
                r = json.loads(ln)
            except Exception:
                continue
            r["_file"] = os.path.basename(f)
            rows.append(r)
    return rows


def main() -> int:
    L = []
    P = L.append
    P("diag_fabric_null.py -- sizing the fabric campaign. No DEM cells.")
    P("=" * 78)

    # --- 1. the closed form, checked against Monte Carlo -------------------------
    P("")
    P("1. THE FINITE-SAMPLE FLOOR.  Closed form E[A^2] = 1/n (derivation in the")
    P(f"   docstring), checked by Monte Carlo, {MC_REPS} reps, seed {MC_SEED}+n.")
    P("")
    P("     n     MC mean A    MC sd A    0.95/sqrt(n)    MC mean T = n*E[A^2]   MC sd T")
    worst = 0.0
    for n in NULL_NS:
        ma, sa, mt, st = mc_null(n)
        worst = max(worst, abs(mt - 1.0))
        P(f"  {n:6d}    {ma:.6f}   {sa:.6f}      {0.95/np.sqrt(n):.6f}"
          f"            {mt:.4f}          {st:.4f}")
    P("")
    P(f"   worst |MC E[T] - 1| over the table: {worst:.4f}")
    if worst > 0.05:
        P("   FAIL: the Monte Carlo does not reproduce the closed form.")
        return 1
    P("   PASS: the Monte Carlo reproduces E[T] = 1 to better than 5 %.")
    sd_t_null = float(np.mean([mc_null(n)[3] for n in (2000, 3000, 4000)]))
    P(f"   sd(T) under the null, averaged over n = 2000/3000/4000: {sd_t_null:.4f}")
    P("   (T is a quadratic form in the 5 independent deviator components, so it is")
    P("    chi2_5/5-like: predicted sd = sqrt(2/5) = 0.6325.)")

    # --- 2. the confound, measured -----------------------------------------------
    rows = load_shards()
    P("")
    P("2. THE CONTACT-COUNT CONFOUND, from committed shards.")
    P(f"   N before any filter: {len(rows)} records in out/phase2_*.jsonl")
    ok = [r for r in rows if not r.get("failed")]
    P(f"   N after dropping failed cells: {len(ok)}")
    sel = [r for r in ok if r.get("n_requested") == N_BODIES and r.get("mu_p") == MU]
    P(f"   N after N={N_BODIES}, mu_p={MU}: {len(sel)}")
    for p in (P_PRIMARY, P_SECOND):
        sub = [r for r in sel if r.get("p_hat_target") == p]
        P(f"   N after p_hat={p:g}: {len(sub)}")
        P(f"   mean n_contacts and the floor it implies, p_hat = {p:g}:")
        P("     family        " + "  ".join(f"CV{cv:.1f}" for cv in CVS)
          + "     floor CV0 -> CV0.5    slope/CV")
        for fam in ("lognormal", "uniform", "bidisperse"):
            vals = []
            for cv in CVS:
                v = [r["n_contacts"] for r in sub
                     if (r.get("family") or "lognormal") == fam and r["cv"] == cv]
                vals.append(np.mean(v) if v else np.nan)
            f0 = 0.95 / np.sqrt(vals[0]) if vals[0] == vals[0] else np.nan
            f5 = 0.95 / np.sqrt(vals[-1]) if vals[-1] == vals[-1] else np.nan
            P(f"     {fam:12s} " + " ".join(f"{v:6.0f}" for v in vals)
              + f"     {f0:.4f} -> {f5:.4f}     {(f5-f0)/0.5:+.4f}")
    P("")
    P("   READ THIS BEFORE READING ANY `aniso` vs CV PLOT: the floor alone rises by")
    P("   roughly a third across the CV ladder, in the SAME direction a restructuring")
    P("   signal would. T = n*A^2 removes it by construction; `aniso` does not.")

    # --- 3. the design ------------------------------------------------------------
    P("")
    P("3. RESOLVING POWER of the slope of T against CV.")
    P(f"   6 CV rungs {CVS}, sd(CV) = {np.std(CVS, ddof=0):.4f}, sd(T) = {sd_t_null:.4f}")
    P(f"   pre-registered effect size: DELTA_T over CV 0 -> 0.5 of {DELTA_T_TARGET}")
    P("")
    P("     seeds   points   se(slope)   2*se(slope)   2*se as DELTA_T over the range"
      "   resolves?")
    sdcv = float(np.std(CVS, ddof=0))
    chosen = None
    for s in (3, 6, 9, 12, 15):
        npts = len(CVS) * s
        se = sd_t_null / (np.sqrt(npts) * sdcv)
        dt = 2 * se * 0.5
        ok_ = dt <= DELTA_T_TARGET
        if ok_ and chosen is None:
            chosen = s
        P(f"     {s:5d}   {npts:6d}   {se:9.4f}   {2*se:11.4f}   {dt:29.4f}"
          f"   {'YES' if ok_ else 'no'}")
    P("")
    P(f"   SMALLEST SEED COUNT THAT RESOLVES THE REGISTERED EFFECT: {chosen}")
    P("   This is the number the prereg registers, and it assumes the measured sd(T)")
    P("   equals the null's. It will NOT, in general: if the packs carry real")
    P("   orientational structure the scatter is larger. The grader therefore")
    P("   RECOMPUTES the achieved resolving power from the measured sd and reports the")
    P("   verdict against that, so a null can come back UNINFORMATIVE rather than pass.")

    txt = "\n".join(L) + "\n"
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "diag_fabric_null.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

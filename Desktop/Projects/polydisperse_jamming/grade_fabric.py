"""grade_fabric.py -- grade the fabric campaign against PREREG_FABRIC_2026-09-18.md.

Reads out/fabric/phase2_fab*.jsonl ONLY. Writes out/grade_fabric.txt. Runs no cell.

Order is the prereg's order and it is load-bearing: F0 (reproduction) runs FIRST and
exits non-zero on a mismatch, because a diagnostic that cannot reproduce the numbers it
audits has not audited them. Every gate prints its registered SENTENCE beside its number,
so the sentence and the line of code that grades it can be read together.

THE OBSERVABLE, and why it is not `aniso`:

    T = n_contacts * aniso^2 ,   E[T] = 1 under the independent-normal null.

`n_contacts` falls by 1.8x across the CV ladder, so the finite-sample floor on `aniso`
rises by a third with CV all on its own, in the same direction a restructuring signal
would. T divides that out by construction. See diag_fabric_null.py for the closed form
E[aniso^2] = 1/n and its Monte-Carlo check through the shipped `fabric_tensor`.
"""

from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FAB = os.path.join(OUT, "fabric")
PREREG = os.path.join(HERE, "PREREG_FABRIC_2026-09-18.md")

CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
FAMILIES = ("lognormal", "uniform", "bidisperse")
P_PRIMARY = 1e-3
P_SECOND = 1e-4
DELTA_T_TARGET = 0.5           # registered effect size, as a change across CV 0 -> 0.5
K_PRIMARY = 2
K_REPORTED = (2, 3, 4)

# The floor model: E[aniso] = C_FLOOR / sqrt(n) for n independent unit normals. The
# 1/sqrt(n) is exact (E[aniso^2] = 1/n); C_FLOOR is the E[chi]/E[chi^2]^(1/2) correction
# and is MEASURED here by Monte Carlo rather than assumed, with a fixed seed.
FLOOR_MC_SEED = 20260918
FLOOR_MC_REPS = 2000

# F1: the three cells re-run under the fabric code, and the committed shard each one's
# reference `phi` is READ FROM (never a literal typed twice -- the literal in the prereg
# is cross-checked against the shard below).
F1_CELLS = {
    "cv0.300_mu0.300_s1_p1.000e-03_N1500_bidisperse": "phase2_bidisperse_shard00of12.jsonl",
    "cv0.300_mu0.300_s1_p1.000e-03_N1500_uniform": "phase2_uniform_shard00of12.jsonl",
    "cv0.400_mu0.300_s1_p1.000e-03_N1500": "phase2_shard00of12.jsonl",
}

BANNED_TAGS = ("mixr", "extr", "tern", "kl", "ug")

L: list[str] = []


def P(s=""):
    L.append(s)


def ols(x, y):
    """Slope, intercept and the slope's standard error, from scratch."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 3:
        return float("nan"), float("nan"), float("nan")
    sxx = float(np.sum((x - x.mean()) ** 2))
    slope = float(np.sum((x - x.mean()) * (y - y.mean())) / sxx)
    icpt = float(y.mean() - slope * x.mean())
    resid = y - (icpt + slope * x)
    s2 = float(np.sum(resid ** 2) / (n - 2))
    return slope, icpt, float(math.sqrt(s2 / sxx))


def c_floor():
    rng = np.random.default_rng(FLOOR_MC_SEED)
    out = []
    for n in (1000, 2000, 3000, 4000):
        v = rng.normal(size=(FLOOR_MC_REPS, n, 3))
        v /= np.linalg.norm(v, axis=2, keepdims=True)
        F = np.einsum("rni,rnj->rij", v, v) / n
        dev = F - np.eye(3) / 3.0
        a = np.sqrt(1.5 * np.einsum("rij,rij->r", dev, dev))
        out.append(float(a.mean()) * math.sqrt(n))
    return float(np.mean(out)), float(np.std(out, ddof=1))


# ------------------------------------------------------------------ loading + guards
def load() -> list[dict]:
    files = sorted(glob.glob(os.path.join(FAB, "phase2_fab*.jsonl")))
    stray = [f for f in sorted(glob.glob(os.path.join(FAB, "*.jsonl")))
             if os.path.basename(f) not in [os.path.basename(x) for x in files]]
    P("TAG GUARD (structural, not by eye)")
    P(f"  files matching out/fabric/phase2_fab*.jsonl : {len(files)}")
    for f in files:
        P(f"    {os.path.basename(f)}")
    if stray:
        P(f"  FAIL: {len(stray)} file(s) in out/fabric/ do NOT match the tag: {stray}")
        sys.exit(1)
    for f in files:
        base = os.path.basename(f)
        for bad in BANNED_TAGS:
            if re.search(rf"phase2_{bad}", base):
                P(f"  FAIL: banned tag {bad!r} in {base}")
                sys.exit(1)
    rows = []
    for f in files:
        for ln in open(f, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            r = json.loads(ln)
            r["_file"] = os.path.basename(f)
            rows.append(r)
    P(f"  N records read: {len(rows)}")
    nofab = [r for r in rows if not r.get("fabric") or not r.get("fabric_nr")]
    if nofab:
        P(f"  FAIL: {len(nofab)} record(s) lack a non-empty fabric/fabric_nr")
        sys.exit(1)
    P("  every record carries a non-empty `fabric` AND `fabric_nr`: PASS")
    return rows


# ----------------------------------------------------------------------------- F0
def f0(rows) -> bool:
    P("")
    P("F0  REPRODUCTION -- runs FIRST, exits non-zero on mismatch.")
    P("    SENTENCE: every record's fabric['n'] equals its own n_contacts, and at each")
    P("    k its fabric_nr['n_k{k}'] equals z_nr_k{k}*(1-rattler_frac_k{k})*n_placed/2,")
    P("    to within 0.5 of an integer; where the backbone is EMPTY both sides must be")
    P("    empty together.")
    bad = 0
    worst = 0.0
    empty_k = {k: 0 for k in K_REPORTED}
    for r in rows:
        if r["fabric"]["n"] != r["n_contacts"]:
            bad += 1
            continue
        rat = r["rattler"]
        for k in K_REPORTED:
            lhs = r["fabric_nr"][f"n_k{k}"]
            z = rat[f"z_nr_k{k}"]
            fr = rat[f"rattler_frac_k{k}"]
            if lhs == 0 or z != z:            # NaN test, no numpy needed
                if not (lhs == 0 and z != z):
                    bad += 1
                else:
                    empty_k[k] += 1
                continue
            rhs = z * (1.0 - fr) * r["n_placed"] / 2.0
            worst = max(worst, abs(lhs - rhs))
            if abs(lhs - rhs) >= 0.5:
                bad += 1
    P(f"    records checked: {len(rows)};  worst |n_k - identity|: {worst:.4f}")
    P("    empty-backbone cells by k: " + ", ".join(f"k{k}={v}" for k, v in empty_k.items())
      + "   (k=4 is registered DEGENERATE at mu_p = 0.3)")
    P(f"    VERDICT: {'PASS' if bad == 0 else f'FAIL on {bad} check(s)'}")
    return bad == 0


# ----------------------------------------------------------------------------- F1
def f1(rows) -> bool:
    P("")
    P("F1  BIT IDENTITY -- the fabric read must not move the solver.")
    P("    SENTENCE: three cells that already exist on disk return `phi` BIT-IDENTICAL")
    P("    under the fabric code. `phi` is the gate because it is the CONTINUOUS")
    P("    channel; z_all, z_nr and rattler_frac are grid-quantized and are printed.")
    prereg = open(PREREG, encoding="utf-8").read()
    ok = True
    seen = 0
    for key, srcfile in F1_CELLS.items():
        ref = None
        for ln in open(os.path.join(OUT, srcfile), encoding="utf-8"):
            r = json.loads(ln)
            if r.get("key") == key:
                ref = r
                break
        if ref is None:
            P(f"    FAIL: reference cell {key} not found in {srcfile}")
            ok = False
            continue
        # The prereg quotes the same number. Pin the document to the artifact, so a
        # hardcoded value that goes stale is caught here and not by a reader.
        if repr(ref["phi"]) not in prereg:
            P(f"    FAIL: PREREG does not quote {srcfile}'s phi = {ref['phi']!r}")
            ok = False
        new = [r for r in rows if r.get("key") == key and "fabgate" in r["_file"]]
        if not new:
            P(f"    FAIL: no fabgate re-run of {key}")
            ok = False
            continue
        n = new[0]
        seen += 1
        same = n["phi"] == ref["phi"]
        ok = ok and same
        P(f"    {key}")
        P(f"      phi  new {n['phi']!r}   ref {ref['phi']!r}   BIT-IDENTICAL {same}")
        P(f"      quantized channels, printed not read: z_all {n['z_all'] == ref['z_all']}"
          f"  z_nr {n['z_nr'] == ref['z_nr']}  rattler_frac"
          f" {n['rattler_frac'] == ref['rattler_frac']}")
    P(f"    cells checked: {seen}/{len(F1_CELLS)}")
    P(f"    VERDICT: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------------------------------------------- production set
def production(rows):
    P("")
    P("PRODUCTION SET -- N printed before and after every filter.")
    P(f"  N all records: {len(rows)}")
    r1 = [r for r in rows if "fabgate" not in r["_file"]]
    P(f"  N after dropping the F1 gate cells: {len(r1)}")
    r2 = [r for r in r1 if not r.get("failed")]
    P(f"  N after dropping failed cells: {len(r2)}")
    r3 = [r for r in r2 if r.get("converged")]
    P(f"  N after requiring converged: {len(r3)}   (dropped {len(r2)-len(r3)})")
    P("")
    P("  the d_max/L > 0.25 finite-size flag, BY FAMILY (never pooled):")
    for fam in FAMILIES:
        sub = [r for r in r3 if (r.get("family") or "lognormal") == fam]
        flag = [r for r in sub if any("FINITE SIZE" in s for s in r.get("notes", []))]
        cvs = sorted({r["cv"] for r in flag})
        P(f"    {fam:12s} {len(flag):3d}/{len(sub):3d} flagged   at CV {cvs}")
    return r3


def tstat(rows, channel, k=K_PRIMARY):
    """(cv, T) pairs for one channel: 'all' or 'nr'."""
    out = []
    for r in rows:
        if channel == "all":
            a, n = r["fabric"]["aniso"], r["fabric"]["n"]
        else:
            a, n = r["fabric_nr"][f"aniso_k{k}"], r["fabric_nr"][f"n_k{k}"]
        if a != a or not n:
            continue
        out.append((r["cv"], n * a * a))
    return out


def slope_block(rows, channel, k, title):
    P("")
    P(f"  {title}")
    P("    family        slope dT/dCV      se     2*se    achieved DELTA_T over range"
      "   |slope|>2se?")
    res = {}
    for fam in FAMILIES:
        sub = [r for r in rows if (r.get("family") or "lognormal") == fam]
        pts = tstat(sub, channel, k)
        if len(pts) < 6:
            P(f"    {fam:12s} only {len(pts)} points -- not graded")
            continue
        x = [p[0] for p in pts]
        y = [p[1] for p in pts]
        s, _, se = ols(x, y)
        res[fam] = (s, se, len(pts), float(np.mean(y)), float(np.std(y, ddof=1)))
        P(f"    {fam:12s} {s:+12.4f}  {se:7.4f}  {2*se:7.4f}   {2*se*0.5:25.4f}"
          f"   {'YES' if abs(s) > 2*se else 'no':>12s}")
    return res


# ----------------------------------------------------------------------------- main
def main() -> int:
    if not os.path.isdir(FAB):
        print(f"no {FAB}; nothing to grade")
        return 1
    P("grade_fabric.py -- PREREG_FABRIC_2026-09-18.md")
    P("=" * 78)
    rows = load()

    ok0 = f0(rows)
    ok1 = f1(rows)
    if not (ok0 and ok1):
        write()
        return 1

    prod = production(rows)
    prim = [r for r in prod if r["p_hat_target"] == P_PRIMARY]
    sec = [r for r in prod if r["p_hat_target"] == P_SECOND]
    P("")
    P(f"  N at the primary rung p_hat = {P_PRIMARY:g}: {len(prim)}")
    P(f"  N at the secondary rung p_hat = {P_SECOND:g}: {len(sec)}")

    # --- F2 ---------------------------------------------------------------------
    P("")
    P("F2  PRIMARY -- the registered prediction is the NULL.")
    P("    SENTENCE: if polydispersity only depopulates the load-bearing network, the")
    P("    orientational correlation of the surviving backbone does not change with CV:")
    P("    the per-family OLS slope of T_nr against CV is consistent with zero at 2*se,")
    P("    on all three families.")
    P(f"    registered effect size: DELTA_T = {DELTA_T_TARGET} across CV 0 -> 0.5.")
    P("    VERDICT GRAMMAR, registered in advance:")
    P("      |slope| > 2*se            -> DRIFT (the null is refuted on that family)")
    P("      |slope| <= 2*se and achieved 2*se*0.5 <= 0.5 -> NULL, INFORMATIVE")
    P("      |slope| <= 2*se and achieved 2*se*0.5 >  0.5 -> UNINFORMATIVE, not a pass")
    verdicts = {}
    res_k = {}
    for k in (2, 3):
        res = slope_block(prim, "nr", k, f"T_nr at k = {k}"
                          + ("  (PRIMARY)" if k == K_PRIMARY else "  (reported)"))
        res_k[k] = res
        if k == K_PRIMARY:
            for fam, (s, se, n, m, sd) in res.items():
                if abs(s) > 2 * se:
                    verdicts[fam] = "DRIFT"
                elif 2 * se * 0.5 <= DELTA_T_TARGET:
                    verdicts[fam] = "NULL, INFORMATIVE"
                else:
                    verdicts[fam] = "UNINFORMATIVE"
    P("")
    P("    mean T_nr (k = 2) and its scatter -- T = 1 is the independent-normal null,")
    P("    and the ABSOLUTE value carries no claim (a contact network is not")
    P("    independent; see the prereg's honest limit).")
    for fam, (s, se, n, m, sd) in res_k[K_PRIMARY].items():
        P(f"      {fam:12s} n = {n:3d}   mean T = {m:.4f}   sd = {sd:.4f}"
          f"   (null sd 0.6375)")
    P("")
    for fam in FAMILIES:
        P(f"    F2 VERDICT {fam:12s}: {verdicts.get(fam, 'no data')}")
    drift = [f for f, v in verdicts.items() if v == "DRIFT"]

    # --- F3 ---------------------------------------------------------------------
    P("")
    P("F3  THE ALL-CONTACT CHANNEL, reported beside it. No prediction registered.")
    slope_block(prim, "all", K_PRIMARY, "T_all (every real contact)")
    P("")
    P("    paired per-cell contrast T_all - T_nr(k=2), by CV -- the size of the")
    P("    conflation the split exists to avoid:")
    for fam in FAMILIES:
        sub = [r for r in prim if (r.get("family") or "lognormal") == fam]
        line = []
        for cv in CVS:
            d = [r["fabric"]["n"] * r["fabric"]["aniso"] ** 2
                 - r["fabric_nr"]["n_k2"] * r["fabric_nr"]["aniso_k2"] ** 2
                 for r in sub if r["cv"] == cv and r["fabric_nr"]["n_k2"]]
            line.append(f"{np.mean(d):+7.4f}" if d else "    n/a")
        P(f"      {fam:12s} " + " ".join(line))

    # --- F4 ---------------------------------------------------------------------
    cf, cf_sd = c_floor()
    P("")
    P("F4  THE CONFOUND MUST BE DEMONSTRATED, NOT ASSERTED.")
    P("    SENTENCE: the raw slope of `aniso` against CV is statistically")
    P("    indistinguishable from the slope the finite-sample floor produces on its own,")
    P("    on all three families, at 2*se of the PAIRED difference.")
    P(f"    floor model E[aniso] = C/sqrt(n), C measured by Monte Carlo (seed"
      f" {FLOOR_MC_SEED}, {FLOOR_MC_REPS} reps): C = {cf:.4f} +- {cf_sd:.4f}")
    P("")
    P("    family        slope(aniso)   slope(floor)   slope(aniso-floor)      se"
      "    consistent?")
    f4_ok = True
    for fam in FAMILIES:
        sub = [r for r in prim if (r.get("family") or "lognormal") == fam]
        if len(sub) < 6:
            continue
        cv = [r["cv"] for r in sub]
        a = [r["fabric"]["aniso"] for r in sub]
        fl = [cf / math.sqrt(r["fabric"]["n"]) for r in sub]
        d = [ai - fi for ai, fi in zip(a, fl)]
        sa, _, _ = ols(cv, a)
        sf, _, _ = ols(cv, fl)
        sd_, _, sed = ols(cv, d)
        cons = abs(sd_) <= 2 * sed
        f4_ok = f4_ok and cons
        P(f"    {fam:12s} {sa:+12.5f}   {sf:+12.5f}   {sd_:+18.5f}  {sed:.5f}"
          f"    {'YES' if cons else 'NO'}")
    P(f"    F4 VERDICT: {'the raw channel is a readout of the contact count'  if f4_ok else 'A REAL DEVIATORIC SIGNAL SURVIVES THE FLOOR -- read this before anything else'}")

    # --- F5 ---------------------------------------------------------------------
    P("")
    P("F5  RUNG CONSISTENCY -- registered as a CONSISTENCY CHANNEL ONLY; at 6 seeds")
    P("    2*se*0.5 = 0.622 > 0.5, so it cannot decide F2.")
    if sec:
        pairs = []
        idx = {(r.get("family") or "lognormal", r["cv"], r["seed"]):
               r["fabric_nr"]["n_k2"] * r["fabric_nr"]["aniso_k2"] ** 2
               for r in sec if r["fabric_nr"]["n_k2"]}
        for r in prim:
            key = ((r.get("family") or "lognormal"), r["cv"], r["seed"])
            if key in idx and r["fabric_nr"]["n_k2"]:
                pairs.append(r["fabric_nr"]["n_k2"] * r["fabric_nr"]["aniso_k2"] ** 2
                             - idx[key])
        if pairs:
            m = float(np.mean(pairs))
            se = float(np.std(pairs, ddof=1) / math.sqrt(len(pairs)))
            P(f"    paired T_nr(1e-3) - T_nr(1e-4) on {len(pairs)} shared cells:"
              f" {m:+.4f} +- {se:.4f}")
        else:
            P("    no shared cells")
    else:
        P("    the secondary rung has not been run")

    # --- F7 ---------------------------------------------------------------------
    P("")
    P("F7  CONDITIONAL -- readable only if F2 finds a drift.")
    if not drift:
        P("    NOT READ: F2 found no drift, so there is nothing to attribute.")
    else:
        P(f"    F2 found a drift on {drift}. Does T_nr collapse on z_nr across families")
        P("    at matched CV? Within a family z and CV are collinear by construction.")
        P("      family        CV      mean z_nr    mean T_nr")
        for cv in CVS:
            for fam in FAMILIES:
                sub = [r for r in prim if (r.get("family") or "lognormal") == fam
                       and r["cv"] == cv and r["fabric_nr"]["n_k2"]]
                if not sub:
                    continue
                z = float(np.mean([r["z_nr"] for r in sub]))
                t = float(np.mean([r["fabric_nr"]["n_k2"]
                                   * r["fabric_nr"]["aniso_k2"] ** 2 for r in sub]))
                P(f"      {fam:12s} {cv:.1f}   {z:9.4f}    {t:9.4f}")

    P("")
    P("=" * 78)
    P("F0 " + ("PASS" if ok0 else "FAIL") + " | F1 " + ("PASS" if ok1 else "FAIL")
      + " | F2 " + " ".join(f"{f}={verdicts.get(f, 'n/a')}" for f in FAMILIES))
    write()
    return 0


def write():
    os.makedirs(OUT, exist_ok=True)
    txt = "\n".join(L) + "\n"
    with open(os.path.join(OUT, "grade_fabric.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    raise SystemExit(main())

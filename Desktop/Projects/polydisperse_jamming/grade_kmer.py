"""grade_kmer.py -- grades PREREG_KMER_2026-09-22.md. Writes out/grade_kmer.txt.

Written before the mixture cells landed. Reads out/shape/phase2_kmer_*.jsonl, the shapechi
endpoint records, and out/kmer_dump/*.npz.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

import actshape as A

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PHI_SCALE = 0.0021
MUS = (0.0, 0.3, 3.0)
RUNGS = (1e-3, 3e-5)
COMPS = ((0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0), (0.0, 0.5), (0.0, 1.0))
T_CRIT = {1: 12.706, 2: 4.303}


def load(pattern):
    recs = {}
    for f in sorted(glob.glob(os.path.join(OUT, "shape", pattern))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            recs[r["key"]] = r
    return recs


def hmix(p1, p2, w):
    return 1.0 / ((1 - w) / p1 + w / p2)


def main() -> int:
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    km = load("phase2_kmer_*.jsonl")
    sc = {k: r for k, r in load("phase2_shapechi*.jsonl").items() if r.get("cv") == 0.0}
    say("grade_kmer -- PREREG_KMER_2026-09-22.md")
    say(f"records: {len(km)} of 126")
    fails = [r for r in km.values() if r.get("failed")]
    for r in fails:
        say(f"  FAILED CELL {r['key']}: {r.get('error')}")
    good = {k: r for k, r in km.items() if not r.get("failed")}

    # ---- K3 health --------------------------------------------------------------
    bad = [k for k, r in good.items()
           if not (r["converged"] and abs(r["mass_over_vol_solid"] - 1) < 1e-9
                   and r["n_contacts_internal"] == 0)]
    say(f"K3 HEALTH: {len(good) - len(bad)}/{len(good)} converged, mass/vol = 1, no internal")
    for k in bad:
        say(f"  EXCLUDED {k}: conv {good[k]['converged']} mass {good[k]['mass_over_vol_solid']:.9f}"
            f" internal {good[k]['n_contacts_internal']} notes {good[k]['notes']}")
    ok = {k: r for k, r in good.items() if k not in bad}

    # ---- K4 (POST-HOC, PREREG_KMER Amendment 2): overlap sanity ---------------------
    # Linear FrictPhys between equal primaries: kn = 2E r1 r2/(r1+r2) = E r, so a contact's
    # overlap is u = fn/(E r) and u/r = fn/(E r^2) = fn/0.25 (E = 1, r = 1/2). Frictional
    # triplet packs were found with 7-13 % of contacts above u/r = 0.05 and some above 1 --
    # interpenetrated primaries, not a rigid-body packing, and an inflated phi. Stricter, not
    # looser: a cell whose dump shows > 0.1 % of contacts above u/r = 0.05 is EXCLUDED.
    k4_bad = {}
    for k in list(ok):
        p = os.path.join(OUT, "kmer_dump", k + ".npz")
        if not os.path.exists(p):
            continue
        fn = np.load(p)["fn"]
        frac = float((fn / 0.25 > 0.05).mean()) if fn.size else 0.0
        if frac > 1e-3:
            k4_bad[k] = frac
    say(f"K4 OVERLAP SANITY (post-hoc, Amendment 2): {len(k4_bad)} cells excluded "
        f"(> 0.1 % of contacts with overlap > 5 % of the radius)")
    for k, frac in sorted(k4_bad.items()):
        say(f"  EXCLUDED {k[:64]}: {frac:.2%} of contacts deep-overlapped")
    ok = {k: r for k, r in ok.items() if k not in k4_bad}
    fs = [k for k, r in ok.items() if r["d_max_over_L"] > 0.25]
    say(f"finite-size flag d_max/L > 0.25 on {len(fs)} cells "
        f"(max {max((r['d_max_over_L'] for r in ok.values()), default=float('nan')):.3f})")

    # ---- K0' singlet path identity; K2 scaled copy --------------------------------
    k0, k2 = [], []
    for k, r in ok.items():
        w2, w3 = r["doublet_w"], r["triplet_w"]
        if w2 == 0 and w3 == 0:
            ref = sc.get(k)
            if ref is not None:
                k0.append(r["phi"] == ref["phi"] and r["iters"] == ref["iters"])
        if w2 == 1.0 and w3 == 0:
            ref = next((x for x in sc.values() if x["alpha"] == 2.0 and x["mu_p"] == r["mu_p"]
                        and x["seed"] == r["seed"] and x["p_hat_target"] == r["p_hat_target"]),
                       None)
            if ref is not None:
                k2.append(abs(r["phi"] - ref["phi"]))
    say(f"K0' singlet path == shapechi alpha 1 (phi and iters exact): {sum(k0)}/{len(k0)}"
        + ("  PASS" if k0 and all(k0) else "  FAIL/INCOMPLETE"))
    say(f"K2 w2 = 1 scaled copy of shapechi alpha 2: worst |dphi| = "
        f"{max(k2) if k2 else float('nan'):.2e} over {len(k2)} cells (bar 1e-4)"
        + ("  PASS" if k2 and max(k2) < 1e-4 else "  FAIL/INCOMPLETE"))

    # ---- group ----------------------------------------------------------------------
    grp = defaultdict(list)
    for r in ok.values():
        grp[(r["doublet_w"], r["triplet_w"], r["mu_p"], r["p_hat_target"])].append(r)

    def m(comp, mu, rung, field="phi"):
        v = [x[field] for x in grp.get((comp[0], comp[1], mu, rung), [])]
        return (float(np.mean(v)), len(v)) if v else (float("nan"), 0)

    # endpoints from THIS campaign (w2 = 0 and 1) -- K0'/K2 say they equal shapechi's
    say("")
    say("=== H_mix: phi(w2) vs ideal specific-volume mixing of the campaign's own endpoints ===")
    say("   bar |phi - H_mix| < 2 PHI_SCALE = 0.0042 at every (mu, rung, w2 in .25/.5/.75)")
    hm_ok = []
    for mu in MUS:
        for rung in RUNGS:
            p0, n0 = m((0.0, 0.0), mu, rung)
            p1, n1 = m((1.0, 0.0), mu, rung)
            row = []
            for w in (0.25, 0.5, 0.75):
                pw, nw = m((w, 0.0), mu, rung)
                pred = hmix(p0, p1, w)
                dev = pw - pred
                if nw and n0 and n1:
                    hm_ok.append(abs(dev) < 2 * PHI_SCALE)
                row.append(f"w {w:.2f}: {pw:.4f} vs {pred:.4f} ({dev:+.4f}, n {nw})")
            say(f"  mu {mu:.1f} p {rung:.0e}: ends {p0:.4f}/{p1:.4f} | " + " | ".join(row))
    say(f"H_mix: {sum(hm_ok)}/{len(hm_ok)} within bar"
        + ("  PASS" if hm_ok and all(hm_ok) and len(hm_ok) == 18 else
           "  FAIL" if hm_ok and not all(hm_ok) else "  INCOMPLETE"))
    excess = []
    for mu in MUS:
        for rung in RUNGS:
            p0, _ = m((0.0, 0.0), mu, rung)
            p1, _ = m((1.0, 0.0), mu, rung)
            pw, nw = m((0.5, 0.0), mu, rung)
            if nw:
                excess.append((mu, rung, pw - hmix(p0, p1, 0.5)))
    say("   excess at w2 = 0.5: " + ", ".join(f"mu {a:.1f} p {b:.0e}: {c:+.4f}" for a, b, c in excess))

    # ---- H_chi ------------------------------------------------------------------------
    say("")
    say("=== H_chi: chi(w2) = phi(mu 3)/phi(mu 0) and phi(mu .3)/phi(mu 0), seed-paired ===")
    reg = {0.0: 0.886, 0.25: 0.859, 0.5: 0.834, 0.75: 0.811, 1.0: 0.788}
    chi_ok = []
    for rung in RUNGS:
        for w in (0.0, 0.25, 0.5, 0.75, 1.0):
            vals = {}
            for mm in (3.0, 0.3):
                a = {x["seed"]: x["phi"] for x in grp.get((w, 0.0, mm, rung), [])}
                b = {x["seed"]: x["phi"] for x in grp.get((w, 0.0, 0.0, rung), [])}
                s = sorted(set(a) & set(b))
                vals[mm] = (float(np.mean([a[x] / b[x] for x in s])) if s else float("nan"),
                            len(s))
            line = (f"  p {rung:.0e} w2 {w:.2f}: chi(3) {vals[3.0][0]:.4f} (n {vals[3.0][1]})"
                    f"  chi(.3) {vals[0.3][0]:.4f} (n {vals[0.3][1]})")
            if rung == 3e-5 and vals[3.0][1] == 0:
                line += f"   registered {reg[w]:.3f} +- 0.01 (no data)"
            elif rung == 3e-5:
                good_c = abs(vals[3.0][0] - reg[w]) <= 0.01
                chi_ok.append(good_c)
                line += f"   registered {reg[w]:.3f} +- 0.01 {'ok' if good_c else 'MISS'}"
            say(line)
    say(f"H_chi (p 3e-5, mu_m 3): {sum(chi_ok)}/{len(chi_ok)} within +-0.01"
        + ("  PASS" if len(chi_ok) == 5 and all(chi_ok) else
           "  FAIL" if chi_ok and not all(chi_ok) else "  INCOMPLETE"))

    # ---- H_Z --------------------------------------------------------------------------
    say("")
    say("=== H_Z: isostatic counting ===")
    z_ok = []
    for comp in COMPS:
        for rung in RUNGS:
            z3, n3 = m(comp, 3.0, rung, "z_nr")
            if n3:
                z_ok.append(abs(z3 - 4.0) <= 0.25 if rung == 3e-5 else True)
            say(f"  comp w2 {comp[0]:.2f} w3 {comp[1]:.2f} p {rung:.0e}: Z_nr(mu 3) {z3:.3f}"
                f"  Z_nr(mu 0) {m(comp, 0.0, rung, 'z_nr')[0]:.3f}")
    say(f"H_Z rough (p 3e-5): {sum(z_ok)}/{len(z_ok)} within 4 +- 0.25"
        + ("  PASS" if z_ok and all(z_ok) else "  FAIL"))
    zl = []
    z0s, _ = m((0.0, 0.0), 0.0, 3e-5, "z_nr")
    z0d, _ = m((1.0, 0.0), 0.0, 3e-5, "z_nr")
    for w in (0.25, 0.5, 0.75):
        n2 = (w / 2) / ((1 - w) + w / 2)
        pred = z0s + n2 * (z0d - z0s)
        zw, nw = m((w, 0.0), 0.0, 3e-5, "z_nr")
        if nw:
            zl.append(abs(zw - pred) <= 0.30)
        say(f"  frictionless w2 {w:.2f} (n2 {n2:.3f}): Z {zw:.3f} vs linear-in-n2 {pred:.3f}")
    say(f"H_Z frictionless linear-in-n2: {sum(zl)}/{len(zl)} within 0.30"
        + ("  PASS" if zl and all(zl) else "  FAIL"))

    # ---- H_sigma*, H_delta from dumps --------------------------------------------------
    say("")
    say("=== H_sigma* and H_delta from per-contact dumps (R* = 1/4 for every contact) ===")
    dm = defaultdict(dict)
    for k, r in ok.items():
        p = os.path.join(OUT, "kmer_dump", k + ".npz")
        if not os.path.exists(p):
            continue
        d = np.load(p)
        fn = d["fn"][d["fn"] > 0]
        g = fn / (0.25 * float(d["p_hat"]))
        dm[(r["doublet_w"], r["triplet_w"], r["mu_p"], r["p_hat_target"])][r["seed"]] = {
            "sstar": A.sigma_star_over_kappa(g), "d": A.delta_from_quantiles(g, 0.1, 0.9),
            "d2": A.delta_from_quantiles(g, 0.25, 0.75), "fcv": float(fn.std() / fn.mean()),
            "q99": float(np.quantile(fn, 0.99) / fn.mean()), "fmax": float(fn.max() / fn.mean())}
    for mu in MUS:
        for rung in RUNGS:
            for comp in COMPS:
                c = dm.get((comp[0], comp[1], mu, rung))
                if not c:
                    continue
                vs = list(c.values())
                say(f"  mu {mu:.1f} p {rung:.0e} w2 {comp[0]:.2f} w3 {comp[1]:.2f}: sigma*/kappa "
                    f"{np.mean([v['sstar'] for v in vs]):.4g}  delta10/90 "
                    f"{np.mean([v['d'] for v in vs]):.3f}  delta25/75 "
                    f"{np.mean([v['d2'] for v in vs]):.3f}  force cv {np.mean([v['fcv'] for v in vs]):.3f}"
                    f"  f99/<f> {np.mean([v['q99'] for v in vs]):.2f}  fmax/<f> "
                    f"{np.mean([v['fmax'] for v in vs]):.1f}  (n {len(vs)})")
    for rung in RUNGS:
        a = dm.get((0.0, 0.0, 0.3, rung), {})
        b = dm.get((1.0, 0.0, 0.3, rung), {})
        s = sorted(set(a) & set(b))
        if len(s) >= 2:
            rat = [b[x]["sstar"] / a[x]["sstar"] for x in s]
            dd = [b[x]["d"] - a[x]["d"] for x in s]
            tcrit = T_CRIT[len(s) - 1]
            se = np.std(dd, ddof=1) / math.sqrt(len(dd))
            t = np.mean(dd) / se if se > 0 else float("inf")
            tag = " [GRADED]" if rung == 3e-5 else " [reported]"
            say(f"  p {rung:.0e} H_sigma*: ratio doublet/singlet {np.mean(rat):.3f} "
                f"(registered [0.55, 0.95]) "
                f"{'PASS' if 0.55 <= np.mean(rat) <= 0.95 else 'FAIL'}{tag}")
            say(f"  p {rung:.0e} H_delta: delta(doublet) - delta(singlet) {np.mean(dd):+.3f}, "
                f"t {t:.2f} (t_crit {tcrit}) "
                f"{'PASS' if np.mean(dd) < 0 and abs(t) >= tcrit else 'FAIL'}{tag}")
    with open(os.path.join(OUT, "grade_kmer.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

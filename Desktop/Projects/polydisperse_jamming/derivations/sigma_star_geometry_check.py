"""sigma_star_geometry_check.py -- checks for derivations/sigma_star_geometry_2026-09-24.md.

    python derivations/sigma_star_geometry_check.py   ->  out/sigma_star_geometry.txt

1. SymPy: R*(a_i+a_j) = a_i a_j; the pair-model closed form for <a_i a_j>_c; the
   monodisperse limit of the decomposition.
2. Monte Carlo: the pair-model closed form against brute-force weighted pair sampling.
3. Data: the exact Love-Weber decomposition of sigma*/kappa on every actshape dump
   (sheared constant-pressure Yade packs, per-contact fn, r1, r2):
      sigma*/kappa = (N_c/V) <aa>_c C / (3 q LW)                               (eq. 7)
   with C = <g aa>/(<g><aa>) and q = Q_{1-1/e}(g)/<g>. The identity is checked
   (LW = sum fn (r1+r2) / (3 V p) must be 1 up to overlap), then the CV 0 -> 0.5
   ratio of sigma*/kappa is split into the four factors, seed-paired.
"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import actshape as A  # noqa: E402

INV_E = math.exp(-1.0)
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


# ---------------------------------------------------------------- 1. SymPy
def sympy_checks() -> None:
    say("=== 1. SymPy ===")
    ai, aj = sp.symbols("a_i a_j", positive=True)
    Rs = ai * aj / (ai + aj)
    say(f"  R*(a_i+a_j) - a_i a_j = {sp.simplify(Rs*(ai+aj) - ai*aj)}   (eq. 3)")
    # pair model A: P(i,j) ~ n_i n_j (a_i+a_j)^2, a_i, a_j iid from the number PSD
    m1, m2, m3 = sp.symbols("m1 m2 m3", positive=True)
    num = sp.expand(ai * aj * (ai + aj) ** 2)          # a_i^3 a_j + 2 a_i^2 a_j^2 + a_i a_j^3
    den = sp.expand((ai + aj) ** 2)
    # replace monomials a_i^p a_j^q by m_p m_q (independence)
    def avg(expr):
        expr = sp.Poly(expr, ai, aj)
        tot = 0
        mm = {0: 1, 1: m1, 2: m2, 3: m3}
        for (p, q), c in expr.terms():
            tot += c * mm[p] * mm[q]
        return sp.simplify(tot)
    aaA = sp.simplify(avg(num) / avg(den))
    target = (m3 * m1 + m2 ** 2) / (m2 + m1 ** 2)
    say(f"  <aa>_c (model A) - (m3 m1 + m2^2)/(m2 + m1^2) = {sp.simplify(aaA - target)}   (eq. 9)")
    # closure C: P(i,j) ~ a_i a_j (a_i + a_j)  (referee pass 1, m2: check it too)
    mm4 = sp.symbols("m4", positive=True)
    numC = sp.expand(ai * aj * ai * aj * (ai + aj))
    denC = sp.expand(ai * aj * (ai + aj))

    def avg4(expr):
        expr = sp.Poly(expr, ai, aj)
        mm = {0: 1, 1: m1, 2: m2, 3: m3, 4: mm4}
        return sp.simplify(sum(c * mm[p] * mm[q] for (p, q), c in expr.terms()))
    say(f"  <aa>_c (closure C) - m3/m1 = {sp.simplify(avg4(numC) / avg4(denC) - m3 / m1)}   (eq. 12, any PSD)")
    # monodisperse: m_n = a^n
    a = sp.symbols("a", positive=True)
    say(f"  monodisperse <aa>_c = {sp.simplify(target.subs({m1: a, m2: a**2, m3: a**3}))}  (expect a^2)")
    # lognormal at fixed number mean 1: m_n = s^(n(n-1)/2), s = 1 + CV^2
    s = sp.symbols("s", positive=True)
    ln = {m1: 1, m2: s, m3: s ** 3}
    geo = sp.simplify(target.subs(ln) / s ** 3)
    say(f"  lognormal <aa>_c/<a^3> (model A) = {sp.factor(geo)}   (eq. 10)")
    return None


# ---------------------------------------------------------------- 2. Monte Carlo
def mc_check() -> None:
    say("\n=== 2. Monte Carlo, eq. (9) vs weighted pair sampling (lognormal, <a> = 1) ===")
    rng = np.random.default_rng(20260924)
    for cv in (0.2, 0.35, 0.5):
        sl = math.sqrt(math.log1p(cv * cv))
        a = rng.lognormal(-0.5 * sl * sl, sl, (2, 2_000_000))
        w = (a[0] + a[1]) ** 2
        mc = float(np.sum(w * a[0] * a[1]) / np.sum(w))
        s = 1 + cv * cv
        cf = (s ** 3 + s ** 2) / (s + 1)
        say(f"  CV {cv:.2f}: MC {mc:.5f}  closed {cf:.5f}  rel err {mc/cf-1:+.2e}")
        wC = a[0] * a[1] * (a[0] + a[1])
        mcC = float(np.sum(wC * a[0] * a[1]) / np.sum(wC))
        say(f"           closure C: MC {mcC:.5f}  m3/m1 = s^3 = {s**3:.5f}  rel err {mcC/s**3-1:+.2e}")


# ---------------------------------------------------------------- 3. data
def load_recs() -> dict:
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase9_actshape_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            recs[r["key"]] = r
    return recs


def cell_factors(r: dict, d) -> dict | None:
    """Per cell, pooled over the 10 steady snapshots. Lengths in the dump's units."""
    N = r["n_placed"]
    d10, cv = r["d10_realised"], r["cv_realised"]
    m2d = d10 ** 2 * (1 + cv ** 2)           # <d^2>, number CV of the realised draw
    m3d = r["d32_realised"] * m2d            # <d^3>
    V = N * math.pi / 6.0 * m3d / r["phi_ss"]
    E = float(d["young"])
    edges = np.concatenate([[0], np.cumsum(d["n"])])
    lw, Nc, aa, C, q, ss, cfg = [], [], [], [], [], [], []
    for k in range(len(d["n"])):
        sl = slice(edges[k], edges[k + 1])
        fn, r1, r2 = d["fn"][sl], d["r1"][sl], d["r2"][sl]
        ok = fn > 0
        p = float(d["p_hat"][k]) * E
        lw.append(float(np.sum(fn[ok] * (r1[ok] + r2[ok])) / (3 * V * p)))
        g = A.g_values(fn, r1, r2, p, "derjaguin")
        prod = (r1 * r2)[ok]
        ends = np.concatenate([r1[ok], r2[ok]])        # contact-end radii: mean = sum k a / sum k
        cfg.append(float(prod.mean() / ends.mean() ** 2))
        Nc.append(int(ok.sum()))
        aa.append(float(prod.mean()))
        C.append(float(np.mean(g * prod) / (g.mean() * prod.mean())))
        q.append(float(np.quantile(g, 1 - INV_E) / g.mean()))
        ss.append(A.sigma_star_over_kappa(g))
    out = {"LW": float(np.mean(lw)), "LWmin": float(np.min(lw)), "LWmax": float(np.max(lw)),
           "ncV": float(np.mean(Nc)) / V, "aa": float(np.mean(aa)), "C": float(np.mean(C)),
           "q": float(np.mean(q)), "sstar": float(np.mean(ss)), "phi": r["phi_ss"],
           "z": r.get("z_ss"), "cv": r["cv"], "cvr": cv, "a3": m3d / 8.0,
           "N": N, "Z": 2.0 * float(np.mean(Nc)) / N, "cfg": float(np.mean(cfg))}
    # realised radius moments (radius = d/2): m1..m4
    m1, m2, m3 = d10 / 2, m2d / 4, m3d / 8
    out["aaC"] = m3 / m1                                 # eq. (12): P(i,j) ~ a_i a_j (a_i+a_j)
    out["aaA"] = (m3 * m1 + m2 ** 2) / (m2 + m1 ** 2)   # eq. (9): P(i,j) ~ (a_i+a_j)^2
    out["aaB"] = (m3 / m2) ** 2                          # eq. (11): P(i,j) ~ a_i^2 a_j^2
    # BMCSL (Santos, Yuste & Lopez de Haro 2002, doi:10.1063/1.1502247, Eq. 7 with lambda = 1/2):
    # w = (a+b)^2 g_ab, g_ab = 1/(1-phi) + 3 xi2/(1-phi)^2 (2R*) + 2 xi2^2/(1-phi)^3 (2R*)^2, xi2 = phi m2/(2 m3)
    ph = r["phi_ss"]
    xi2 = ph * m2 / (2 * m3)
    c0, c1, c2 = 1 / (1 - ph), 3 * xi2 / (1 - ph) ** 2, 2 * xi2 ** 2 / (1 - ph) ** 3
    num = c0 * (2 * m3 * m1 + 2 * m2 ** 2) + 2 * c1 * (2 * m3 * m2) + 4 * c2 * m3 ** 2
    den = c0 * (2 * m2 + 2 * m1 ** 2) + 2 * c1 * (2 * m2 * m1) + 4 * c2 * m2 ** 2
    out["aaBM"] = num / den                              # eq. (12b): hard-sphere-liquid weighting, no parameter
    # reconstructed sigma*/kappa from eq. (7), per snapshot means (not exact for the mean of a ratio)
    out["recon"] = out["ncV"] * out["aa"] * out["C"] / (3 * out["q"] * out["LW"])
    return out


def data_checks() -> None:
    say("\n=== 3. Love-Weber decomposition on the actshape dumps (Derjaguin g, snapshot means) ===")
    recs = load_recs()
    cells = {}
    for key, r in recs.items():
        if "error" in r:
            continue
        path = os.path.join(ROOT, "out", "actshape_dump", key + ".npz")
        if not os.path.exists(path):
            continue
        if not (r["steady"] and r["servo_held"] and r["quasistatic"]):
            continue
        d = np.load(path)
        cells[key] = (r, cell_factors(r, d))
    say(f"  healthy cells with dumps: {len(cells)}")
    say("   mu   p        CV  seed   LW(mean,min,max)        Nc/V      <aa>_c   C       q       s*/k(meas)  s*/k(eq7)")
    rows = sorted(cells.values(), key=lambda t: (t[0]["mu_p"], t[0]["p_hat_target"], t[0]["cv"], t[0]["seed"]))
    for r, f in rows:
        say(f"  {r['mu_p']:.1f} {r['p_hat_target']:.0e} {r['cv']:.1f} {r['seed']:3d}   "
            f"{f['LW']:.4f} [{f['LWmin']:.4f},{f['LWmax']:.4f}]  {f['ncV']:.5f}  {f['aa']:.5f}  "
            f"{f['C']:.4f}  {f['q']:.4f}  {f['sstar']:.5g}  {f['recon']:.5g}")
    # configuration-model test (5b): <r1 r2>_c / <r_end>^2 = 1 iff partner radii are uncorrelated given degree
    say("  configuration-model test <r1 r2>_c/<r_end>^2 by (mu, p, CV): mean [min, max] over seeds")
    cg = {}
    for r, f in rows:
        cg.setdefault((r["mu_p"], r["p_hat_target"], r["cv"]), []).append(f["cfg"])
    for kk in sorted(cg):
        v = np.array(cg[kk])
        say(f"    {kk}: {v.mean():.4f} [{v.min():.4f}, {v.max():.4f}]  n {v.size}")
    # referee M4: count reconstruction exceptions instead of quoting them
    exc = [(r["key"], f["recon"] / f["sstar"] - 1) for r, f in rows if abs(f["recon"] / f["sstar"] - 1) > 0.003]
    say()
    say(f"  eq. 7 reconstruction |recon/meas - 1| > 0.3 %: {len(exc)} of {len(rows)} cells")
    for k, e in exc:
        say(f"    {k}: {e:+.4f}")
    # referee M3: the monodisperse limiting check, printed
    say("  limiting case CV 0 (eq. 13, C = 1): Z phi/(8 pi a) / (q Lambda) vs measured sigma*/kappa")
    for r, f in rows:
        if r["cv"] == 0.0:
            pred = f["Z"] * f["phi"] / (8 * math.pi * r["d10_realised"] / 2) / (f["q"] * f["LW"])
            say(f"    mu {r['mu_p']} p {r['p_hat_target']:.0e} seed {r['seed']}: {pred:.5f} vs {f['sstar']:.5f}  ({pred/f['sstar']-1:+.4f})")
    # referee M1: attrition by condition
    say("  health-filter attrition (attempted -> healthy with dump), by (mu, p, CV):")
    att = {}
    for key, r in recs.items():
        kk = (r.get("mu_p"), r.get("p_hat_target"), r.get("cv"))
        att.setdefault(kk, [0, 0])[0] += 1
    for r, f in rows:
        att[(r["mu_p"], r["p_hat_target"], r["cv"])][1] += 1
    for kk in sorted(att, key=lambda t: tuple(-1 if v is None else v for v in t)):
        say(f"    {kk}: {att[kk][0]} -> {att[kk][1]}")
    # ratios CV -> 0, seed-paired where seeds match
    say("\n  CV ratio of sigma*/kappa split into factors (seed-paired means; lengths rescaled to d10 = 1)")
    say("  sigma*/kappa ~ length^-1 (kappa is force/length), so at fixed d10: multiply by d10(CV)/d10(0)")
    say("   mu   p      CV   n  R_meas   = [Nc/V ratio x <aa> ratio]_geom  x  C ratio  x  1/q ratio   (L^-1 rescale folded into geom)")
    groups = {}
    for r, f in rows:
        groups.setdefault((r["mu_p"], r["p_hat_target"]), {}).setdefault(r["cv"], {})[(r["seed"], f["N"])] = (r, f)
    geo_rows, zrows = [], []
    for (mu, p), byc in sorted(groups.items()):
        if 0.0 not in byc:
            continue
        for cv in (0.2, 0.5):
            if cv not in byc:
                continue
            seeds = sorted(set(byc[0.0]) & set(byc[cv]))   # (seed, N) pairs: same seed AND same N
            if not seeds:
                continue
            Rm, G, Cr, Qr, ZP, AA, MA, MB, RC, FR, ZN, PH, BM = [], [], [], [], [], [], [], [], [], [], [], [], []
            for s in seeds:
                r0, f0 = byc[0.0][s]
                r1_, f1 = byc[cv][s]
                lam = r1_["d10_realised"] / r0["d10_realised"]   # length rescale to d10 = 1
                # a length rescale by 1/lam multiplies sigma*/kappa by lam (dimension L^-1)
                Rm.append(f1["sstar"] / f0["sstar"] * lam)
                # geometric part: Nc/V ~ L^-3, <aa> ~ L^2 -> product ~ L^-1, rescale the same way
                G.append(f1["ncV"] * f1["aa"] / (f0["ncV"] * f0["aa"]) * lam * f0["LW"] / f1["LW"])
                Cr.append(f1["C"] / f0["C"])
                Qr.append(f0["q"] / f1["q"])
                ZP.append(f1["Z"] * f1["phi"] / (f0["Z"] * f0["phi"]))
                AA.append((f1["aa"] / f1["a3"]) / (f0["aa"] / f0["a3"]) * lam)
                MA.append((f1["aaA"] / f1["a3"]) / (f0["aaA"] / f0["a3"]) * lam)
                MB.append((f1["aaB"] / f1["a3"]) / (f0["aaB"] / f0["a3"]) * lam)
                BM.append((f1["aaBM"] / f1["a3"]) / (f0["aaBM"] / f0["a3"]) * lam)
                RC.append(f1["aa"] / f1["aaC"])
                FR.append((1 - r1_["rattler_frac_ss"]) / (1 - r0["rattler_frac_ss"]))
                ZN.append(r1_["z_nr_ss"] / r0["z_nr_ss"])
                PH.append(r1_["phi_ss"] / r0["phi_ss"])
            Rm, G, Cr, Qr = map(np.array, (Rm, G, Cr, Qr))
            geo_rows.append((mu, p, cv, len(seeds), np.mean(ZP), np.mean(AA), np.mean(MA), np.mean(MB), G.mean(), np.mean(RC), np.min(RC), np.max(RC), np.mean(BM)))
            zrows.append((mu, p, cv, len(seeds), np.mean(FR), np.mean(ZN), np.mean(PH), np.mean(ZP), Rm.mean()))
            say(f"  {mu:.1f} {p:.0e}  {cv:.1f}  {len(seeds)}  {Rm.mean():.4f}   = {G.mean():.4f}  x  "
                f"{Cr.mean():.4f}  x  {Qr.mean():.4f}   (product {np.mean(G*Cr*Qr):.4f})")
    say()
    say("  geometric factor split: G = [Z phi ratio] x [(<aa>_c/<a^3>) ratio, d10 = 1]  (eq. 8)")
    say("  pair closures for <aa>_c/<a^3>: A = (a_i+a_j)^2-weighted (eq. 9), B = a_i^2 a_j^2-weighted (eq. 11)")
    say("   mu   p      CV   n  Zphi     <aa>/<a3>(meas)  A        B        | product  G(direct)   (ranges quoted in the .md are ACROSS conditions, not error bars)")
    for mu, p, cv, n, zp, aa, ma, mb, g, rc, rcmin, rcmax, bm in geo_rows:
        say(f"  {mu:.1f} {p:.0e}  {cv:.1f}  {n}  {zp:.4f}   {aa:.4f}           {ma:.4f}   {mb:.4f}   | {zp*aa:.4f}   {g:.4f}"
            f"   C-closure: <aa>_meas/<aa>_C = {rc:.4f} [{rcmin:.4f},{rcmax:.4f}]   BMCSL ratio {bm:.4f}")
    say("  (closure C, eq. 12, predicts (<aa>_c/<a^3>) ratio = 1 exactly at fixed <a>, for ANY PSD)")
    say()
    say("  contact density from the RECORDS: Z = (1 - f_r) Z_nr (eq. 14); ratios CV / CV 0")
    say("   mu   p      CV   n  (1-f_r)   Z_nr     phi      | product  Zphi(dump)  R_meas")
    for mu, p, cv, n, fr, zn, ph, zp, rm in zrows:
        say(f"  {mu:.1f} {p:.0e}  {cv:.1f}  {n}  {fr:.4f}   {zn:.4f}   {ph:.4f}   | {fr*zn*ph:.4f}   {zp:.4f}      {rm:.4f}")


def main() -> int:  # noqa: D103
    sympy_checks()
    mc_check()
    data_checks()
    os.makedirs(os.path.join(ROOT, "out"), exist_ok=True)
    with open(os.path.join(ROOT, "out", "sigma_star_geometry.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

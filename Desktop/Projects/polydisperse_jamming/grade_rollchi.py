"""grade_rollchi.py -- grade PREREG_ROLLCHI_2026-09-25.md. Writes out/grade_rollchi.txt.

Written after the registration and before any graded record was read. Every registered test
(R0-R7, A1-A4) and gate (G0, G1, K4) is computed here, in the registration's words; post-hoc
readings, if any, are labelled POST HOC.

Inputs: out/rollchi/phase2_{rchi,rchi0,achi,achi5,achig}_lognormal*.jsonl (this campaign),
out/shape/phase2_{kmer,shapechi}_*.jsonl (the doublet and alpha = 1 / 2 references).
Seed 99 (the disclosed shakeout) is never read: it lives in out/rollchi_smoke/.
"""
from __future__ import annotations

import glob
import json
import math
import os
from collections import defaultdict

import numpy as np
from scipy import stats

from shape import sub_sphere_diameter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PHI_SCALE = 2.1e-3          # the protocol's own re-run scatter (mixture_model.PHI_SCALE)
SQ3 = 2.0 * math.sqrt(3.0)  # the 2*sqrt(3) of phi = Z/(Z + 2 sqrt 3)
RUNGS = (1e-3, 3e-5)
ETAS = (0.0, 0.01, 0.03, 0.1, 0.3, 1.0)
MUS_R = (0.0, 0.1, 0.3, 1.0, 3.0)
ALPHAS = (1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0)
K4_DEPTH, K4_FRAC = 0.05, 1e-3
L: list[str] = []


def say(s: str = "") -> None:
    print(s)
    L.append(s)


def load(pattern: str) -> list[dict]:
    rows = []
    for f in sorted(glob.glob(pattern)):
        for line in open(f, encoding="utf-8"):
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def tcrit(df: int) -> float:
    return float(stats.t.ppf(0.975, max(df, 1)))


def mean_se(v):
    v = np.asarray(v, float)
    if len(v) == 0:
        return float("nan"), float("nan"), 0
    return float(v.mean()), (float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float("nan")), len(v)


def welch(a, b):
    """(diff, se, t, t_crit) for mean(a) - mean(b), Welch-Satterthwaite dof."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1)) if se > 0 else 1
    d = float(a.mean() - b.mean())
    return d, se, (d / se if se > 0 else float("inf")), tcrit(int(math.floor(df)))


def main() -> int:
    say("grade_rollchi -- PREREG_ROLLCHI_2026-09-25.md")
    raw = load(os.path.join(OUT, "rollchi", "phase2_*_lognormal*.jsonl"))
    say(f"N raw records in out/rollchi: {len(raw)}")
    by_key = {}
    for r in raw:
        by_key[r["key"]] = r                       # a key is written once by construction
    say(f"N unique keys: {len(by_key)} (duplicates dropped: {len(raw) - len(by_key)})")
    recs = list(by_key.values())
    assert all(r.get("seed") != 99 for r in recs), "seed 99 is the disclosed peek"
    failed = [r for r in recs if r.get("failed")]
    notconv = [r for r in recs if not r.get("failed") and not r.get("converged")]
    ok = [r for r in recs if not r.get("failed") and r.get("converged")]
    say(f"G1 HEALTH: {len(ok)} converged, {len(notconv)} not converged (EXCLUDED), {len(failed)} failed (EXCLUDED)")
    for r in notconv + failed:
        say(f"   excluded {r['key']}  {r.get('error', 'not converged')}")
    by_tag = defaultdict(int)
    for r in ok:
        by_tag[r["tag"]] += 1
    say(f"   converged per tag: {dict(sorted(by_tag.items()))}  (registered: rchi 432, rchi0 36, achi 90, achi5 60, achig 1)")

    # ---------------- references -----------------------------------------------------
    ref = [r for r in load(os.path.join(OUT, "shape", "phase2_*_lognormal*.jsonl"))
           if not r.get("failed") and r.get("converged") and r.get("tag") in ("kmer", "shapechi")
           and r.get("n_requested") == 1700 and r.get("cv") == 0.0 and r.get("triplet_w", 0) == 0]

    def refphi(alpha, mu, rung, w2=0.0, tag=None):
        v = {}
        for r in ref:
            if (abs(r["alpha"] - alpha) < 1e-9 and r["mu_p"] == mu and r["p_hat_target"] == rung
                    and abs(r.get("doublet_w", 0) - w2) < 1e-9 and (tag is None or r["tag"] == tag)):
                v[r["seed"]] = r                   # kmer and shapechi alpha-1 records are bit-identical (K0')
        return v

    # ---------------- arm R index --------------------------------------------------------
    R = defaultdict(dict)   # (mu, eta, cv, rung) -> {seed: rec}
    for r in ok:
        if r["tag"].startswith("rchi"):
            R[(r["mu_p"], r["eta_roll"], r["cv"], r["p_hat_target"])][r["seed"]] = r

    def phis(mu, eta, cv, rung):
        d = R.get((mu, eta, cv, rung), {})
        return {s: d[s]["phi"] for s in d}

    # ---------------- R0 ------------------------------------------------------------
    say("\n=== R0 (control, two-sided): mu_p 0 inert, mu_p 0.3 live ===")
    n_pair, n_bit = 0, 0
    for cv in (0.0, 0.2, 0.5):
        for rung in RUNGS:
            a, b = R.get((0.0, 0.0, cv, rung), {}), R.get((0.0, 1.0, cv, rung), {})
            for s in sorted(set(a) & set(b)):
                n_pair += 1
                same = a[s]["phi"] == b[s]["phi"] and a[s]["iters"] == b[s]["iters"]
                n_bit += same
                if not same:
                    say(f"   NOT bit-identical cv {cv} p {rung:.0e} s{s}: phi {a[s]['phi']!r} vs {b[s]['phi']!r}, it {a[s]['iters']} vs {b[s]['iters']}")
    say(f"   inert half: {n_bit}/{n_pair} pairs bit-identical in (phi, iters)  (registered 18)")
    live_ok, live_n = 0, 0
    for cv in (0.0, 0.2, 0.5):
        for rung in RUNGS:
            a, b = phis(0.3, 1.0, cv, rung), phis(0.3, 0.0, cv, rung)
            if a and b:
                live_n += 1
                d = np.mean(list(a.values())) - np.mean(list(b.values()))
                live_ok += abs(d) > 2 * PHI_SCALE
                say(f"   live half cv {cv} p {rung:.0e}: phi(eta 1) - phi(eta 0) = {d:+.4f}  (bar |d| > {2*PHI_SCALE:.4f})")
    r0 = n_pair == 18 and n_bit == n_pair and live_n == 6 and live_ok == live_n
    say(f"R0: {'PASS' if r0 else ('INCOMPLETE' if n_pair < 18 or live_n < 6 else 'FAIL')}")

    # ---------------- R1 ------------------------------------------------------------
    say("\n=== R1 (path control): CohFrictMat eta 0 vs FrictMat records, CV 0, N 1700 ===")
    r1 = []
    for mu in (0.0, 0.3, 3.0):
        for rung in RUNGS:
            a = phis(mu, 0.0, 0.0, rung)
            b = {s: r["phi"] for s, r in refphi(1.0, mu, rung).items()}
            if not a or not b:
                say(f"   mu {mu} p {rung:.0e}: missing (coh n {len(a)}, frict n {len(b)})")
                r1.append(None)
                continue
            d = np.mean(list(a.values())) - np.mean(list(b.values()))
            r1.append(abs(d) < 2 * PHI_SCALE)
            say(f"   mu {mu} p {rung:.0e}: coh {np.mean(list(a.values())):.4f} (n {len(a)}) frict {np.mean(list(b.values())):.4f} (n {len(b)})  d {d:+.4f}  {'ok' if r1[-1] else 'OUT'}")
    say(f"R1: {'PASS' if all(x is True for x in r1) else ('INCOMPLETE' if None in r1 else 'FAIL')}")

    # ---------------- R2 ------------------------------------------------------------
    say("\n=== R2 (saturation at mu_p 0.3, CV 0): |phi(eta) - phi(1.0)| < PHI_SCALE for eta 0.1, 0.3 ===")
    r2 = []
    for rung in RUNGS:
        ref1 = phis(0.3, 1.0, 0.0, rung)
        for eta in ETAS:
            a = phis(0.3, eta, 0.0, rung)
            common = sorted(set(a) & set(ref1))
            if not common:
                continue
            dd = [a[s] - ref1[s] for s in common]         # same seed = same initial state: paired
            m, se, n = mean_se(dd)
            tag = ""
            if eta in (0.1, 0.3):
                r2.append(abs(m) < PHI_SCALE and n == 3)
                tag = "  <- graded " + ("ok" if r2[-1] else "OUT")
            say(f"   p {rung:.0e} eta {eta:<5}: phi {np.mean(list(a.values())):.4f}  paired d vs eta 1 = {m:+.4f} +- {se:.4f} (n {n}){tag}")
    say(f"R2: {'PASS' if len(r2) == 4 and all(r2) else ('INCOMPLETE' if len(r2) < 4 else 'FAIL')}")

    # ---------------- R3 ------------------------------------------------------------
    say("\n=== R3 (H_roll, CONFIRMATORY after the seed-99 peek): saturated rolling sphere vs pure doublet, mu_p 0.3, CV 0 ===")
    r3 = []
    for rung in RUNGS:
        a = list(phis(0.3, 1.0, 0.0, rung).values())
        b = [r["phi"] for r in refphi(1.0, 0.3, rung, w2=1.0, tag="kmer").values()]
        sph = [r["phi"] for r in refphi(1.0, 0.3, rung).values()]
        if len(a) < 2 or len(b) < 2:
            r3.append(None)
            say(f"   p {rung:.0e}: missing")
            continue
        d, se, t, tc = welch(a, b)
        frac = (np.mean(a) - np.mean(sph)) / (np.mean(b) - np.mean(sph))
        r3.append(abs(d) < 2 * PHI_SCALE)
        say(f"   p {rung:.0e}: rolling {np.mean(a):.4f}  doublet {np.mean(b):.4f}  sphere {np.mean(sph):.4f}  "
            f"d {d:+.4f} +- {se:.4f} (t {t:+.2f}, t_crit {tc:.2f})  fraction of doublet shift carried by rolling = {frac:.3f}  "
            f"{'ok' if r3[-1] else 'OUT'} (bar {2*PHI_SCALE:.4f})")
    say(f"R3: {'PASS' if all(x is True for x in r3) else ('INCOMPLETE' if None in r3 else 'FAIL')}")

    # ---------------- R4 ------------------------------------------------------------
    say("\n=== R4 (friction dependence): Delta_roll(mu) = phi(eta 1) - phi(eta 0) monotone non-increasing in mu, CV 0 ===")
    r4 = []
    for rung in RUNGS:
        seq = []
        for mu in MUS_R:
            a, b = phis(mu, 1.0, 0.0, rung), phis(mu, 0.0, 0.0, rung)
            common = sorted(set(a) & set(b))
            dd = [a[s] - b[s] for s in common]
            m, se, n = mean_se(dd)
            seq.append(m)
            say(f"   p {rung:.0e} mu {mu:<4}: Delta_roll {m:+.4f} +- {se:.4f} (n {n})")
        mono = all(seq[i + 1] <= seq[i] for i in range(len(seq) - 1)) and not any(math.isnan(x) for x in seq)
        r4.append(mono and seq[0] == 0.0)
        say(f"   p {rung:.0e}: monotone non-increasing {mono}, Delta_roll(0) == 0 exactly {seq[0] == 0.0}")
    say(f"R4: {'PASS' if len(r4) == 2 and all(r4) else 'FAIL'}")

    # ---------------- R5 / R6 -------------------------------------------------------
    say("\n=== R5 (geometric Z -> phi map): phi = Z/(Z + 2 sqrt 3), Z = z_nr, all converged CV-0 sphere cells of arm R ===")
    res5 = []
    for (mu, eta, cv, rung), d in sorted(R.items()):
        if cv != 0.0:
            continue
        for s, r in d.items():
            res5.append((mu, eta, rung, r["phi"] - r["z_nr"] / (r["z_nr"] + SQ3), r["z_nr"], r["phi"]))
    res = np.array([x[3] for x in res5])
    within = float(np.mean(np.abs(res) < 0.015)) if len(res) else float("nan")
    say(f"   N cells {len(res)}: residual phi - SWM mean {res.mean():+.4f}, sd {res.std(ddof=1):.4f}, "
        f"min {res.min():+.4f}, max {res.max():+.4f}; within 0.015: {within:.3f} (bar 0.90)")
    for key, idx in (("mu", 0), ("eta", 1), ("rung", 2)):
        groups = defaultdict(list)
        for x in res5:
            groups[x[idx]].append(x[3])
        say("   by " + key + ": " + ", ".join(f"{k}: {np.mean(v):+.4f}" for k, v in sorted(groups.items())))
    zz = np.array([x[4] for x in res5])
    say(f"   Z range {zz.min():.2f} - {zz.max():.2f}; spearman(residual, Z) = {stats.spearmanr(res, zz)[0]:+.3f}")
    say(f"R5: {'PASS' if within >= 0.9 else 'FAIL'}")

    say("\n=== R6 (separability): one filling factor w(CV) from mu_p 0 cells, phi = Z/(Z + 2 sqrt3 w), predicts the rest ===")
    r6 = []
    for cv in (0.0, 0.2, 0.5):
        w0 = [r["z_nr"] * (1 - r["phi"]) / (SQ3 * r["phi"])
              for (mu, eta, c, rung), d in R.items() if mu == 0.0 and c == cv for r in d.values()]
        if not w0:
            continue
        w = float(np.mean(w0))
        pr = [(r["phi"] - r["z_nr"] / (r["z_nr"] + SQ3 * w))
              for (mu, eta, c, rung), d in R.items() if mu > 0.0 and c == cv for r in d.values()]
        pr = np.array(pr)
        fr = float(np.mean(np.abs(pr) < 0.015))
        r6.append(fr >= 0.9)
        say(f"   CV {cv}: w = {w:.4f} (sd over {len(w0)} mu-0 cells {np.std(w0, ddof=1):.4f}); N predicted {len(pr)}, "
            f"residual mean {pr.mean():+.4f} sd {pr.std(ddof=1):.4f} range [{pr.min():+.4f}, {pr.max():+.4f}]; within 0.015: {fr:.3f}")
    say(f"R6: {'PASS' if len(r6) == 3 and all(r6) else 'FAIL'}")

    # ---------------- R7 ------------------------------------------------------------
    say("\n=== R7 (size x constraint interaction, the discriminator): rho = D(eta 1)/D(eta 0), D = phi(CV .5) - phi(CV 0), mu_p 0.3 ===")
    for rung in RUNGS:
        D = {}
        for eta in (0.0, 1.0):
            a, b = phis(0.3, eta, 0.5, rung), phis(0.3, eta, 0.0, rung)
            D[eta] = {s: a[s] - b[s] for s in sorted(set(a) & set(b))}
        common = sorted(set(D[0.0]) & set(D[1.0]))
        if len(common) < 2:
            say(f"   p {rung:.0e}: INCOMPLETE (n {len(common)})")
            continue
        d1 = np.array([D[1.0][s] for s in common])
        d0 = np.array([D[0.0][s] for s in common])
        n = len(common)
        m1, m0 = d1.mean(), d0.mean()
        rho = m1 / m0
        c = np.cov(d1, d0, ddof=1) / n
        se = abs(rho) * math.sqrt(c[0, 0] / m1 ** 2 + c[1, 1] / m0 ** 2 - 2 * c[0, 1] / (m1 * m0))
        tc = tcrit(n - 1)
        lo, hi = rho - tc * se, rho + tc * se
        sep = 0.33
        inform = tc * se < sep / 2
        if lo >= 0.45 and hi <= 0.80:
            v = "CLUMP-LIKE"
        elif lo > 0.80:
            v = "NOT CLUMP-LIKE"
        else:
            v = "UNDETERMINED"
        phim = np.mean(list(phis(0.3, 1.0, 0.0, rung).values())) / np.mean(list(phis(0.3, 0.0, 0.0, rung).values()))
        say(f"   p {rung:.0e}: D(eta 0) {m0:+.4f}  D(eta 1) {m1:+.4f}  rho {rho:.3f}, 95% CI [{lo:.3f}, {hi:.3f}] (n {n}, delta method, t_crit {tc:.2f}); "
            f"separable prediction phi_m(1)/phi_m(0) = {phim:.3f}; informative {inform} (half-width {tc*se:.3f} vs 0.165) -> {v if inform else 'NOT INFORMATIVE (' + v + ')'}")

    # ---------------- arm A ---------------------------------------------------------
    say("\n=== ARM A: two-sphere clumps ===")
    A = defaultdict(dict)   # (alpha, mu, cv, rung) -> {seed: rec}
    for r in ok:
        if r["tag"].startswith("achi") and r["tag"] != "achig":
            A[(r["alpha"], r["mu_p"], r["cv"], r["p_hat_target"])][r["seed"]] = r
    for alpha in (1.0, 2.0):
        for mu in (0.0, 0.3, 3.0):
            for rung in RUNGS:
                for s, r in refphi(alpha, mu, rung, tag="shapechi").items():
                    A[(alpha, mu, 0.0, rung)][s] = r
    # G0
    g = [r for r in ok if r["tag"] == "achig"]
    gr = refphi(2.0, 0.3, 1e-3, tag="shapechi").get(1)
    if g and gr:
        same = all(g[0][k] == gr[k] for k in ("phi", "iters", "n_contacts"))
        say(f"G0 wiring (alpha 2, mu .3, p 1e-3, s1 vs shapechi): {'BIT-IDENTICAL' if same else 'DIFFERS'} "
            f"phi {g[0]['phi']!r} vs {gr['phi']!r}, iters {g[0]['iters']} vs {gr['iters']}")
    else:
        say("G0: INCOMPLETE")
    # K4
    k4_bad = set()
    n_k4 = 0
    for (alpha, mu, cv, rung), d in A.items():
        if cv != 0.0 or alpha in (1.0, 2.0):
            continue
        rr = 0.5 * sub_sphere_diameter(1.0, alpha)
        for s, r in d.items():
            p = os.path.join(OUT, "achi_dump", r["key"] + ".npz")
            if not os.path.exists(p):
                continue
            n_k4 += 1
            fn = np.load(p)["fn"]
            frac = float(np.mean(fn / (1.0 * rr ** 2) > K4_DEPTH)) if len(fn) else 0.0
            if frac > K4_FRAC:
                k4_bad.add(r["key"])
                say(f"   K4 EXCLUDED a {alpha} mu {mu} p {rung:.0e} s{s}: {100*frac:.2f}% of contacts deeper than {K4_DEPTH} r")
    say(f"K4: {n_k4} CV-0 clump cells checked, {len(k4_bad)} excluded (CV 0.5 clump cells: member radii not in the dump -> K4 not computable, reported only)")

    def aphi(alpha, mu, cv, rung):
        d = A.get((alpha, mu, cv, rung), {})
        return {s: r["phi"] for s, r in d.items() if r["key"] not in k4_bad}

    say("\n--- phi(alpha) table, CV 0 (mean over seeds, n) ---")
    for rung in RUNGS:
        for mu in (0.0, 0.3, 3.0):
            row = []
            for a in ALPHAS:
                v = aphi(a, mu, 0.0, rung)
                row.append(f"{a}: {np.mean(list(v.values())):.4f}({len(v)})" if v else f"{a}: --")
            say(f"   p {rung:.0e} mu {mu}: " + "  ".join(row))

    a1 = False
    say("\n=== A1 (orientational filling): phi_0(alpha) (mu 0) exceeds the sphere by > 2 PHI_SCALE somewhere in 1.1-1.75, p 3e-5 ===")
    for rung in RUNGS:
        base = np.mean(list(aphi(1.0, 0.0, 0.0, rung).values()))
        vals = {a: np.mean(list(aphi(a, 0.0, 0.0, rung).values())) for a in ALPHAS if aphi(a, 0.0, 0.0, rung)}
        amax = max(vals, key=vals.get)
        say(f"   p {rung:.0e}: sphere {base:.4f}; max {vals[amax]:.4f} at alpha {amax} (excess {vals[amax]-base:+.4f})"
            + ("  <- graded" if rung == 3e-5 else ""))
        if rung == 3e-5:
            a1 = max(v - base for a, v in vals.items() if 1.05 < a < 1.8) > 2 * PHI_SCALE
    say(f"A1: {'PASS' if a1 else 'FAIL'}")

    say("\n=== A2 (where the window opens): chi(alpha) = phi(mu .3)/phi(mu 0) ===")
    chi = {}
    for rung in RUNGS:
        for a in ALPHAS:
            m, f0 = aphi(a, 0.3, 0.0, rung), aphi(a, 0.0, 0.0, rung)
            if m and f0:
                chi[(a, rung)] = np.mean(list(m.values())) / np.mean(list(f0.values()))
        say(f"   p {rung:.0e}: " + "  ".join(f"{a}: {chi[(a, rung)]:.4f}" for a in ALPHAS if (a, rung) in chi))
        mu3 = {a: np.mean(list(aphi(a, 3.0, 0.0, rung).values())) / np.mean(list(aphi(a, 0.0, 0.0, rung).values()))
               for a in ALPHAS if aphi(a, 3.0, 0.0, rung) and aphi(a, 0.0, 0.0, rung)}
        say(f"   p {rung:.0e} (mu_m 3, K4-clean only, secondary): " + "  ".join(f"{a}: {v:.4f}" for a, v in mu3.items()))
    if all((a, 3e-5) in chi for a in (1.0, 1.2, 2.0)):
        line = chi[(1.0, 3e-5)] - 0.5 * (chi[(1.0, 3e-5)] - chi[(2.0, 3e-5)])
        a2 = chi[(1.2, 3e-5)] <= line
        say(f"   graded (p 3e-5): chi(1.2) = {chi[(1.2, 3e-5)]:.4f} vs early-opening line chi(1) - (chi(1)-chi(2))/2 = {line:.4f}")
        say(f"A2: {'PASS (early opening)' if a2 else 'FAIL (late or non-monotone opening)'}")
    else:
        say("A2: INCOMPLETE")

    say("\n=== A3 (eccentricity map): chi_clump(alpha) vs chi_roll(eta_eff(alpha)), alpha >= 1.2, p 3e-5 ===")
    a3 = []
    for rung in RUNGS:
        f0 = np.mean(list(phis(0.0, 0.0, 0.0, rung).values())) if phis(0.0, 0.0, 0.0, rung) else float("nan")
        eta_ax = [e for e in ETAS if e > 0]
        croll = [np.mean(list(phis(0.3, e, 0.0, rung).values())) / f0 for e in eta_ax]
        for a in ALPHAS[2:-1] + (2.0,):
            if (a, rung) not in chi:
                continue
            eta_eff = (math.pi / 4) * (a - 1) * sub_sphere_diameter(1.0, a)
            cr = float(np.interp(math.log(eta_eff), np.log(eta_ax), croll))  # clamps at the ends: eta_eff <= 1 here
            dlt = chi[(a, rung)] - cr
            if rung == 3e-5 and a < 2.0:
                a3.append(abs(dlt) < 0.015)
            say(f"   p {rung:.0e} alpha {a}: eta_eff {eta_eff:.3f}  chi_clump {chi[(a, rung)]:.4f}  chi_roll {cr:.4f}  d {dlt:+.4f}")
    say(f"A3: {'PASS' if a3 and all(a3) else ('INCOMPLETE' if not a3 else 'FAIL')} ({sum(a3)}/{len(a3)})")

    say("\n=== A4 (size x aspect, report only): D(alpha)/D(1), D = phi(CV .5) - phi(CV 0); alpha-1 CV .5 reference = arm R CohFrict eta 0 (R1-path) ===")
    for rung in RUNGS:
        for mu in (0.0, 0.3):
            d1 = np.mean(list(phis(mu, 0.0, 0.5, rung).values())) - np.mean(list(phis(mu, 0.0, 0.0, rung).values())) \
                if phis(mu, 0.0, 0.5, rung) and phis(mu, 0.0, 0.0, rung) else float("nan")
            row = []
            for a in ALPHAS[1:-1]:
                h, l0 = aphi(a, mu, 0.5, rung), aphi(a, mu, 0.0, rung)
                if h and l0:
                    row.append(f"{a}: {(np.mean(list(h.values())) - np.mean(list(l0.values()))) / d1:.3f}")
            say(f"   p {rung:.0e} mu {mu}: D(1) {d1:+.4f}; " + "  ".join(row))

    with open(os.path.join(OUT, "grade_rollchi.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

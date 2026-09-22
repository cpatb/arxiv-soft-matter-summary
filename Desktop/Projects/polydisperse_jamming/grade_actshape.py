"""grade_actshape.py -- grades PREREG_ACTSHAPE_2026-09-22.md. Writes out/grade_actshape.txt.

Written before the cells landed. Reads out/phase9_actshape_*.jsonl and
out/actshape_dump/*.npz, and the rate-matched ladder's records for W1.

    python grade_actshape.py [--out DIR]
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

import actshape as A
from measure import force_moments

HERE = os.path.dirname(os.path.abspath(__file__))
T_CRIT_DF2 = 4.303                       # two-sided 95 %, df = 2 (n = 3 seeds)
NORMS = ("constant", "derjaguin", "area")
FCS = (0.10, 0.25, 0.50, 0.75, 0.90)
PRIMARY = dict(mu=0.3, p=3e-4)
CVS = (0.0, 0.2, 0.5)
MUS = (0.0, 0.3)
PS = (1e-4, 3e-4, 1e-3)


def load_records(out: str, pattern: str) -> dict[str, dict]:
    recs = {}
    for f in sorted(glob.glob(os.path.join(out, pattern))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                recs[r["key"]] = r
    return recs


def snapshots(d) -> list[dict]:
    """Split a dump into per-snapshot arrays."""
    edges = np.concatenate([[0], np.cumsum(d["n"])])
    out = []
    for k in range(len(d["n"])):
        s = slice(edges[k], edges[k + 1])
        out.append({"p": float(d["p_hat"][k]), "fn": d["fn"][s], "r1": d["r1"][s],
                    "r2": d["r2"][s]})
    return out


def cell_metrics(snaps: list[dict], mu_m: float) -> dict:
    m = {}
    for norm in NORMS:
        g = np.concatenate([A.g_values(s["fn"], s["r1"], s["r2"], s["p"], norm)
                            for s in snaps])
        e = {"n": int(g.size),
             "d1090": A.delta_from_quantiles(g, 0.10, 0.90),
             "d2575": A.delta_from_quantiles(g, 0.25, 0.75),
             "sstar": A.sigma_star_over_kappa(g),
             "sstar_tau": A.sigma_star_over_kappa(g) * mu_m}
        for fc in FCS:
            e[f"sjam{fc:.2f}"] = A.sigma_jam_over_kappa(g, fc)
        e["curve_fc"] = np.linspace(0.10, 0.90, 81)
        e["curve"] = np.array([A.sigma_jam_over_kappa(g, fc) for fc in e["curve_fc"]])
        m[norm] = e
    return m


def paired(a: list[float], b: list[float]) -> tuple[float, float, float]:
    d = np.asarray(b) - np.asarray(a)
    if d.size < 2:
        return float(d.mean()) if d.size else float("nan"), float("nan"), float("nan")
    se = float(d.std(ddof=1) / math.sqrt(d.size))
    return float(d.mean()), se, float(d.mean() / se) if se > 0 else float("inf")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "out"))
    ap.add_argument("--report", default=os.path.join(HERE, "out", "grade_actshape.txt"))
    args = ap.parse_args()
    lines: list[str] = []

    def say(s: str = "") -> None:
        print(s)
        lines.append(s)

    recs = load_records(args.out, "phase9_actshape_*.jsonl")
    rm = load_records(os.path.join(HERE, "out"), "phase9_rm_*.jsonl")
    say("grade_actshape -- PREREG_ACTSHAPE_2026-09-22.md")
    say(f"records: {len(recs)} of 54 registered")

    # ---- W4 health + load dumps ------------------------------------------------------
    cells = {}
    n_err = 0
    for key, r in sorted(recs.items()):
        if "error" in r:
            n_err += 1
            say(f"  ERROR CELL {key}: {r['error']}")
            continue
        path = os.path.join(args.out, "actshape_dump", key + ".npz")
        if not os.path.exists(path):
            say(f"  NO DUMP for {key}")
            continue
        d = np.load(path)
        healthy = bool(r["steady"] and r["servo_held"] and r["quasistatic"])
        cells[key] = {"rec": r, "dump": d, "healthy": healthy,
                      "cv": r["cv"], "mu": r["mu_p"], "p": r["p_hat_target"],
                      "seed": r["seed"]}
    bad = [k for k, c in cells.items() if not c["healthy"]]
    say(f"W4 HEALTH: {len(cells) - len(bad)} healthy of {len(cells)} with dumps; "
        f"{n_err} error cells")
    for k in bad:
        say(f"  EXCLUDED {k}: {cells[k]['rec'].get('notes')}")

    # ---- W1 identity against the rate-matched ladder --------------------------------
    say("")
    say("=== W1 IDENTITY vs phase9_rm (dump must not perturb the run) ===")
    w1 = []
    for key, c in sorted(cells.items()):
        if key not in rm:
            continue
        a, b = c["rec"], rm[key]
        same = {f: a[f] == b[f] for f in ("p_hat_series", "z_series", "force")}
        dphi = max(abs(x - y) for x, y in zip(a["phi_series"], b["phi_series"]))
        ok = all(same.values())
        w1.append(ok)
        say(f"  {key[:44]}  p/z/force identical={same}  max|dphi|={dphi:.1e}  "
            f"{'PASS' if ok else 'FAIL'}")
    say(f"W1: {sum(w1)}/{len(w1)} identical (9 registered)"
        + ("  PASS" if w1 and all(w1) and len(w1) == 9 else "  NOT (YET) PASSED"))

    # ---- W2 the dump is the record's network ----------------------------------------
    w2 = []
    for key, c in cells.items():
        last = snapshots(c["dump"])[-1]
        fm = force_moments(last["fn"])
        w2.append(fm == c["rec"]["force"])
    say(f"W2 DUMP == RECORD force block on the last snapshot: {sum(w2)}/{len(w2)}"
        + ("  PASS" if w2 and all(w2) else "  FAIL") + "   [as registered, compared on Windows]")
    # W2' / W2'' added 2026-09-22 after W2 failed on tail_slope alone at ~1e-15: n, mean and
    # cv were bit-identical, and tail_slope is np.polyfit (LAPACK), which differs between the
    # WSL numpy that wrote the record and the Windows one grading it. W2 stays printed.
    w2p = []
    for key, c in cells.items():
        last = snapshots(c["dump"])[-1]
        fm, rf = force_moments(last["fn"]), c["rec"]["force"]
        w2p.append(all(fm[k] == rf[k] for k in ("n", "mean", "cv", "n_tail"))
                   and abs(fm["tail_slope"] - rf["tail_slope"]) < 1e-12)
    say(f"W2' n/mean/cv/n_tail bit-identical and |d tail_slope| < 1e-12: "
        f"{sum(w2p)}/{len(w2p)}" + ("  PASS" if w2p and all(w2p) else "  FAIL"))
    wsl = os.path.join(args.out, "actshape_w2wsl.txt")
    if os.path.exists(wsl):
        with open(wsl, encoding="utf-8") as fh:
            tail = [ln.strip() for ln in fh if ln.startswith("W2''")]
        say(f"{tail[-1] if tail else 'W2 WSL file has no summary line'}"
            "   [recomputed on the platform that wrote the record]")
    else:
        say("W2'' not run (scripts/w2_same_platform.py in WSL)")

    # ---- metrics ------------------------------------------------------------------
    for c in cells.values():
        r = c["rec"]
        mu_m = r["tau_ss_hat"] / r["p_ss_hat"] if r["p_ss_hat"] > 0 else float("nan")
        c["mu_M"] = mu_m
        c["m"] = cell_metrics(snapshots(c["dump"]), mu_m)
        c["d50"] = r["d50_realised"]
        c["d32"] = r["d32_realised"]

    # ---- W3 the anchor that can fail -------------------------------------------------
    w3a = [max(c["m"][n]["d1090"] for n in NORMS) - min(c["m"][n]["d1090"] for n in NORMS)
           for c in cells.values() if c["cv"] == 0.0]
    w3b = [abs(c["m"]["constant"]["d1090"] - c["m"]["derjaguin"]["d1090"])
           for c in cells.values() if c["cv"] == 0.5]
    ok3 = bool(w3a) and max(w3a) < 1e-12 and bool(w3b) and min(w3b) > 1e-3
    say(f"W3 CV=0 norms identical: worst spread {max(w3a) if w3a else float('nan'):.1e} "
        f"(bar 1e-12); CV=0.5 constant vs derjaguin differ: min "
        f"{min(w3b) if w3b else float('nan'):.3e} (bar 1e-3)  {'PASS' if ok3 else 'FAIL'}")

    good = {k: c for k, c in cells.items() if c["healthy"]}
    grp: dict[tuple, dict[int, dict]] = defaultdict(dict)
    for c in good.values():
        grp[(c["cv"], c["mu"], c["p"])][c["seed"]] = c

    def col(cv, mu, p, norm, field):
        g = grp.get((cv, mu, p), {})
        return {s: c["m"][norm][field] for s, c in g.items()}

    # ---- the tables ---------------------------------------------------------------
    for norm in NORMS:
        say("")
        say(f"=== norm = {norm}   (seed mean +- sd; n = seeds) ===")
        say("   mu     p      CV   n   d10/90          d25/75          sigma*/kappa(d10)   "
            "mu_M")
        for mu in MUS:
            for p in PS:
                for cv in CVS:
                    g = grp.get((cv, mu, p), {})
                    if not g:
                        continue
                    vals = {f: np.array([c["m"][norm][f] for c in g.values()])
                            for f in ("d1090", "d2575", "sstar")}
                    mum = np.array([c["mu_M"] for c in g.values()])
                    say(f"  {mu:.1f}  {p:.0e}  {cv:.1f}  {len(g)}   "
                        f"{vals['d1090'].mean():.3f} +- {vals['d1090'].std(ddof=1) if len(g) > 1 else 0:.3f}   "
                        f"{vals['d2575'].mean():.3f} +- {vals['d2575'].std(ddof=1) if len(g) > 1 else 0:.3f}   "
                        f"{vals['sstar'].mean():.4g} +- {vals['sstar'].std(ddof=1) if len(g) > 1 else 0:.2g}   "
                        f"{mum.mean():.3f}")

    # ---- hypotheses -----------------------------------------------------------------
    mu0, p0 = PRIMARY["mu"], PRIMARY["p"]

    def contrast(field, norm, cv_a, cv_b, mu, p):
        a, b = col(cv_a, mu, p, norm, field), col(cv_b, mu, p, norm, field)
        seeds = sorted(set(a) & set(b))
        return [a[s] for s in seeds], [b[s] for s in seeds], seeds

    say("")
    say("=== H1  delta_10/90 (derjaguin), CV 0 -> 0.5, mu 0.3, p 3e-4, seed-paired ===")
    a, b, s = contrast("d1090", "derjaguin", 0.0, 0.5, mu0, p0)
    if len(s) >= 2:
        mean, se, t = paired(a, b)
        ratio = float(np.mean(b) / np.mean(a))
        h1 = mean < 0 and abs(t) >= T_CRIT_DF2 and 0.60 <= ratio <= 0.90
        say(f"  seeds {s}: diff {mean:+.4f} +- {se:.4f}, t = {t:.2f} (t_crit {T_CRIT_DF2}); "
            f"ratio {ratio:.3f} (registered [0.60, 0.90])  -> {'PASS' if h1 else 'FAIL'}")
        if mean < 0 and abs(t) >= T_CRIT_DF2 and not 0.60 <= ratio <= 0.90:
            say("  (sign and significance as predicted; the MAGNITUDE band is missed)")
    else:
        say("  not gradable yet")

    say("")
    say("=== H1b  ordering at CV 0.5: delta(constant) < delta(derjaguin) < delta(area) ===")
    for mu in MUS:
        for p in PS:
            g = grp.get((0.5, mu, p), {})
            if not g:
                continue
            m = {n: np.mean([c["m"][n]["d1090"] for c in g.values()]) for n in NORMS}
            okb = m["constant"] < m["derjaguin"] < m["area"]
            say(f"  mu {mu:.1f} p {p:.0e}: {m['constant']:.3f} / {m['derjaguin']:.3f} / "
                f"{m['area']:.3f}  {'PASS' if okb else 'FAIL'}")

    say("")
    say("=== H2  sigma*/kappa (derjaguin), CV 0 -> 0.5, at fixed d10 / d50 / d32 ===")
    for norm, mexp in (("derjaguin", 1), ("constant", 0)):
        for mu in MUS:
            a, b, s = contrast("sstar", norm, 0.0, 0.5, mu, p0)
            if len(s) < 2:
                continue
            ga, gb = grp[(0.0, mu, p0)], grp[(0.5, mu, p0)]
            rat10 = [bb / aa for aa, bb in zip(a, b)]
            # rescaling lengths so d50 (d32) = 1 multiplies sigma* by d50^(2-m)
            rat50 = [gb[x]["m"][norm]["sstar"] * gb[x]["d50"] ** (2 - mexp)
                     / (ga[x]["m"][norm]["sstar"] * ga[x]["d50"] ** (2 - mexp)) for x in s]
            rat32 = [gb[x]["m"][norm]["sstar"] * gb[x]["d32"] ** (2 - mexp)
                     / (ga[x]["m"][norm]["sstar"] * ga[x]["d32"] ** (2 - mexp)) for x in s]
            mean, se, t = paired(list(np.log(a)), list(np.log(b)))
            graded = norm == "derjaguin"
            verdict = ("PASS" if (mean < 0 and abs(t) >= T_CRIT_DF2) else "FAIL") \
                if graded else "reported"
            say(f"  {norm:9s} mu {mu:.1f}: sigma*(0.5)/sigma*(0) at d10 {np.mean(rat10):.3f} "
                f"(ln diff t = {t:.2f}), at d50 {np.mean(rat50):.3f}, at d32 {np.mean(rat32):.3f}"
                f"   [{verdict if mu == mu0 else 'reported'}]")

    say("")
    say("=== H3  friction: delta_10/90 (derjaguin) mu 0 vs 0.3, p 3e-4, seed-paired ===")
    for cv in CVS:
        a = col(cv, 0.0, p0, "derjaguin", "d1090")
        b = col(cv, 0.3, p0, "derjaguin", "d1090")
        s = sorted(set(a) & set(b))
        if len(s) >= 2:
            mean, se, t = paired([a[x] for x in s], [b[x] for x in s])
            say(f"  CV {cv:.1f}: mu0 {np.mean([a[x] for x in s]):.3f}  mu0.3 "
                f"{np.mean([b[x] for x in s]):.3f}  diff {mean:+.3f} +- {se:.3f} t {t:.2f}")

    say("")
    say("=== H4  is delta a material property? span over p < half the CV 0->0.5 effect ===")
    for norm in ("derjaguin", "constant"):
        for mu in MUS:
            eff = [np.mean(list(col(0.5, mu, p, norm, "d1090").values()))
                   - np.mean(list(col(0.0, mu, p, norm, "d1090").values())) for p in PS
                   if col(0.5, mu, p, norm, "d1090") and col(0.0, mu, p, norm, "d1090")]
            eff_primary = (np.mean(list(col(0.5, mu, p0, norm, "d1090").values()))
                           - np.mean(list(col(0.0, mu, p0, norm, "d1090").values()))
                           if col(0.5, mu, p0, norm, "d1090") and col(0.0, mu, p0, norm, "d1090")
                           else float("nan"))
            for cv in CVS:
                vals = [np.mean(list(col(cv, mu, p, norm, "d1090").values()))
                        for p in PS if col(cv, mu, p, norm, "d1090")]
                if len(vals) < 3:
                    continue
                span = max(vals) - min(vals)
                okh = span < 0.5 * abs(eff_primary)
                tag = ("PASS" if okh else "FAIL") if norm == "derjaguin" else "reported"
                say(f"  {norm:9s} mu {mu:.1f} CV {cv:.1f}: delta over p "
                    f"{' / '.join(f'{v:.3f}' for v in vals)}  span {span:.3f} vs half-effect "
                    f"{0.5 * abs(eff_primary):.3f}  [{tag}]")
            if eff:
                say(f"      CV effect at each p: {' / '.join(f'{e:+.3f}' for e in eff)}")

    say("")
    say("=== H5  Weibull adequacy: delta_25/75 > delta_10/90 on CV 0 ===")
    for mu in MUS:
        for p in PS:
            g = grp.get((0.0, mu, p), {})
            if not g:
                continue
            a = np.mean([c["m"]["derjaguin"]["d1090"] for c in g.values()])
            b = np.mean([c["m"]["derjaguin"]["d2575"] for c in g.values()])
            say(f"  mu {mu:.1f} p {p:.0e}: d10/90 {a:.3f}  d25/75 {b:.3f}  ratio {b / a:.3f}  "
                f"{'PASS' if b > a else 'FAIL'}")

    say("")
    say("=== H6  sigma_jam/kappa vs f_c, CV 0 vs 0.5, derjaguin, fixed d10: do they cross? ===")
    for mu in MUS:
        g0, g5 = grp.get((0.0, mu, p0), {}), grp.get((0.5, mu, p0), {})
        if not (g0 and g5):
            continue
        c0 = np.mean([np.log(c["m"]["derjaguin"]["curve"]) for c in g0.values()], axis=0)
        c5 = np.mean([np.log(c["m"]["derjaguin"]["curve"]) for c in g5.values()], axis=0)
        fcg = next(iter(g0.values()))["m"]["derjaguin"]["curve_fc"]
        diff = c5 - c0
        cross = bool(np.any(np.sign(diff[:-1]) != np.sign(diff[1:])))
        where = fcg[np.flatnonzero(np.sign(diff[:-1]) != np.sign(diff[1:]))]
        say(f"  mu {mu:.1f}: ln[sjam(0.5)/sjam(0)] at f_c 0.1/0.25/0.5/0.75/0.9 = "
            + " / ".join(f"{np.interp(f, fcg, diff):+.3f}" for f in FCS)
            + f"   cross={cross} at f_c {np.round(where, 3).tolist()}"
            + ("   [H6 " + ("PASS" if cross else "FAIL") + "]" if mu == mu0 else ""))

    say("")
    say("solver presets for scale (NOT graded): DELTA si_fit 1.3852 (withdrawn), cs_fit 1.3169")
    with open(args.report, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""grade_actshape_fs.py -- PREREG_ACTSHAPE Amendment 1, the N = 4000 finite-size control.
Writes out/grade_actshape_fs.txt.

F1  |delta(N=4000) - delta(N=1000)| < 0.05 at CV 0.5, mu_p 0.3 (unpaired seed means)
F2  H1 redone: CV 0.5 at N = 4000 vs CV 0 at N = 1000, Welch t, delta falls, |t| >= t_crit,
    ratio in [0.60, 0.90]
F3  sigma*/kappa (derjaguin, d10) ratio CV 0.5 (N 4000) / CV 0 (N 1000) < 1
Every verdict is printed twice: all healthy N = 4000 seeds, and only the seeds whose
d_max/L clears 0.25 (seed 2 is 0.269, still flagged).
"""
import os
import sys

import numpy as np
from scipy import stats

import grade_actshape as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def cells(pattern, cv, mu, n):
    res = {}
    for key, r in G.load_records(OUT, pattern).items():
        if "error" in r or r["cv"] != cv or r["mu_p"] != mu or r["p_hat_target"] != 3e-4:
            continue
        if r["n_requested"] != n or not (r["steady"] and r["servo_held"] and r["quasistatic"]):
            continue
        d = np.load(os.path.join(OUT, "actshape_dump", key + ".npz"))
        m = G.cell_metrics(G.snapshots(d), r["tau_ss_hat"] / r["p_ss_hat"])["derjaguin"]
        res[r["seed"]] = {"d": m["d1090"], "s": m["sstar"], "fs": r["d_max_over_L"] > 0.25}
    return res


def welch(a, b):
    t, p = stats.ttest_ind(b, a, equal_var=False)
    va, vb = np.var(a, ddof=1) / len(a), np.var(b, ddof=1) / len(b)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    return float(t), float(df), float(stats.t.ppf(0.975, df))


def main():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    c0 = cells("phase9_actshape_c*.jsonl", 0.0, 0.3, 1000)
    c5 = cells("phase9_actshape_c*.jsonl", 0.5, 0.3, 1000)
    c4k = cells("phase9_actshape_fs*.jsonl", 0.5, 0.3, 4000)
    m0 = cells("phase9_actshape_fs*.jsonl", 0.5, 0.0, 4000)
    say("grade_actshape_fs -- PREREG_ACTSHAPE Amendment 1 (mu_p 0.3, p 3e-4, derjaguin)")
    say(f"N=1000 CV0 delta {[round(v['d'], 4) for v in c0.values()]}")
    say(f"N=1000 CV0.5 delta {[round(v['d'], 4) for v in c5.values()]} (all finite-size flagged)")
    say(f"N=4000 CV0.5 delta {[(s, round(v['d'], 4), 'FLAG' if v['fs'] else 'ok') for s, v in sorted(c4k.items())]}")
    for label, sel in (("all healthy N=4000 seeds", lambda v: True),
                       ("unflagged N=4000 seeds only", lambda v: not v["fs"])):
        k = [v for v in c4k.values() if sel(v)]
        say("")
        say(f"--- {label}: n = {len(k)} ---")
        if len(k) < 2:
            say("   fewer than 2 seeds: F1-F3 not gradable on this subset")
            continue
        d4 = [v["d"] for v in k]
        d1 = [v["d"] for v in c5.values()]
        f1 = abs(np.mean(d4) - np.mean(d1))
        say(f"F1 |delta(4000) - delta(1000)| = |{np.mean(d4):.4f} - {np.mean(d1):.4f}| = {f1:.4f}"
            f" (bar 0.05)  {'PASS' if f1 < 0.05 else 'FAIL'}")
        dz = [v["d"] for v in c0.values()]
        t, df, tc = welch(dz, d4)
        ratio = np.mean(d4) / np.mean(dz)
        ok2 = np.mean(d4) < np.mean(dz) and abs(t) >= tc and 0.60 <= ratio <= 0.90
        say(f"F2 CV0 (N1000) -> CV0.5 (N4000): ratio {ratio:.3f}, Welch t {t:.2f}, df {df:.2f}, "
            f"t_crit {tc:.3f}  {'PASS' if ok2 else 'FAIL'}")
        s4 = np.mean([v["s"] for v in k]) / np.mean([v["s"] for v in c0.values()])
        say(f"F3 sigma*/kappa ratio CV0.5(N4000)/CV0(N1000) at d10 = {s4:.3f}  "
            f"{'PASS' if s4 < 1 else 'FAIL'}")
    say("")
    say(f"mu_p 0 at N=4000 (printed, not graded): healthy seeds {sorted(m0)}, delta "
        f"{[round(v['d'], 4) for v in m0.values()]}")
    with open(os.path.join(OUT, "grade_actshape_fs.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""actshape_cv02.py -- the CV 0 -> 0.2 activation-shape contrast, the part of PREREG_ACTSHAPE
that carries NO finite-size flag (every CV 0.5 cell is flagged d_max/L > 0.25; Amendment 1).

Reported, not a registered gate: H1 was registered for CV 0 -> 0.5. This prints the same
estimators, seed-paired, for CV 0 -> 0.2 at every rung and both frictions, so the clean part
of the result can be quoted while the N = 4000 control runs. Writes out/actshape_cv02.txt.
"""
import os
import sys

import numpy as np

import grade_actshape as G

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> int:
    recs = G.load_records(os.path.join(HERE, "out"), "phase9_actshape_*.jsonl")
    cells = {}
    for key, r in recs.items():
        if "error" in r or r.get("n_requested") != 1000 or r["cv"] not in (0.0, 0.2):
            continue
        if not (r["steady"] and r["servo_held"] and r["quasistatic"]):
            continue
        if any("FINITE" in n for n in r["notes"]):
            print("unexpected finite-size flag", key)
            return 1
        d = np.load(os.path.join(HERE, "out", "actshape_dump", key + ".npz"))
        mu_m = r["tau_ss_hat"] / r["p_ss_hat"]
        cells[(r["cv"], r["mu_p"], r["p_hat_target"], r["seed"])] = G.cell_metrics(
            G.snapshots(d), mu_m)
    lines = ["actshape CV 0 -> 0.2 (no finite-size flag), seed-paired, t_crit(df 2) = 4.303",
             "  mu   p       norm       field   CV0      CV0.2    diff     se      t     ratio"]
    for mu in G.MUS:
        for p in G.PS:
            for norm in ("derjaguin", "constant"):
                for field in ("d1090", "sstar"):
                    s = sorted(k[3] for k in cells if k[:3] == (0.0, mu, p)
                               and (0.2, mu, p, k[3]) in cells)
                    if len(s) < 2:
                        continue
                    a = [cells[(0.0, mu, p, x)][norm][field] for x in s]
                    b = [cells[(0.2, mu, p, x)][norm][field] for x in s]
                    if field == "sstar":
                        m, se, t = G.paired(list(np.log(a)), list(np.log(b)))
                    else:
                        m, se, t = G.paired(a, b)
                    lines.append(f"  {mu:.1f}  {p:.0e}  {norm:9s}  {field:6s}  {np.mean(a):.4f}  "
                                 f"{np.mean(b):.4f}  {m:+.4f}  {se:.4f}  {t:+6.2f}  "
                                 f"{np.mean(b) / np.mean(a):.3f}  n={len(s)}"
                                 + ("  (ln)" if field == "sstar" else ""))
    print("\n".join(lines))
    with open(os.path.join(HERE, "out", "actshape_cv02.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

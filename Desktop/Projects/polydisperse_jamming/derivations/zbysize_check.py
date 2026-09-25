"""zbysize_check.py -- direct per-grain test of rattler_closure_2026-09-24.md eq. 3.

    python derivations/zbysize_check.py  ->  out/zbysize_check.txt

For each out/zbysize/*.npz (per grain: diameter d, contact count k, iterated rattler flag), in 8 quantile
bins of d: measured <k | d> and P(rattler | d) against three shapes of lambda(a), each normalised so that
its PSD mean equals the measured Z (so the only thing tested is the SIZE DEPENDENCE):
   C:  lambda ~ a (a m1 + m2)   (closure C, eq. 3)
   a1: lambda ~ a
   a2: lambda ~ a^2
The shape comparison uses <k|d>/Z (no Poisson assumption); the rattler comparison additionally assumes
Poisson counts. A weighted rms over bins (weights = grains per bin) summarises each.
"""
from __future__ import annotations

import glob
import os

import numpy as np
from scipy.stats import poisson

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def shapes(a: np.ndarray) -> dict[str, np.ndarray]:
    m1, m2 = a.mean(), (a * a).mean()
    raw = {"C": a * (a * m1 + m2), "a1": a, "a2": a * a}
    return {k: v / v.mean() for k, v in raw.items()}


def main() -> int:
    files = sorted(glob.glob(os.path.join(ROOT, "out", "zbysize", "*.npz")))
    say(f"zbysize_check -- {len(files)} cells")
    for f in files:
        z = np.load(f)
        d, k, rat = z["d"], z["k"].astype(float), z["rattler"]
        a = d / 2
        Z = k.mean()
        sh = shapes(a)
        say(f"\n{os.path.basename(f)}: N {d.size}, phi {float(z['phi']):.6f}, converged {bool(z['converged'])}, "
            f"Z {Z:.3f}, f_r {rat.mean():.3f}")
        edges = np.unique(np.quantile(d, np.linspace(0, 1, 9)))
        idx = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, len(edges) - 2)
        say("   d/<d>   n    <k|d>   C      a1     a2   |  P(rat|d)  C      a1     a2")
        err = {m: [] for m in sh}
        errr = {m: [] for m in sh}
        w = []
        for b in range(len(edges) - 1):
            m = idx == b
            if not m.any():
                continue
            kk = k[m].mean()
            pr = rat[m].mean()
            pred = {mm: float(Z * sh[mm][m].mean()) for mm in sh}
            prr = {mm: float(np.mean(poisson.cdf(1, Z * sh[mm][m]))) for mm in sh}
            for mm in sh:
                err[mm].append(pred[mm] - kk)
                errr[mm].append(prr[mm] - pr)
            w.append(m.sum())
            say(f"   {d[m].mean() / d.mean():.3f}  {m.sum():4d}  {kk:.3f}  {pred['C']:.3f}  {pred['a1']:.3f}  "
                f"{pred['a2']:.3f} |  {pr:.3f}     {prr['C']:.3f}  {prr['a1']:.3f}  {prr['a2']:.3f}")
        w = np.array(w, float)
        for mm in sh:
            e, er = np.array(err[mm]), np.array(errr[mm])
            say(f"   {mm:3s}: <k|d> weighted rms err {np.sqrt(np.sum(w * e**2) / w.sum()):.3f};  "
                f"P(rat|d) weighted rms err {np.sqrt(np.sum(w * er**2) / w.sum()):.3f};  "
                f"f_r pred {np.mean(poisson.cdf(1, Z * sh[mm])):.3f} (meas {rat.mean():.3f})")
    with open(os.path.join(ROOT, "out", "zbysize_check.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

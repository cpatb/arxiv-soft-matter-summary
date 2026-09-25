"""degree_exponent_probe.py -- EXPLORATORY: which power p in <k|a> ~ a^p fits the per-grain MEAN degree
(out/zbysize, isotropic), and what configuration-model <a_i a_j>_c = (m_{p+1}/m_p)^2 it implies compared with
the measured one. Hypothesis-generating for sigma_star §5c; not a test of the SHEARED degree law.

    python derivations/degree_exponent_probe.py  ->  out/degree_exponent_probe.txt
"""
from __future__ import annotations

import glob
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def main() -> int:
    ps = np.linspace(0.0, 3.0, 61)
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "zbysize", "*.npz"))):
        z = np.load(f)
        a, k = z["d"] / 2, z["k"].astype(float)
        Z = k.mean()
        rms = [np.sqrt(np.mean((Z * a**p / np.mean(a**p) - k) ** 2)) for p in ps]
        pb = ps[int(np.argmin(rms))]
        ak_meas = np.sum(k * a) / np.sum(k)
        ak_p = np.mean(a ** (pb + 1)) / np.mean(a ** pb)
        ak_2 = np.mean(a**3) / np.mean(a**2)
        say(f"{os.path.basename(f)}: best p (per-grain rms of k) {pb:.2f} (rms {min(rms):.3f}; p=2 rms "
            f"{rms[int(np.argmin(abs(ps - 2)))]:.3f})  a_k: measured {ak_meas:.4f}, p_best {ak_p:.4f}, p=2 {ak_2:.4f}")
    with open(os.path.join(ROOT, "out", "degree_exponent_probe.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

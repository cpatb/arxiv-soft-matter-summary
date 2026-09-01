"""calib_cost.py — measure per-cell wall clock and probe the Phase-1 anchors.

Sizes the Phase-2/3 sweeps from a MEASURED cost, and checks the frictionless and
high-friction coordination anchors before any production sweep is launched.

N is constrained from below by the periodic collider: the largest grain needs a final
cell edge >= 2.5*d_max, and at CV=0.5 with truncate=(0.2,5.0) that means N >= 1299.
See `packing.build_and_compress`.
"""
import sys
import time

from packing import build_and_compress


def main() -> int:
    print("=== cost vs N (cv=0.3, mu=0.3, p_hat=1e-3) ===", flush=True)
    for n in (500, 1000, 1500):
        build_and_compress(cv=0.3, mu_p=0.3, p_hat=1e-3, n=n, seed=1)

    print("\n=== cost vs cv (N=1500, mu=0.3, p_hat=1e-3) ===", flush=True)
    for cv in (0.0, 0.3, 0.5):
        build_and_compress(cv=cv, mu_p=0.3, p_hat=1e-3, n=1500, seed=1)

    print("\n=== ANCHOR A1/A2: frictionless monodisperse, p_hat ladder ===", flush=True)
    print("A1 expects phi -> ~0.64 as p_hat -> 0 ; A2 expects Z_nr -> 6", flush=True)
    for p_hat in (1e-2, 1e-3, 1e-4, 1e-5):
        build_and_compress(cv=0.0, mu_p=0.0, p_hat=p_hat, n=1000, seed=1)

    print("\n=== ANCHOR A3/A4: friction ladder at p_hat=1e-4 ===", flush=True)
    print("A3 expects Z_nr -> 4 at large mu ; A4 expects phi to FALL with mu", flush=True)
    for mu in (0.0, 0.1, 0.3, 1.0, 10.0):
        build_and_compress(cv=0.0, mu_p=mu, p_hat=1e-4, n=1000, seed=1)
    return 0


if __name__ == "__main__":
    t = time.time()
    rc = main()
    print(f"\ntotal wall clock: {time.time() - t:.1f} s")
    sys.exit(rc)

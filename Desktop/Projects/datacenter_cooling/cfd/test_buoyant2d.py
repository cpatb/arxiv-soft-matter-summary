"""Validate the buoyant solver against the de Vahl Davis (1983) benchmark:
average Nusselt number on the hot wall of a differentially heated square cavity."""
import numpy as np
from buoyant2d import Cavity

# de Vahl Davis (1983) benchmark, Pr=0.71, average hot-wall Nu
BENCH = {1e3: (40, 1.118), 1e4: (56, 2.243)}   # (grid N, reference Nu)
TOL = 0.10   # 10% on a coarse grid


def test_nusselt_benchmarks():
    ok = True
    for Ra, (N, ref) in BENCH.items():
        c = Cavity(N=N, H=1.0, Ra=Ra, Pr=0.71)
        c.solve()
        nu = c.nusselt_hot()
        err = abs(nu - ref) / ref
        flag = "PASS" if err < TOL else "FAIL"
        print(f"Ra={Ra:.0e}: Nu={nu:.3f} (ref {ref:.3f}), err {100*err:4.1f}%  [{flag}]")
        ok = ok and err < TOL
    assert ok, "one or more Nusselt benchmarks outside tolerance"


def test_no_convection_gives_conduction():
    # Ra->0: pure conduction, Nu -> 1.0
    c = Cavity(N=32, H=1.0, Ra=1.0, Pr=0.71)
    c.solve()
    assert abs(c.nusselt_hot() - 1.0) < 0.03, "conduction limit Nu must be ~1"


if __name__ == "__main__":
    for fn in (test_no_convection_gives_conduction, test_nusselt_benchmarks):
        try:
            fn(); print(f"[PASS] {fn.__name__}")
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")

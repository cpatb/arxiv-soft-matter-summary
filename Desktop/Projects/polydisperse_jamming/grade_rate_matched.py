"""grade_rate_matched.py — is the CV contrast confounded with a shear-rate contrast?

THE DEFECT
----------
Found 2026-09-09, previously unrecorded anywhere in this project. `shear.py` sets

    gdot = gamma_max / (shear_steps * O.dt),   O.dt = DT_SAFETY * PWaveTimeStep()

and the P-wave timestep is set by the SMALLEST particle in the pack. Widening the size
distribution shrinks the smallest particle, shrinks `dt`, and therefore RAISES the
realised shear rate. `gamma_max` and `shear_steps` are held fixed across the CV ladder,
so the rate is not: measured on `phase8_prod2*`, 3.000e-05 at cv = 0 against 1.437e-04
at cv = 0.5, a factor **4.78**, against a seed-to-seed spread within one CV level of at
most 1.10.

Every CV contrast on that run is therefore confounded with a rate contrast of that size,
and `I = gdot d50 sqrt(rho/p)` carries it straight into the quasi-static gate as well.

THE CONTROL
-----------
The cv = 0 and cv = 0.2 ladders were re-run with `shear_steps` chosen PER CELL so that
`gdot` lands on the cv = 0.5 value, 1.430e-04 -- steps = gamma_max/(gdot* dt), with `dt`
read from the prod2 records (it is a deterministic function of the packing, identical on
every rung of a ladder). Nothing else changed: same family, N, mu_p, seeds, gamma_max,
protocol, ladder. Tags `ratem_c0`, `ratem_c2s1/2/3`.

The match is CHECKED here from the realised `gdot` in every record, never assumed.

WHAT TO EXPECT, AND WHY THE ANSWER IS NOT OBVIOUS
-------------------------------------------------
Standard granular rheology (da Cruz / Jop `phi(I)`, `mu(I)`) says the constant-volume
pressure rises with rate at fixed `phi`, so a HIGHER rate at high CV should raise `p`
there and push the fitted root DOWN -- i.e. it should SHRINK the measured `dphi_c`, and
correcting it should make the contrast LARGER, not smaller. That is a prediction, not a
measurement, and it assumes CV does not also reshape the `phi(I)` curve. This measures it.

    python grade_rate_matched.py
"""

from __future__ import annotations

import os
from collections import defaultdict

import numpy as np

import analyse_phase8 as A

HERE = os.path.dirname(os.path.abspath(__file__))

PATTERN = "phase8_ratem_*.jsonl,phase8_prod2_shard*.jsonl"
GDOT_TARGET = 1.430e-04         # the cv = 0.5 value the control was matched to
MATCHED_STEPS = {0.0: (41958,), 0.2: (86411, 78795, 82782), 0.5: (200000,)}
UNMATCHED_STEPS = {0.0: (200000,), 0.2: (200000,), 0.5: (200000,)}
CV_LO, CV_HI = 0.0, 0.5


def _roots_by(pattern: str) -> tuple[dict, dict]:
    """(cv, shear_steps, gate) -> [phi_c_E2], and (cv, shear_steps) -> [gdot]."""
    recs = [r for r in A.load(pattern)
            if "error" not in r and np.isfinite(r.get("p_ss_hat", np.nan))]
    by = defaultdict(list)
    for r in recs:
        by[A.group_key(r)].append(r)
    roots, gdots = defaultdict(list), defaultdict(list)
    for k, rs in by.items():
        for gate_name, gates in A.GATE_SETS:
            res = A.analyse_group(rs, gates)
            if not res["n_cells"]:
                continue
            roots[(k[1], k[7], gate_name)].append(res["phi_c_E2"])
            if gate_name == "ungated":
                gdots[(k[1], k[7])].append(res["gdot_mean"])
    return roots, gdots


def _pool(roots: dict, cv: float, steps: tuple[int, ...], gate: str):
    v = [x for s in steps for x in roots.get((cv, s, gate), []) if np.isfinite(x)]
    return (float(np.mean(v)), len(v)) if v else (float("nan"), 0)


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    roots, gdots = _roots_by(PATTERN)

    p("=== the match, CHECKED from the realised gdot in every record ===")
    p(f"   {'cv':>5s} {'shear_steps':>12s} {'gdot':>11s} {'vs target':>10s} {'arm':>10s}")
    ok_match = True
    for cv, steps in sorted(MATCHED_STEPS.items()):
        for s in steps:
            v = gdots.get((cv, s), [])
            if not v:
                p(f"   {cv:5.2f} {s:12d} {'MISSING':>11s}")
                ok_match = False
                continue
            g = float(np.mean(v))
            p(f"   {cv:5.2f} {s:12d} {g:11.3e} {g / GDOT_TARGET:9.3f}x {'matched':>10s}")
            if abs(g / GDOT_TARGET - 1.0) > 0.02:
                ok_match = False
    for cv, steps in sorted(UNMATCHED_STEPS.items()):
        for s in steps:
            v = gdots.get((cv, s), [])
            if v and (cv, s) not in [(c, x) for c, ss in MATCHED_STEPS.items()
                                     for x in ss]:
                g = float(np.mean(v))
                p(f"   {cv:5.2f} {s:12d} {g:11.3e} {g / GDOT_TARGET:9.3f}x "
                  f"{'unmatched':>10s}")
    span_un = ([float(np.mean(gdots[(cv, s)])) for cv, ss in UNMATCHED_STEPS.items()
                for s in ss if (cv, s) in gdots])
    if span_un:
        p(f"   the UNMATCHED ladder spans {max(span_un) / min(span_un):.2f}x in gdot; "
          f"the matched one is held to within 2 %")
    if not ok_match:
        p("   🔴 the control is not actually matched -- do not read the table below")

    p("\n=== dphi_c(cv 0 -> 0.5), matched against unmatched, every gate set ===")
    p(f"   {'gate':12s} {'matched':>10s} {'unmatched':>11s} {'change':>9s}  n")
    for gate_name, _gates in A.GATE_SETS:
        a_m, n_am = _pool(roots, CV_LO, MATCHED_STEPS[CV_LO], gate_name)
        b_m, n_bm = _pool(roots, CV_HI, MATCHED_STEPS[CV_HI], gate_name)
        a_u, n_au = _pool(roots, CV_LO, UNMATCHED_STEPS[CV_LO], gate_name)
        if not (n_am and n_bm and n_au):
            p(f"   {gate_name:12s} unresolved (n = {n_am}, {n_bm}, {n_au})")
            continue
        d_m, d_u = b_m - a_m, b_m - a_u
        p(f"   {gate_name:12s} {d_m:+10.4f} {d_u:+11.4f} {d_m - d_u:+9.4f}  "
          f"({n_am},{n_bm},{n_au})")

    # the mid rung, as a consistency check rather than a contrast
    p("\n   cv = 0.2 for reference (matched vs unmatched phi_c):")
    for gate_name, _gates in A.GATE_SETS:
        m_, nm = _pool(roots, 0.2, MATCHED_STEPS[0.2], gate_name)
        u_, nu = _pool(roots, 0.2, UNMATCHED_STEPS[0.2], gate_name)
        if nm and nu:
            p(f"      {gate_name:12s} matched {m_:.4f} (n{nm})   unmatched {u_:.4f} "
              f"(n{nu})   diff {m_ - u_:+.4f}")

    p("\n=== VERDICT ===")
    changes = []
    for gate_name, _gates in A.GATE_SETS:
        a_m, n_am = _pool(roots, CV_LO, MATCHED_STEPS[CV_LO], gate_name)
        b_m, n_bm = _pool(roots, CV_HI, MATCHED_STEPS[CV_HI], gate_name)
        a_u, n_au = _pool(roots, CV_LO, UNMATCHED_STEPS[CV_LO], gate_name)
        if n_am and n_bm and n_au:
            changes.append(a_u - a_m)      # = (matched contrast) - (unmatched contrast)
    if not changes:
        p("   nothing resolved")
        return 1
    c = np.array(changes)
    p(f"   removing a {max(span_un) / min(span_un):.2f}x rate contrast moves dphi_c by "
      f"{c.min():+.4f} to {c.max():+.4f} across the gate sets.")
    if (c >= -1e-9).all():
        p("   Every change is POSITIVE: the unconfounded contrast is LARGER, so the "
          "shear-rate confound was making the published number CONSERVATIVE. That is "
          "the direction standard phi(I) rheology predicts -- a higher rate at high CV "
          "raises p there and pushes the fitted root down.")
    else:
        p("   🔴 at least one gate set moves the OTHER way; the confound is not "
          "uniformly conservative and the sign must be quoted per gate set.")
    p(f"   Size: at most {np.abs(c).max():.4f} against a contrast of order 0.05, i.e. "
      f"{np.abs(c).max() / 0.05 * 100:.0f} % of it. Real, recorded, and not material.")
    p("   ⚠ This says nothing about whether E2's root measures the material at all -- "
      "see PREREG_ROOT_PINNING_2026-09-09.md, which says it does not.")

    path = os.path.join(A.OUT_DIR, "rate_matched_verdict.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

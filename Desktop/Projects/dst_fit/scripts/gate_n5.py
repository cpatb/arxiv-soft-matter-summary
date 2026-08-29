"""Gate N5 / G6: does one shared phi_m describe the whole corn-starch series?

G3 failed: six independent Wyart-Cates fits returned six different phi_m, from
0.492 to 0.580, where the model says it should be one number. Four candidate
explanations were written into the spec. This tests the cheapest one.

The single-curve fits pinned ``phi_0`` at 0.64, random close packing for
frictionless monodisperse spheres. Corn-starch granules are none of those
things, and because ``phi_J = (1-f) phi_0 + f phi_m`` a wrong ``phi_0`` is
absorbed straight into ``phi_m``. So: free ``phi_0``, force a single shared
``phi_m`` and ``sigma*`` across all six curves, give each curve only its own
``phi``, and let AIC decide whether that is a better description than six
independent fits.

Pass conditions, fixed before the run:

  G6a  the joint fit converges with finite standard errors on the shared parameters
  G6b  AIC(joint) < AIC(independent). The joint model has 9 parameters against
       18, so this asks whether sharing survives after paying the parameter
       count. AIC rather than a residual ratio, so the bar is not a taste
       judgement
  G6c  the fitted phi_0 is physically admissible, 0.55 <= phi_0 <= 0.71, between
       roughly random loose packing and a little short of FCC. A fit that needs
       phi_0 outside that is absorbing something else and the answer is no
  G6d  the external check survives: per-curve phi from the joint fit still tracks
       nominal phi with r2 >= 0.99. The independent fits managed 0.99925, and a
       shared model that buys its consistency by destroying that agreement has
       not explained anything

G6b is the one that decides the product. If sharing wins, joint fitting is the
package's headline feature. If it loses, corn starch is not a Wyart-Cates
material at these concentrations and the package ships against silica with this
series as a worked example of a model being refused.

Run:  python scripts/gate_n5.py
Exit code is 0 only if every gate passes.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dst_fit import fit_model, read_series, thickening_branch  # noqa: E402
from dst_fit.joint import aic_from_residuals, fit_series  # noqa: E402

DATA_DIR = Path(r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM")

PHI0_MIN, PHI0_MAX = 0.55, 0.71  # G6c
R2_MIN = 0.99  # G6d


def main() -> int:
    failures: list[str] = []
    sweeps = read_series(sorted(DATA_DIR.glob("*_e1.xls")))
    nominal = np.array([s.phi for s in sweeps])

    curves = []
    indep_res = []
    for s in sweeps:
        br = thickening_branch(s.viscosity)
        curves.append((s.stress[br], s.viscosity[br]))
        indep_res.append(fit_model("wyart_cates", s.stress, s.viscosity,
                                   branch=br, n_starts=12).residuals)

    indep_residuals = np.concatenate(indep_res)
    indep_k = 3 * len(sweeps)
    indep_aic = aic_from_residuals(indep_residuals, indep_k)
    indep_rmse = float(np.sqrt(np.mean(indep_residuals**2)))
    print(f"independent fits: k={indep_k}  n={indep_residuals.size}  "
          f"rmse={indep_rmse:.4f}  AIC={indep_aic:.1f}\n")

    for shared_ss in (True, False):
        label = "shared sigma*" if shared_ss else "per-curve sigma*"
        try:
            jf = fit_series(curves, shared_sigma_star=shared_ss, n_starts=16)
        except RuntimeError as exc:
            print(f"  {label}: FAILED to converge ({exc})")
            if shared_ss:
                failures.append(f"G6a: joint fit did not converge: {exc}")
            continue

        print(f"  {label}")
        for line in jf.summary().splitlines():
            print("    " + line)
        d_aic = jf.aic - indep_aic
        print(f"    AIC vs independent: {d_aic:+.1f} "
              f"({'joint wins' if d_aic < 0 else 'independent wins'})")

        slope, intercept = np.polyfit(nominal, jf.phi, 1)
        r2 = float(np.corrcoef(nominal, jf.phi)[0, 1] ** 2)
        print(f"    phi vs nominal: slope {slope:.4f}, intercept {intercept:+.4f}, "
              f"r2 {r2:.5f}\n")

        if not shared_ss:
            continue  # gates apply to the fully shared model

        if not (np.isfinite(jf.stderr["phi_0"]) and np.isfinite(jf.stderr["phi_m"])):
            failures.append("G6a: shared-parameter standard errors are not finite")
        if d_aic >= 0:
            failures.append(f"G6b: joint AIC is worse by {d_aic:+.1f}")
        if not (PHI0_MIN <= jf.phi_0 <= PHI0_MAX):
            failures.append(
                f"G6c: fitted phi_0 = {jf.phi_0:.4f}, outside [{PHI0_MIN}, {PHI0_MAX}]"
            )
        if r2 < R2_MIN:
            failures.append(f"G6d: phi vs nominal r2 = {r2:.5f}, below {R2_MIN}")

    if failures:
        print(f"GATE N5 FAIL ({len(failures)} problem(s)):")
        for f in failures:
            print("  " + f)
        return 1
    print("GATE N5 PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""phase2_phiJ.py — measure phi_J(CV, mu_p) and settle mechanism M2.

phi_J is NEVER read off a single run. For each (cv, mu_p, seed) a LADDER of reduced
pressures p_hat is compressed, and phi(p_hat) is extrapolated to p_hat -> 0.

The extrapolation carries its own anchor. Yade's `Law2_ScGeom_FrictPhys_CundallStrack`
is a LINEAR (harmonic) contact law, and for harmonic contacts the pressure near jamming
scales as

    p ~ (phi - phi_J)^1

so a free-exponent fit  phi = phi_J + c * p_hat^(1/alpha)  must return alpha ~ 1. That
exponent is a prediction of the contact law, not a fitted convenience, and it is
reported beside every phi_J: an alpha far from 1 means the ladder is not in the
critical regime and the phi_J from it is not trustworthy.

Results are appended as one JSON object per line to out/phase2_*.jsonl, and the script
is RESUMABLE -- a cell already present in the file is skipped. That matters because the
sweep is hours long and must survive an interruption.

Usage
-----
    ~/rl.sh phase2_phiJ.py                 # full sweep
    ~/rl.sh phase2_phiJ.py --quick         # small, for shakeout
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

import packing
from packing import build_and_compress

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# --- the ladder ------------------------------------------------------------------
# p_hat spans 1.5 decades. The 3e-3 rung was dropped after a timing measurement:
# the sweep ran at 0.43 cells/min across 8 workers, so the original 270-cell grid
# was a 10.5 h job. 3e-3 is the cheapest cell but `analyse_phase2.P_HAT_FIT_MAX`
# excludes it from the extrapolation anyway (it is out of the critical regime), so
# dropping it costs the fit nothing. All four remaining rungs are INSIDE the fit
# window, which the free-exponent check needs (it requires >= 4 points).
P_HAT_LADDER = (1e-3, 3e-4, 1e-4, 3e-5)

# CV ladder. Capped at 0.5: at CV = 0.8 the (0.2, 5.0) truncation needed to keep the
# largest grain inside the periodic cell distorts the distribution badly (realised
# CV 0.72, d32 21% below analytic), so a nominal 0.8 would not be a 0.8. Everything
# downstream is reported against the REALISED cv, never the nominal one.
CV_LADDER = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)

# mu_p = 1.0 dropped from the production grid for the same timing reason. The A4
# anchor (phi_J falls with friction) is still testable on 0.0 vs 0.3, and a
# supplementary mu = 1.0 run can be added with --tag if time allows.
MU_LADDER = (0.0, 0.3)
SEEDS = (1, 2, 3)

# N is bounded below by the periodic collider: at CV = 0.5 with truncate=(0.2, 5.0)
# the largest grain needs N >= 1299. 1500 clears every cell in CV_LADDER.
N_PROD = 1500
N_QUICK = 600


def fit_phi_J(p_hat: np.ndarray, phi: np.ndarray) -> dict[str, float]:
    """Extrapolate phi(p_hat) -> 0 two ways, and report both.

    `phi_J_lin`   : linear fit, phi = phi_J + c*p_hat. This is what a harmonic contact
                    law predicts and is the PRIMARY estimate.
    `phi_J_free`  : free-exponent fit, phi = phi_J + c*p_hat^b. `b` is the check: for a
                    linear contact law it must come back near 1.

    Both are least squares. With <= 3 points the free fit is under-determined and is
    returned as NaN rather than a number with no support.
    """
    from scipy.optimize import curve_fit

    p = np.asarray(p_hat, float)
    f = np.asarray(phi, float)
    ok = np.isfinite(p) & np.isfinite(f) & (p > 0)
    p, f = p[ok], f[ok]
    out = {"n_points": int(p.size), "phi_J_lin": float("nan"),
           "slope_lin": float("nan"), "r2_lin": float("nan"),
           "phi_J_free": float("nan"), "b_free": float("nan"),
           "phi_J_lin_se": float("nan")}
    if p.size < 2:
        return out

    A = np.vstack([np.ones_like(p), p]).T
    coef, res, *_ = np.linalg.lstsq(A, f, rcond=None)
    out["phi_J_lin"] = float(coef[0])
    out["slope_lin"] = float(coef[1])
    pred = A @ coef
    ss_tot = float(np.sum((f - f.mean()) ** 2))
    out["r2_lin"] = float(1 - np.sum((f - pred) ** 2) / ss_tot) if ss_tot > 0 else np.nan
    if p.size > 2:
        dof = p.size - 2
        s2 = float(np.sum((f - pred) ** 2) / dof)
        cov = s2 * np.linalg.inv(A.T @ A)
        out["phi_J_lin_se"] = float(np.sqrt(cov[0, 0]))

    if p.size >= 4:
        try:
            def model(x, phij, c, b):
                return phij + c * np.power(x, b)
            popt, _ = curve_fit(model, p, f,
                                p0=[out["phi_J_lin"], max(out["slope_lin"], 1e-6), 1.0],
                                maxfev=20000)
            out["phi_J_free"], out["b_free"] = float(popt[0]), float(popt[2])
        except Exception as exc:            # noqa: BLE001 - report, do not hide
            out["fit_error"] = repr(exc)
    return out


def fs_safe(family: str) -> str:
    """Make a family string usable as part of a FILENAME, without touching the key.

    Parameterised families are encoded in the family STRING (`binary:r3.0:x0.30`,
    `multiplet:f1.8:n0.70/0.20/0.10`), which is right for the record and wrong for a
    path: `:` is illegal on Windows and `/` would silently create a DIRECTORY on Linux,
    so a sharded run would scatter its output instead of failing.

    ⚠ THE CELL KEY IS NOT SANITISED and must not be -- it is what `load_done` resumes
    on and what every record on disk carries. Only the filename is.

    Every family that has ever been run through this script (`lognormal`, `uniform`,
    `bidisperse`) contains none of these characters, so this is byte-identical for all
    existing output names; `test_multiplet.py::test_fs_safe_does_not_move_existing_names`
    pins that.

    >>> fs_safe("multiplet:f1.8:n0.70/0.20/0.10")
    'multiplet-f1.8-n0.70_0.20_0.10'
    """
    return family.replace(":", "-").replace("/", "_")


def cell_key(cv: float, mu: float, seed: int, p_hat: float, n: int,
             family: str = "lognormal") -> str:
    # The family suffix is omitted for lognormal so that keys written before
    # the uniform family existed still match and stay resumable.
    base = f"cv{cv:.3f}_mu{mu:.3f}_s{seed}_p{p_hat:.3e}_N{n}"
    key = base if family == "lognormal" else f"{base}_{family}"
    # The iteration cap is part of the protocol, so it is part of the key -- but ONLY
    # when it is non-default, so every key written before 2026-09-09 still matches and
    # stays resumable. Without this a re-run at a higher cap would collide with the
    # capped record it exists to replace, and analyse_phase2's dedup keeps whichever
    # shard it read last.
    cap = packing.max_iter()
    return key if cap == packing.MAX_ITER else f"{key}_it{cap}"


def load_done(path: str) -> set[str]:
    done = set()
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["key"])
                except Exception:
                    continue
    return done


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--tag", default="")
    ap.add_argument("--shard", type=int, default=0,
                    help="index of this worker, 0-based")
    ap.add_argument("--nshards", type=int, default=1,
                    help="total number of workers; cells are dealt round-robin")
    ap.add_argument("--family", default="lognormal",
                    choices=["lognormal", "uniform", "bidisperse"],
                    help="PSD family; uniform has skewness 0 at every CV")
    ap.add_argument("--n", type=int, default=0, help="override N")
    ap.add_argument("--mus", default="",
                    help="comma-separated mu_p values, overriding MU_LADDER")
    ap.add_argument("--cvs", default="",
                    help="comma-separated CV values, overriding CV_LADDER")
    ap.add_argument("--seeds", default="",
                    help="comma-separated seeds, overriding SEEDS")
    ap.add_argument("--phats", default="",
                    help="comma-separated p_hat rungs, overriding P_HAT_LADDER. A "
                         "SINGLE rung is not a phi_J measurement -- phi_J is always "
                         "extrapolated from the ladder -- so any campaign that passes "
                         "this must say in its prereg that it is reading an observable "
                         "AT a pressure, not a jamming point.")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    n = args.n or (N_QUICK if args.quick else N_PROD)
    cvs = (0.0, 0.3) if args.quick else CV_LADDER
    mus = (0.0, 0.3) if args.quick else MU_LADDER
    if args.mus:
        mus = tuple(float(x) for x in args.mus.split(','))
    if args.cvs:
        cvs = tuple(float(x) for x in args.cvs.split(','))
    # --seeds AFTER the default, not before: until 2026-09-18 these two lines were the
    # other way round, so the flag was accepted and silently overwritten. Every campaign
    # that passed it passed `1,2,3`, which is SEEDS, so no record on disk is affected --
    # but a 6-seed run would have quietly been a 3-seed one.
    seeds = (1,) if args.quick else SEEDS
    if args.seeds:
        seeds = tuple(int(x) for x in args.seeds.split(','))
    ladder = (1e-3, 3e-4, 1e-4) if args.quick else P_HAT_LADDER
    if args.phats:
        ladder = tuple(float(x) for x in args.phats.split(','))

    suffix = ("_quick" if args.quick else "") + (f"_{args.tag}" if args.tag else "")
    if args.family != "lognormal":
        suffix += f"_{fs_safe(args.family)}"
    if args.nshards > 1:
        suffix += f"_shard{args.shard:02d}of{args.nshards:02d}"
    path = os.path.join(OUT_DIR, f"phase2{suffix}.jsonl")
    done = load_done(path)
    total = len(cvs) * len(mus) * len(seeds) * len(ladder)
    mine = (total + args.nshards - 1 - args.shard) // args.nshards
    print(f"phase2: {total} cells total, {mine} on shard {args.shard}/"
          f"{args.nshards}, N={n}, {len(done)} already done -> {path}", flush=True)

    t0 = time.time()
    idx = 0
    for cv in cvs:
        for mu in mus:
            for seed in seeds:
                for p_hat in ladder:
                    idx += 1
                    # Round-robin deal so every worker gets a mix of cheap and
                    # expensive cells; blocking by cv would leave one worker with
                    # all the CV=0.5 cells and the rest idle.
                    if (idx - 1) % args.nshards != args.shard:
                        continue
                    key = cell_key(cv, mu, seed, p_hat, n, args.family)
                    if key in done:
                        continue
                    try:
                        st = build_and_compress(cv=cv, mu_p=mu, p_hat=p_hat, n=n,
                                                seed=seed, verbose=True,
                                                family=args.family)
                        rec = json.loads(st.to_json())
                    except Exception as exc:      # noqa: BLE001
                        print(f"CELL FAILED {key}: {exc!r}", flush=True)
                        rec = {"failed": True, "error": repr(exc), "cv": cv,
                               "mu_p": mu, "seed": seed, "p_hat_target": p_hat,
                               "n_requested": n, "family": args.family}
                    rec["key"] = key
                    with open(path, "a", encoding="utf-8") as fh:
                        fh.write(json.dumps(rec, sort_keys=True) + "\n")
                    el = time.time() - t0
                    print(f"  [{idx}/{total}] elapsed {el/60:.1f} min", flush=True)
    print(f"phase2 done in {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
fit_closure.py -- post-process one impinging-jet CFD case:
  * read the sampled plate wallShearStress (kinematic, m^2/s^2) and yPlus,
  * convert to tau_w(r) in Pa (multiply by rho),
  * fit the Phares (2000) turbulent closure over the developed wall-jet band:
        (tau_w / (rho V0^2)) * Re_d^{1/2} = C_P (r/d)^{-n_t}
    -> log-log linear fit gives n_t (slope) and C_P (intercept),
  * report y+ range on the plate and the tau_w peak,
  * write tau_w_r.csv  (r_mm, r_over_d, tau_w_Pa, phares_LHS).

Usage: python3 fit_closure.py <case_dir> [fit_lo_rOverd] [fit_hi_rOverd]
Prints a one-line JSON summary consumed by run_sweep.py.
"""
import sys, os, glob, json
import numpy as np

RHO = 1000.0
D = 1.0e-3
NU = 8.9e-7


def read_info(case):
    info = {}
    for line in open(os.path.join(case, "case_info.txt")):
        kv = line.strip().split("=")
        if len(kv) == 2:
            try: info[kv[0]] = float(kv[1])
            except ValueError: info[kv[0]] = kv[1]
    return info


def latest_raw(case, field):
    # postProcessing/sampleWSS/<time>/<field>_plateSurf.raw
    pats = glob.glob(os.path.join(case, "postProcessing", "sampleWSS", "*", f"{field}_*.raw"))
    if not pats:
        return None
    # pick the highest time directory
    def t_of(p):
        try: return float(os.path.basename(os.path.dirname(p)))
        except ValueError: return -1
    return sorted(pats, key=t_of)[-1]


def load_raw(path, ncomp):
    """raw: columns x y z [field comps]. Returns (xyz Nx3, field Nx ncomp)."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    xyz = data[:, :3]
    fld = data[:, 3:3 + ncomp]
    return xyz, fld


def main():
    case = sys.argv[1].rstrip("/")
    fit_lo = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
    fit_hi = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    info = read_info(case)
    V0, Re = info["V0"], info["Re_d"]

    wss_path = latest_raw(case, "wallShearStress")
    if wss_path is None:
        print(json.dumps({"case": case, "error": "no wallShearStress sample"}))
        return
    xyz, wss = load_raw(wss_path, 3)
    r = np.hypot(xyz[:, 1], xyz[:, 2])            # radial (y,z); x is axial = H
    tau_kin = np.linalg.norm(wss, axis=1)          # |tau|/rho  [m^2/s^2]
    tau_Pa = tau_kin * RHO
    order = np.argsort(r)
    r, tau_Pa = r[order], tau_Pa[order]
    r_over_d = r / D

    # y+ (optional)
    yp_path = latest_raw(case, "yPlus")
    yplus = None
    if yp_path is not None:
        _, yp = load_raw(yp_path, 1)
        yplus = yp.ravel()

    # Phares LHS
    LHS = (tau_Pa / (RHO * V0**2)) * np.sqrt(Re)

    # fit band: developed wall jet, exclude stagnation & outer boundary
    band = (r_over_d >= fit_lo) & (r_over_d <= fit_hi) & (tau_Pa > 0)
    n_fit = int(band.sum())
    if n_fit >= 4:
        lx = np.log(r_over_d[band]); ly = np.log(LHS[band])
        slope, intercept = np.polyfit(lx, ly, 1)
        n_t = -slope
        C_P = float(np.exp(intercept))
        resid = ly - (slope * lx + intercept)
        ss_res = float(np.sum(resid**2)); ss_tot = float(np.sum((ly - ly.mean())**2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    else:
        n_t = C_P = r2 = float("nan")

    # peak
    ipk = int(np.argmax(tau_Pa))
    summary = dict(
        case=os.path.basename(case), Re_d=Re, H_over_d=info["H_over_d"], V0=V0,
        refine=info.get("refine", 1.0),
        n_t=float(n_t), C_P=float(C_P), r2=float(r2), n_fit=n_fit,
        fit_band=[fit_lo, fit_hi],
        tau_peak_Pa=float(tau_Pa[ipk]), r_peak_mm=float(r[ipk] * 1e3),
        tau_at_15mm_Pa=float(np.interp(0.015, r, tau_Pa)),
        yplus_min=float(np.nanmin(yplus)) if yplus is not None else None,
        yplus_max=float(np.nanmax(yplus)) if yplus is not None else None,
        yplus_median=float(np.nanmedian(yplus)) if yplus is not None else None,
    )

    # write csv
    out = np.column_stack([r * 1e3, r_over_d, tau_Pa, LHS])
    np.savetxt(os.path.join(case, "tau_w_r.csv"), out,
               header="r_mm,r_over_d,tau_w_Pa,phares_LHS", delimiter=",", comments="")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()

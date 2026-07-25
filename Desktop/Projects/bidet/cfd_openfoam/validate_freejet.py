#!/usr/bin/env python3
"""
validate_freejet.py <case_dir> -- validation check (a): the free-jet centerline
velocity decay.  For a round turbulent jet the centerline speed obeys
    Vc(x)/V0 = B * d / (x - x0)   with B ~ 5-6 (Hussein-Capp-George 1994),
in the developed region before the plate influence.  We sample U along a near-axis
line and fit B over x/d in [8, H/d - 4], excluding the nozzle near-field and the
stagnation region just above the plate.

Runs `postProcess -func <set>` in the case, reads the sampled line, prints B.
Only meaningful for tall standoffs (H/d >= 20).
"""
import sys, os, glob, subprocess
import numpy as np

D = 1.0e-3


def read_info(case):
    info = {}
    for line in open(os.path.join(case, "case_info.txt")):
        kv = line.strip().split("=")
        if len(kv) == 2:
            try: info[kv[0]] = float(kv[1])
            except ValueError: pass
    return info


def main():
    case = sys.argv[1].rstrip("/")
    info = read_info(case)
    V0, H, Hd = info["V0"], info["H"], info["H_over_d"]

    # write a sampledSets dict: line along axis (r ~ 0.03 mm) from x=0 to x=H
    roff = 0.03e-3
    setdict = f"""FoamFile {{ version 2.0; format ascii; class dictionary; object sampleAxis; }}
type sets;
libs (sampling);
interpolationScheme cellPoint;
setFormat raw;
fields (U);
sets
(
    axisLine
    {{
        type    uniform;
        axis    x;
        start   (0 {roff} 0);
        end     ({H} {roff} 0);
        nPoints 300;
    }}
);
"""
    with open(os.path.join(case, "system", "sampleAxis"), "w") as f:
        f.write(setdict)

    env = dict(os.environ)
    r = subprocess.run(["postProcess", "-func", "sampleAxis", "-latestTime"],
                       cwd=case, env=env, capture_output=True, text=True)
    # find output
    pats = glob.glob(os.path.join(case, "postProcessing", "sampleAxis", "*", "axisLine_U.*"))
    if not pats:
        print(f"{os.path.basename(case)}: no axis sample produced\n{r.stderr[-400:]}")
        return
    data = np.loadtxt(sorted(pats)[-1])
    x = data[:, 0]; Ux = data[:, 1]           # raw sets: x_coord then U components
    Vc = np.abs(Ux)

    xd = x / D
    fit = (xd >= 8) & (xd <= (Hd - 4)) & (Vc > 0.05 * V0)
    if fit.sum() < 5:
        print(f"{os.path.basename(case)} (H/d={Hd:.0f}): too few free-jet points "
              f"(H/d too small for a developed core); skipping B fit")
        return
    # Vc/V0 = B d/x  ->  V0/Vc = (1/B) (x/d)  ->  slope = 1/B
    slope = np.polyfit(xd[fit], (V0 / Vc[fit]), 1)[0]
    B = 1.0 / slope
    # also report the raw ratio at a couple stations
    print(f"{os.path.basename(case)} (H/d={Hd:.0f}): free-jet decay constant B = {B:.2f} "
          f"(target ~5-6; Hussein-Capp-George)  [npts={fit.sum()}]")


if __name__ == "__main__":
    main()

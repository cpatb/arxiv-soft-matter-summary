"""check_control_matching.py -- is the sigma_u control matched to the cases it divides?

FRISBEE_GENERALIZATION.md section 1 divides all three production cases by a single
control (`abl_control.npz`, sampled from `dir_090_empty`) and reports that buildings
RAISED sigma_u by about half again at 210 and 270 deg. CLAUDE_CODE_FRISBEE_HARDENING.md
section C2 says that comparison is confounded: the control and `dir_090` are fed
z0_inlet = 0.001 while `dir_210` and `dir_270` are fed z0_inlet = 0.5, so the two
directions carrying the claim meet 3x the inlet fluctuation before reaching a building.

This script re-derives that from the case metadata rather than trusting either document
(the standing rule: reproduce the claim's own numbers before extending or repeating it).

    python check_control_matching.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE / "SeanEarl_DiscFlight" / "P3_CFD" / "cases"

CONTROL = "dir_090_empty"
PRODUCTION = ["dir_090", "dir_210", "dir_270"]

# sigma_u ratios (production / control) as reported in the memo, via
# vertical_generated.tex \VertSigRatioH*.
MEMO_SIGMA_RATIO = {"dir_090": 0.91, "dir_210": 1.51, "dir_270": 1.57}

FIELDS = ["z0_inlet", "z0_ground", "u_star", "k_inlet_ref"]


def load(case: str) -> dict:
    p = CASES / case / "case_info.json"
    if not p.is_file():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def get(d: dict, key: str):
    """case_info.json may nest the inflow parameters; search one level down too."""
    if key in d:
        return d[key]
    for v in d.values():
        if isinstance(v, dict) and key in v:
            return v[key]
    return None


def main() -> int:
    ctl = load(CONTROL)
    if not ctl:
        print("control metadata not found:", CASES / CONTROL / "case_info.json")
        return 1

    print("=" * 78)
    print("Is the sigma_u control matched to the cases it is used to divide?")
    print("=" * 78)
    print("%-16s %10s %11s %9s %9s   %s"
          % ("case", "z0_inlet", "z0_ground", "u_star", "k_inlet", "matched?"))

    ctl_vals = {f: get(ctl, f) for f in FIELDS}
    print("%-16s %10s %11s %9s %9s   %s"
          % (CONTROL + " *", ctl_vals["z0_inlet"], ctl_vals["z0_ground"],
             ctl_vals["u_star"], ctl_vals["k_inlet_ref"], "(control)"))

    rows = []
    for case in PRODUCTION:
        d = load(case)
        if not d:
            print("%-16s  METADATA MISSING -- cannot grade" % case)
            return 1
        vals = {f: get(d, f) for f in FIELDS}
        matched = all(vals[f] == ctl_vals[f] for f in ("z0_inlet", "k_inlet_ref"))
        rows.append((case, vals, matched))
        print("%-16s %10s %11s %9s %9s   %s"
              % (case, vals["z0_inlet"], vals["z0_ground"], vals["u_star"],
                 vals["k_inlet_ref"], "YES" if matched else "*** NO ***"))

    print()
    print("CORRECTED sigma_u RATIOS")
    print("------------------------")
    print("An inlet carrying more turbulent kinetic energy delivers more sigma_u before")
    print("the flow meets any building. sigma_u ~ sqrt(k), so the approach-flow ratio is")
    print("sqrt(k_case / k_control); dividing it out gives the buildings' own effect.")
    print()
    print("%-16s %14s %12s %14s   %s"
          % ("case", "memo ratio", "approach", "vs matched", "reading"))
    any_confounded = False
    for case, vals, matched in rows:
        memo = MEMO_SIGMA_RATIO[case]
        kc, k0 = vals["k_inlet_ref"], ctl_vals["k_inlet_ref"]
        if kc is None or k0 is None or k0 <= 0:
            print("%-16s  k_inlet unavailable -- cannot correct" % case)
            return 1
        approach = math.sqrt(kc / k0)
        corrected = memo / approach
        if not matched:
            any_confounded = True
        print("%-16s %14.2f %12.2f %14.2f   %s"
              % (case, memo, approach, corrected,
                 "roughness-clean" if matched
                 else ("buildings REDUCE" if corrected < 1.0 else "buildings raise")))

    print()
    print("VERDICT")
    print("-------")
    if not any_confounded:
        print("=> All production cases share the control's inflow. No roughness confound.")
        return 0

    clean = [(c, MEMO_SIGMA_RATIO[c]) for c, _, m in rows if m]
    print("=> CONFOUNDED. The only roughness-clean comparison is %s"
          % ", ".join("%s (ratio %.2f)" % cv for cv in clean))
    print("   -- and that is the heading where the memo's ratio is already BELOW 1, i.e.")
    print("   where the 'buildings raise sigma_u' claim does not hold anyway.")
    print()
    print("   So the claim rests entirely on the two headings whose inflow is NOT matched")
    print("   to the control, and against matched approach flow every corrected ratio")
    print("   below 1 means the buildings REDUCED sigma_u. The headline cannot be stated")
    print("   from these cases; it needs an empty-domain control run at z0_inlet = 0.5.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

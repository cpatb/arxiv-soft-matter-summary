"""grade_jamwin2.py -- grades PREREG_JAMWIN_2026-09-22.md Amendment 2 (the contingency).

Written before any jamwin2 cell ran. Input out/jamwin2/<label>/ from
`scripts/pull_jamwin.sh jamwin2`. Also grades the jamwin arm-B cells with the same binary
outcome, printed separately and never pooled.
"""

from __future__ import annotations

import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIGS = ("10", "31.6", "100", "316", "1000", "3162", "10000")
DELTAS = ("0.7", "1.3169", "2.6")
PHI_B, JF_B = 0.527, 0.75
T_END_MIN = 5.9
JAM_LEVEL = 0.5


def outcome(root: str, lab: str) -> tuple[str, dict]:
    d = os.path.join(root, lab)
    st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
    ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
    info = {"rc": st["rc"], "t_end": st["last_t_over_tau"]}
    m = re.search(r"phi=([\d.]+) .*sigma\*=([\d.e+]+) Pa delta=([\d.]+)", ban)
    jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
    info["prov"] = bool(m and jb) and math.isclose(float(m.group(1)), PHI_B) \
        and math.isclose(float(jb.group(1)), JF_B)
    dat = os.path.join(d, "spread.dat")
    if not os.path.exists(dat):
        return "UNDETERMINED", info
    a = np.loadtxt(dat, comments="#", ndmin=2)
    jam = a[:, 7]
    info["jmax"] = float(jam.max()) if len(jam) else math.nan
    hit = np.flatnonzero(jam >= JAM_LEVEL)
    if hit.size:
        info["t_jam"] = float(a[hit[0], 1])
        return "JAMMED", info
    if st["rc"] == 0 and a[-1, 1] >= T_END_MIN:
        return "NOT JAMMED", info
    return "UNDETERMINED", info


def main() -> int:
    out = []

    def say(s=""):
        print(s)
        out.append(s)

    root = os.path.join(HERE, "out", "jamwin2")
    say("grade_jamwin2 -- PREREG_JAMWIN Amendment 2 (arm B, x = 0.85, f_c = 0.75)")
    say(f"outcome: JAMMED iff jam_frac(fhat > f_c) >= {JAM_LEVEL} at a written row; "
        f"NOT JAMMED iff survives to t/tau >= {T_END_MIN}; else UNDETERMINED")
    table = {}
    for d in DELTAS:
        for s in SIGS:
            lab = f"B_s{s}_d{d}"
            if os.path.exists(os.path.join(root, lab, "RUN_STATUS.json")):
                table[(d, s)] = outcome(root, lab)
    say(f"cells landed: {len(table)} of 21;  provenance ok: "
        f"{sum(v[1]['prov'] for v in table.values())}/{len(table)}")
    say("")
    say("   delta \\ sigma*  " + "  ".join(f"{s:>7s}" for s in SIGS))
    sc = {}
    for d in DELTAS:
        row = [table.get((d, s), ("-", {}))[0] for s in SIGS]
        say(f"   {d:>7s}         " + "  ".join(
            {"JAMMED": "    JAM", "NOT JAMMED": "   free", "UNDETERMINED": "      ?",
             "-": "      -"}[r] for r in row))
        # sigma*_c: largest JAMMED rung, only if outcomes are monotone (JAM below, free above)
        known = [(float(s), r) for s, r in zip(SIGS, row) if r in ("JAMMED", "NOT JAMMED")]
        jams = [v for v, r in known if r == "JAMMED"]
        frees = [v for v, r in known if r == "NOT JAMMED"]
        monotone = not jams or not frees or max(jams) < min(frees)
        complete = all(r != "-" for r in row)
        sc[d] = max(jams) if (jams and monotone and complete) else None
        say(f"            sigma*_c = {sc[d]}   (monotone={monotone}, complete={complete})")
    say("")
    vals = [sc[d] for d in DELTAS]
    if any(v is None for v in vals):
        verdict = "UNDETERMINED (a delta has no readable boundary)"
    elif vals[0] == vals[1] == vals[2]:
        verdict = "UNDETERMINED (tie at this ladder spacing)"
    elif vals[0] <= vals[1] <= vals[2] and vals[2] > vals[0]:
        verdict = "(S) stress-thresholded: sigma*_c RISES with delta"
    elif vals[0] >= vals[1] >= vals[2] and vals[0] > vals[2]:
        verdict = "(R) rate-driven bottleneck: sigma*_c FALLS with delta"
    else:
        verdict = "NEITHER (non-monotone in delta)"
    say(f"VERDICT: {verdict}")
    with open(os.path.join(HERE, "out", "grade_jamwin2.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

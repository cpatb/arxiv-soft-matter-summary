"""grade_so_triplet.py -- grades PREREG_so_triplet_2026-09-27.md exactly as registered.

Gates S1-S6 plus the Amendment-2 health rule (H1: every initial residual except Uz below
1e-5 at the last iteration; H2: max|Uz|/max|U_inplane| < 1e-10 in the final field).
The original residualControl verdict is printed beside H1/H2 for every case.

Refuses to grade a family until all three members exist and are healthy (a set-bar is
not graded on a partial set). Quantities come from gci.fit_record (fit_closure.py, one
code path for every record); orders from gci.observed_order (Celik 2008 with the
divergence guard).

Writes convergence/so_triplet_grade.md and prints the same. Exit 0 always: the verdicts
are the result.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gci  # noqa: E402

QOI = [q for q, _ in gci.QOI]
BAND = gci.ASYMPTOTIC_BAND
TARGET = 1e-5
UZ_ROUNDOFF = 1e-10
CAUCHY_TOL = 0.005  # S6
S5_BAR = 0.01
FAMILIES = {"so_r15": gci.FAMILIES["so_r15"], "fo_r15": gci.FAMILIES["fo_r15"]}
MAPPED = ("case_Re10000_H6_fo_r2p25", "case_Re10000_H6_so_r2p25", "case_Re10000_H6_fo_r1p5")


def last_residuals(case):
    pat = re.compile(r"Solving for (\w+),\s+Initial residual = ([0-9.eE+-]+)")
    it, cur, done, conv = None, {}, None, False
    with open(os.path.join(HERE, case, "log.simpleFoam"), errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                if it is not None and cur:
                    done = (it, cur)
                it, cur = int(float(line.split("=")[1])), {}
            elif "SIMPLE solution converged" in line:
                conv = True
            else:
                m = pat.search(line)
                if m and m.group(1) not in cur:
                    cur[m.group(1)] = float(m.group(2))
    if cur:
        done = (it, cur)
    return done, conv


def uz_ratio(case):
    """max|Uz| / max|U_inplane| in the final written U field."""
    d = os.path.join(HERE, case)
    times = sorted((t for t in os.listdir(d) if re.fullmatch(r"[0-9]+", t) and t != "0"
                    and os.path.exists(os.path.join(d, t, "U"))), key=int)
    if not times:
        return float("nan"), None
    txt = open(os.path.join(d, times[-1], "U")).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\(", txt)
    n = int(m.group(1))
    body = txt[m.end():]
    body = body[: body.find("\n)")]
    U = np.array(re.findall(r"\(([-0-9.eE+]+) ([-0-9.eE+]+) ([-0-9.eE+]+)\)", body), dtype=float)
    assert U.shape == (n, 3), (case, U.shape, n)
    return float(np.abs(U[:, 2]).max() / np.linalg.norm(U[:, :2], axis=1).max()), times[-1]


def health(case):
    if not os.path.exists(os.path.join(HERE, case, "log.simpleFoam")):
        return None
    (it, res), conv = last_residuals(case)
    nz = {k: v for k, v in res.items() if k != "Uz"}
    h1 = max(nz.values()) < TARGET
    r, t = uz_ratio(case)
    h2 = bool(np.isfinite(r) and r < UZ_ROUNDOFF)
    return dict(iterations=it, worst_nonUz=max(nz.values()), Uz=res.get("Uz", float("nan")),
                h1=h1, uz_ratio=r, final_time=t, h2=h2, residualControl=conv, ok=h1 and h2)


def main():
    L = ["# PREREG_so_triplet grade\n"]
    rec, hl = {}, {}
    allcases = sorted({c for fam in FAMILIES.values() for c, _ in fam} | {m + "_cont" for m in MAPPED})
    L.append("## Health (Amendment 2)\n")
    L.append("| case | iter | worst non-Uz | Uz resid | H1 | max|Uz|/max|U| | H2 | orig. residualControl |")
    L.append("|---|---|---|---|---|---|---|---|")
    for c in allcases:
        h = health(c)
        hl[c] = h
        if h is None:
            L.append(f"| {c} | not run | | | | | | |")
            continue
        L.append(f"| {c} | {h['iterations']} | {h['worst_nonUz']:.2e} | {h['Uz']:.2e} | "
                 f"{'PASS' if h['h1'] else 'FAIL'} | {h['uz_ratio']:.1e} | "
                 f"{'PASS' if h['h2'] else 'FAIL'} | {'met' if h['residualControl'] else 'NOT met'} |")
        if h["ok"]:
            rec[c] = gci.fit_record(c)

    L.append("\n## S6: Cauchy continuation of the mapped-start cases (+400 iterations)\n")
    s6 = {}
    for m in MAPPED:
        if m not in rec or m + "_cont" not in rec:
            s6[m] = None
            L.append(f"- {m}: NOT GRADED (base or continuation missing/unhealthy)")
            continue
        moves = {q: abs(rec[m + "_cont"][q] - rec[m][q]) / abs(rec[m][q]) for q in QOI}
        s6[m] = max(moves.values()) < CAUCHY_TOL
        L.append(f"- {m}: " + ", ".join(f"{q} {100 * v:.3f}%" for q, v in moves.items())
                 + f" -> {'PASS' if s6[m] else 'FAIL'} (bar {100 * CAUCHY_TOL:.1f}%)")

    L.append("\n## Triplets\n")
    verdict = {}
    for fam, members in FAMILIES.items():
        names = [c for c, _ in members]
        refs = [r for _, r in members]
        missing = [c for c in names if c not in rec]
        s6_bad = [c for c in names if c in s6 and s6[c] is not True]
        if missing or s6_bad:
            verdict[fam] = None
            L.append(f"**{fam}: INCOMPLETE** -- missing/unhealthy {missing}, S6 not passed {s6_bad}\n")
            continue
        r21, r32 = refs[0] / refs[1], refs[1] / refs[2]
        L.append(f"**{fam}** ({', '.join(names)}), r21 = {r21:g}, r32 = {r32:g}\n")
        L.append("| quantity | finest | middle | coarsest | R = eps21/eps32 | p | class |")
        L.append("|---|---|---|---|---|---|---|")
        verdict[fam] = {}
        for q in QOI:
            v1, v2, v3 = (rec[c][q] for c in names)
            eps21, eps32 = v2 - v1, v3 - v2
            R = eps21 / eps32 if eps32 != 0 else float("nan")
            p, s, ok, why = gci.observed_order(v1, v2, v3, r21, r32)
            cls = ("monotone converging" if 0 < R < 1 else
                   "oscillatory" if -1 < R < 0 else "DIVERGING" if abs(R) >= 1 else "degenerate")
            verdict[fam][q] = (R, p, cls)
            L.append(f"| {q} | {v1:.5g} | {v2:.5g} | {v3:.5g} | {R:+.3f} | "
                     f"{p:.3f} | {cls} |" if ok else
                     f"| {q} | {v1:.5g} | {v2:.5g} | {v3:.5g} | {R:+.3f} | -- | {cls} |")
        L.append("")

    L.append("## Gates\n")
    so, fo = verdict.get("so_r15"), verdict.get("fo_r15")
    if so:
        R, p, cls = so["C_P"]
        s1 = 0 < R < 1 and BAND[0] <= p <= BAND[1]
        L.append(f"- **S1** SO C_P: R = {R:+.3f}, p = {p:.3f} -> {'PASS' if s1 else 'FAIL'}")
        for q in ("n_t", "tau_at_15mm_Pa"):
            R2, p2, c2 = so[q]
            L.append(f"- **S2** SO {q}: R = {R2:+.3f} ({c2}), p = {p2:.3f} -> "
                     f"{'PASS' if 0 < R2 < 1 else 'FAIL'}")
    else:
        L.append("- S1, S2: INCOMPLETE")
    if fo:
        R, p, cls = fo["C_P"]
        fo_conv = 0 < R < 1
        L.append(f"- **S3** FO C_P on the same grids: R = {R:+.3f} ({cls}) -> "
                 + ("scheme contrast INFORMATIVE (FO not monotone-converging)" if not fo_conv
                    else "scheme contrast NOT INFORMATIVE (FO converges too)"))
    else:
        L.append("- S3: INCOMPLETE")
    for c in ("case_Re10000_H6_fo_r1p5",):
        h = hl.get(c)
        L.append("- **S4** FO refine 1.5 health: " + ("not run" if h is None else
                 f"H1 {'PASS' if h['h1'] else 'FAIL'} (worst non-Uz {h['worst_nonUz']:.2e}); "
                 f"original residualControl {'met' if h['residualControl'] else 'NOT met'}"))
    if "case_Re10000_H6" in rec and "case_Re10000_H6_so_r1" in rec:
        d = abs(rec["case_Re10000_H6_so_r1"]["C_P"] - rec["case_Re10000_H6"]["C_P"]) / rec["case_Re10000_H6"]["C_P"]
        L.append(f"- **S5** |dC_P|/C_P at refine 1, SO vs FO = {100 * d:.2f}% -> "
                 f"{'PASS' if d > S5_BAR else 'FAIL'} (bar {100 * S5_BAR:.0f}%)")
    out = "\n".join(L) + "\n"
    os.makedirs(os.path.join(HERE, "convergence"), exist_ok=True)
    with open(os.path.join(HERE, "convergence", "so_triplet_grade.md"), "w", encoding="utf-8") as fh:
        fh.write(out)
    print(out)


if __name__ == "__main__":
    main()

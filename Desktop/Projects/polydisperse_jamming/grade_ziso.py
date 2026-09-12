"""grade_ziso.py — grade PREREG_ZISO_2026-09-12.md, gates G1-G7.

Written against the bars in that file, which were fixed before any z_iso was
re-measured. Nothing here computes a bar; every threshold is quoted.

    python grade_ziso.py          -> out/ziso_audit.txt

No simulation. Reads out/phase2_*.jsonl (isotropic ladders) and
out/phase8_commonlad*.jsonl (the common shear ladder), plus the COMMITTED artifact
out/z_iso_measured.txt, which G2 reproduces against and which nothing here overwrites.
"""

from __future__ import annotations

import os
import re

import numpy as np

import analyse_phase8 as A
import measure_z_iso as M
from grade_common_ladder import cells_by_group, contrast

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# --- bars, all QUOTED from PREREG_ZISO_2026-09-12.md §5 ----------------------------
G2_BAR = 0.005          # contacts; --cap default must reproduce the committed artifact
G3_SPAN_BAR = 0.10      # contacts; fit-rule systematic on z_iso(mu_p = 0.3)
G3_MU3_BAR = 0.15       # contacts; |z_iso(mu_p = 3) - 4| for a rule to be admissible
G4_DRIFT_BAR = 0.10     # contacts; |lognormal drift| must survive the repair
G5_BAR = 0.010          # in phi; how far delta phi_c may move before requalifying
G6_TOL = 1e-9           # the refactor control on the default path, AS PRE-REGISTERED
G6A_BAR = 5e-5          # its replacement: the artifact prints d to four decimals
PUBLISHED_E4 = +0.0554  # common ladder, ungated, out/common_ladder_grade.txt
PUBLISHED_E4M = +0.0497
Z_ISO_THEORY_FRICTIONLESS = 6.0     # 2d, exact
Z_ISO_THEORY_INFINITE_MU = 4.0      # d + 1, the mu_p -> infinity limit
RULES_TESTED = ("broken", "r2", "r3")
COMMON_PATTERN = "phase8_commonlad*.jsonl"
CONTRAST_CVS = (0.0, 0.2)   # the CV levels the headline contrast needs a target at


def parse_committed(path: str) -> dict:
    """The per-arm means in the COMMITTED out/z_iso_measured.txt.

    Parsed rather than hardcoded, because a number typed into a grader stops tracking
    the artifact it claims to check.
    """
    if not os.path.exists(path):
        return {}
    txt = open(path, encoding="utf-8").read()
    block = txt.split("=== z_iso by friction ===")
    if len(block) < 2:
        return {}
    arms = {}
    for line in block[1].splitlines():
        m = re.match(r"\s*([\d.]+)\s+(\d+)\s+([\d.]+)\s+(\S+)", line)
        if m:
            arms[round(float(m.group(1)), 3)] = {"n": int(m.group(2)),
                                                 "z": float(m.group(3))}
    return arms


def arm_means(res: dict) -> dict:
    return {mu: {"n": len(v), "z": float(np.mean(v)),
                 "sd": float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")}
            for mu, v in res["per_mu"].items()}


def lognormal_table(res: dict) -> dict:
    """(cv -> z_iso) for the lognormal family at mu_p = 0.3, from one measurement."""
    return {cv: float(np.mean(v)) for (fam, cv), v in res["per_cv"].items()
            if fam == "lognormal"}


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    ok: dict[str, bool] = {}

    p("=== grade_ziso.py -- PREREG_ZISO_2026-09-12.md, gates G1-G7 ===")
    p("")

    committed = parse_committed(os.path.join(OUT, "z_iso_measured.txt"))
    p("committed artifact out/z_iso_measured.txt: "
      + ", ".join(f"mu{mu}: {d['z']:.3f} (n{d['n']})"
                  for mu, d in sorted(committed.items())))
    p("")

    # ---------------------------------------------------------------- G1
    p("=== G1  the committed artifact is STALE (predicted: yes) ===")
    stale = M.main(cap="all", rule="broken", dedup=False, quiet=True,
                   out=os.path.join(OUT, "z_iso_STALE_repro.txt"))
    sa = arm_means(stale)
    p("   re-running the SHIPPED behaviour (no cap filter, no dedup) on today's files:")
    for mu in sorted(sa):
        c = committed.get(mu)
        p(f"      mu_p {mu:4.2f}: {sa[mu]['z']:.3f} (n{sa[mu]['n']})"
          + (f"   committed {c['z']:.3f} (n{c['n']})   d = "
             f"{sa[mu]['z'] - c['z']:+.3f}" if c else "   [not in the artifact]"))
    d03 = abs(sa[0.3]["z"] - committed[0.3]["z"]) if 0.3 in committed else np.nan
    n1 = sa.get(1.0, {}).get("n", 0)
    ok["G1"] = bool(d03 > 0 and n1 > 1)
    p(f"   |d(mu 0.3)| = {d03:.3f} (must be > 0), n(mu 1.0) = {n1} (must be > 1) "
      f"-> {'PASS -- the glob now reaches the capfix shards' if ok['G1'] else 'FAIL'}")
    p(f"   and the shipped rule fits {M.N_LOW_RECORDS_BROKEN} RECORDS = 2 distinct "
      f"pressures per ladder (D2); the table above prints pts and rungs per cell.")
    p("")

    # ---------------------------------------------------------------- G2
    p(f"=== G2  reproduce before extending: --cap default + dedup, bar {G2_BAR} ===")
    base = M.main(cap="default", rule="broken", quiet=True,
                  out=os.path.join(OUT, "z_iso_default_broken.txt"))
    ba = arm_means(base)
    worst = 0.0
    for mu in (0.0, 0.3):
        if mu in committed and mu in ba:
            d = abs(ba[mu]["z"] - committed[mu]["z"])
            worst = max(worst, d)
            p(f"      mu_p {mu:4.2f}: {ba[mu]['z']:.4f} (n{ba[mu]['n']}) vs committed "
              f"{committed[mu]['z']:.4f} (n{committed[mu]['n']})   |d| = {d:.4f}")
    ok["G2"] = worst <= G2_BAR
    p(f"   worst |d| = {worst:.4f} (bar {G2_BAR}) -> "
      f"{'PASS' if ok['G2'] else 'FAIL -- the dedup changed the sample'}")
    p(f"   {base['info']['dupes']} duplicate keys dropped, of which "
      f"{len(base['info']['conflicts'])} disagree about phi -- the CV = 0 "
      f"non-determinism, and it lives in the FRICTIONLESS arm the mu_p = 0 move sits in")
    p("")

    # ---------------------------------------------------------------- G3
    p(f"=== G3  the fit rule as a systematic (bar: span < {G3_SPAN_BAR}) ===")
    runs = {"broken": base}
    for r in ("r2", "r3"):
        runs[r] = M.main(cap="default", rule=r, quiet=True,
                         out=os.path.join(OUT, f"z_iso_default_{r}.txt"))
    p(f"   {'rule':8s}{'mu 0':>9s}{'|-6|':>8s}{'mu 0.3':>9s}{'sd':>8s}{'n':>4s}"
      f"{'mu 3':>8s}{'|-4|':>8s}")
    for r in RULES_TESTED:
        a = arm_means(runs[r])
        z0, z3, z300 = a.get(0.0), a.get(0.3), a.get(3.0)
        p(f"   {r:8s}{z0['z']:9.3f}{abs(z0['z'] - Z_ISO_THEORY_FRICTIONLESS):8.3f}"
          f"{z3['z']:9.3f}{z3['sd']:8.3f}{z3['n']:4d}"
          f"{z300['z']:8.3f}{abs(z300['z'] - Z_ISO_THEORY_INFINITE_MU):8.3f}")
    prod = [arm_means(runs[r])[0.3]["z"] for r in RULES_TESTED]
    span = max(prod) - min(prod)
    ok["G3"] = span < G3_SPAN_BAR
    p(f"   span of z_iso(mu_p = 0.3) over the three rules = {span:.4f} "
      f"(bar {G3_SPAN_BAR}) -> {'PASS' if ok['G3'] else 'FAIL'}")
    p(f"   in phi, at the shear ladder's dz/dphi ~ 22: {span / 22:.4f} -- so D2 is a "
      f"real defect of the CODE whose numerical reach is {'small' if span < G3_SPAN_BAR else 'large'}")

    # the pre-registered SELECTION, decided by an anchor outside the code under test
    p("")
    p("   §4 selection by the mu_p = 0 anchor (z_iso = 2d = 6 exactly):")
    adm = {r: abs(arm_means(runs[r])[3.0]["z"] - Z_ISO_THEORY_INFINITE_MU) <= G3_MU3_BAR
           for r in ("r2", "r3")}
    dist = {r: abs(arm_means(runs[r])[0.0]["z"] - Z_ISO_THEORY_FRICTIONLESS)
            for r in ("r2", "r3")}
    for r in ("r2", "r3"):
        p(f"      {r}: |z_iso(0) - 6| = {dist[r]:.3f}, mu_p = 3 admissible = {adm[r]}")
    cands = [r for r in ("r2", "r3") if adm[r]]
    selected = min(cands, key=lambda r: dist[r]) if cands else None
    if selected is None:
        p("   🔴 NEITHER rung-correct rule is admissible at mu_p = 3; §6 applies")
        raise SystemExit(1)
    p(f"   -> SELECTED RULE: {selected!r}")
    p(f"   ⚠ §5's DIRECTIONAL prediction was |R2 - 6| < |R3 - 6| (a narrower pressure "
      f"range should suffer less curvature). Measured: R2 {dist['r2']:.3f} vs R3 "
      f"{dist['r3']:.3f} -> the prediction "
      f"{'HELD' if dist['r2'] < dist['r3'] else 'FAILED'}. With only four rungs the "
      f"two-rung fit is noise-limited, not curvature-limited (its sd at mu_p = 0 is "
      f"{arm_means(runs['r2'])[0.0]['sd']:.3f} against "
      f"{arm_means(runs['r3'])[0.0]['sd']:.3f}). The rule is selected by the anchor, "
      f"not by the prediction.")
    p("")

    # ---------------------------------------------------------------- G4
    p(f"=== G4  the CV drift survives the repair (bar |drift| >= {G4_DRIFT_BAR}, "
      f"sign negative) ===")
    sel = runs[selected]
    p(f"   {'rule':8s}{'family':12s}drift(cv0 -> cv_max)     table")
    drifts = {}
    for r in RULES_TESTED:
        for fam in ("lognormal", "uniform", "bidisperse"):
            t = {cv: float(np.mean(v)) for (f_, cv), v in runs[r]["per_cv"].items()
                 if f_ == fam}
            if len(t) < 2:
                continue
            cvs = sorted(t)
            dr = t[cvs[-1]] - t[cvs[0]]
            drifts[(r, fam)] = (dr, cvs[-1])
            p(f"   {r:8s}{fam:12s}{dr:+.3f} (cv0 -> cv{cvs[-1]:.1f})   "
              + " ".join(f"{cv:.1f}:{t[cv]:.3f}" for cv in cvs))
    dr_sel = drifts[(selected, "lognormal")][0]
    ok["G4"] = dr_sel < 0 and abs(dr_sel) >= G4_DRIFT_BAR
    p(f"   lognormal drift under {selected!r} = {dr_sel:+.3f} -> "
      f"{'PASS' if ok['G4'] else 'FAIL -- the CV-dependent bias was a fit artefact'}")
    p("   ⚠ NOT a general result: the drift is a LOGNORMAL effect. uniform and "
      "bidisperse collapse under the rung-correct rules, and uniform's broken-rule "
      "-0.193 rested entirely on its cv 0.5 cell, a ladder with only TWO pressures "
      "(cv 0 -> 0.4 it is -0.059/-0.058/-0.056, stable across all three rules). That "
      "the drift lives on the tail family and not on the zero-skew ones is the "
      "project's own CV/skewness separation, arrived at from a different observable.")
    p("")

    # ---------------------------------------------------------------- G5, G6
    p(f"=== G5/G6  propagation into the common ladder (bar |d(delta phi_c)| < "
      f"{G5_BAR}) ===")
    recs = [r for r in A.load(COMMON_PATTERN) if not r.get("error")]
    groups = cells_by_group(recs)
    p(f"   {len(recs)} usable common-ladder records, {len(groups)} groups")

    d_e4, *_ = contrast(groups, (), "E4")
    p(f"   E4 (theory target 4.0) through the new code path = {d_e4:+.8f} vs "
      f"published {PUBLISHED_E4:+.4f}; |d| = {abs(d_e4 - PUBLISHED_E4):.2e}")
    p(f"   🔴 G6's bar as pre-registered ({G6_TOL:.0e}) was UNACHIEVABLE, not failed: "
      f"out/common_ladder_grade.txt prints d to FOUR decimals, so no correct "
      f"implementation could match it to 1e-9. The failing number stays printed above.")
    p("   And G6 could not have tested what it claimed either way: E4 uses the THEORY "
      "target and never calls measure_z_iso, so both a correct and a broken refactor "
      "of that module pass it. The real refactor control is G2, which reproduces the "
      "committed artifact through the new selection path.")
    ok["G6a"] = abs(d_e4 - PUBLISHED_E4) <= G6A_BAR
    p(f"   G6a (replacement, at the artifact's own precision, bar {G6A_BAR:.0e}): "
      f"{'PASS' if ok['G6a'] else 'FAIL'}")

    old_tbl = dict(A.Z_ISO_MEASURED)
    d_e4m_old, means_old, _, _ = contrast(groups, (), "E4m")
    new_lognormal = lognormal_table(sel)
    patched = dict(old_tbl)
    for (fam, cv) in list(patched):
        if fam == "lognormal" and round(cv, 3) in \
                {round(c, 3) for c in new_lognormal}:
            patched[(fam, cv)] = new_lognormal[
                [c for c in new_lognormal if abs(c - cv) < 1e-9][0]]
    p("   z_iso table, lognormal, mu_p = 0.3:")
    for cv in sorted(new_lognormal):
        p(f"      cv {cv:.1f}: {old_tbl.get(('lognormal', cv), float('nan')):.3f} "
          f"-> {new_lognormal[cv]:.3f}   "
          f"({new_lognormal[cv] - old_tbl.get(('lognormal', cv), np.nan):+.3f})")
    try:
        A.Z_ISO_MEASURED.clear()
        A.Z_ISO_MEASURED.update(patched)
        d_e4m_new, means_new, _, _ = contrast(groups, (), "E4m")
    finally:
        A.Z_ISO_MEASURED.clear()
        A.Z_ISO_MEASURED.update(old_tbl)
    # G6b, the replacement control, and it CAN fail: the patch must be reversible, or
    # every number computed after this point is silently on the other table.
    d_e4m_back, *_ = contrast(groups, (), "E4m")
    ok["G6b"] = abs(d_e4m_back - d_e4m_old) <= 1e-12 and \
        A.Z_ISO_MEASURED == old_tbl
    p(f"   G6b control: after restoring the table, E4m returns {d_e4m_back:+.8f} "
      f"against {d_e4m_old:+.8f} -> {'PASS' if ok['G6b'] else 'FAIL -- THE PATCH LEAKED'}")
    move = d_e4m_new - d_e4m_old
    ok["G5"] = abs(move) < G5_BAR
    p(f"   E4m delta phi_c: {d_e4m_old:+.4f} (published {PUBLISHED_E4M:+.4f}) "
      f"-> {d_e4m_new:+.4f}   move {move:+.4f} (bar {G5_BAR}) -> "
      f"{'PASS -- the headline stands' if ok['G5'] else 'FAIL -- REQUALIFY the headline'}")
    p("   per level:")
    for cv in sorted(means_old):
        p(f"      cv {cv:.1f}: {means_old[cv]:.4f} -> {means_new[cv]:.4f}  "
          f"({means_new[cv] - means_old[cv]:+.4f})")
    p("   the arithmetic this was predicted from: a UNIFORM shift of the target "
      "cancels out of a contrast exactly, so only the CV-dependence of the shift "
      "survives, and G4 says that dependence is the drift already in the table.")
    p("")

    # ---------------------------------------------------------------- G7
    p("=== G7  the ADOPTION decision, rule fixed in §5 before the run ===")
    rep = M.main(cap="repaired", rule=selected, quiet=True,
                 out=os.path.join(OUT, f"z_iso_repaired_{selected}.txt"))
    cov = sorted({(fam, cv) for (fam, cv) in rep["per_cv"]})
    p(f"   the repaired cap measures z_iso at mu_p = 0.3 for: {cov or 'NOTHING'}")
    missing = [cv for cv in CONTRAST_CVS
               if not any(abs(cv - c) < 1e-9 for (_f, c) in cov)]
    p(f"   the headline contrast needs a target at CV {list(CONTRAST_CVS)} + 0.5; "
      f"missing from the repaired cap: {missing}")
    ok["G7"] = bool(missing)      # the pre-registered outcome is DO NOT ADOPT
    p(f"   -> {'DO NOT ADOPT' if missing else 'ADOPT IS AVAILABLE'}: "
      f"{'a half-adopted table would be default-cap at CV 0/0.2 and repaired-cap at CV 0.4/0.5, putting the ITERATION CAP on the CV axis under test' if missing else 'every CV level is covered'}")
    p("   PHI_J_ISO is not adopted under any outcome: it records where the shear "
      "ladders were physically placed, not an estimate of anything (D4).")
    p("")

    # ---------------------------------------------------------------- summary
    p("=== SUMMARY ===")
    for g in ("G1", "G2", "G3", "G4", "G5", "G6a", "G6b", "G7"):
        p(f"   {g}: {'PASS' if ok.get(g) else 'FAIL'}")
    n_fail = sum(1 for v in ok.values() if not v)
    p(f"   {len(ok) - n_fail} of {len(ok)} gates pass")

    path = os.path.join(OUT, "ziso_audit.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

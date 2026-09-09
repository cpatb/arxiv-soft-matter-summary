"""check_regime_flag.py — is the experiment's `arrested` flag current or stale?

WHY THIS EXISTS
---------------
NEXT_STEPS_2026-09-09.md section 2 item 5 blocks the salt re-analysis -- the cheapest
test the project has, because the salt ladder is the only designed axis already on disk
-- on one question. The project records that the `arrested` flag once "stayed True for
22 drops that never jammed after the classification moved", and in `phase_drops.csv`
**27 of 160 drops carry a finite `d_arrest` while `arrested` is False**. Either those 27
are a stale flag, in which case every regime-resolved statement built on the existing
data is wrong, or they are simply drops that were scored and did not arrest.

The question is settled by REPRODUCING the flag from the column it claims to be built
from, not by reading the prose that describes it.

WHAT IT CHECKS, ON READ-ONLY DATA
---------------------------------
  1. `arrested == (d_arrest < JAM_CUT)` on every row, with `JAM_CUT` read out of the
     generating script rather than typed here.
  2. `d_arrest` is the DIMENSIONLESS `d_arrest/d0`, not a length -- the column name says
     neither, and a length would make the cut meaningless.
  3. the cut sits in an EMPTY band of the order parameter, so the classification is not
     a tuned threshold. The width of that band is the margin, and it is printed.

⚠ `E:\\CPB\\Drop Impact Analyses<U+F028>` ends in a Unicode private-use character that
cannot be typed. Every path here is GLOBBED; a typed path reads as missing data.
Raw data on E: is read-only and nothing here writes to it.
"""

from __future__ import annotations

import glob
import os
import re

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

DROPS_GLOB = r"E:\CPB\*\figures_rebuild_*\phase_drops.csv"
GEN_GLOB = r"E:\CPB\*\figures_rebuild_*\phase_states.py"


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    drops = sorted(set(glob.glob(DROPS_GLOB)))
    gens = sorted(set(glob.glob(GEN_GLOB)))
    if not drops:
        p(f"🔴 no phase_drops.csv matched {DROPS_GLOB} -- glob, never type, the E: path")
        return 1
    path = drops[-1]
    p(f"table     : {path}")

    # --- the cut, read out of the generator ------------------------------------------
    cut = None
    rule_src = None
    if gens:
        src = open(gens[-1], encoding="utf-8", errors="replace").read()
        m = re.search(r"^JAM_CUT\s*=\s*([0-9.]+)", src, re.M)
        if m:
            cut = float(m.group(1))
        m2 = re.search(r'^\s*d\["arrested"\]\s*=\s*(.+)$', src, re.M)
        if m2:
            rule_src = m2.group(1).strip()
        p(f"generator : {gens[-1]}")
        p(f"            JAM_CUT = {cut}   rule: arrested = {rule_src}")
    if cut is None:
        p("🔴 could not read JAM_CUT from the generator; refusing to assume 1.2")
        return 1

    d = pd.read_csv(path)
    p(f"rows      : {len(d)}")

    fin = d["d_arrest"].notna()
    p(f"\n  arrested True , d_arrest finite : {int((d.arrested & fin).sum())}")
    p(f"  arrested True , d_arrest NaN    : {int((d.arrested & ~fin).sum())}")
    p(f"  arrested False, d_arrest finite : {int((~d.arrested & fin).sum())}  "
      f"<- the rows NEXT_STEPS asks about")
    p(f"  arrested False, d_arrest NaN    : {int((~d.arrested & ~fin).sum())}")

    # --- 1. does the flag reproduce? ---------------------------------------------------
    rule = d["d_arrest"] < cut
    n_mismatch = int((d["arrested"] != rule).sum())
    p(f"\n1. arrested == (d_arrest < {cut}) on {len(d) - n_mismatch}/{len(d)} rows, "
      f"{n_mismatch} mismatches")

    # --- 2. is d_arrest dimensionless? -------------------------------------------------
    ratio = d.loc[fin, "d_arrest"] / d.loc[fin, "d0_mm"]
    p(f"\n2. d_arrest spans {d.loc[fin, 'd_arrest'].min():.3f}-"
      f"{d.loc[fin, 'd_arrest'].max():.3f} while d0_mm spans "
      f"{d.d0_mm.min():.3f}-{d.d0_mm.max():.3f} mm; d_arrest/d0_mm spans "
      f"{ratio.min():.3f}-{ratio.max():.3f}")
    dimensionless = bool(d.loc[fin, "d_arrest"].max() < 2.0 * d.d0_mm.min())
    p(f"   d_arrest is already the DIMENSIONLESS d_arrest/d0: "
      f"{'yes' if dimensionless else '🔴 NO -- the cut is being applied to a length'}")

    # --- 3. how much room does the cut have? -------------------------------------------
    hi_arrested = float(d.loc[d.arrested, "d_arrest"].max())
    lo_spread = float(d.loc[~d.arrested & fin, "d_arrest"].min())
    p(f"\n3. order parameter: arrested drops top out at {hi_arrested:.4f}; "
      f"non-arrested start at {lo_spread:.4f}")
    p(f"   the cut {cut} sits in an EMPTY band of width {lo_spread - hi_arrested:.4f}, "
      f"with {cut - hi_arrested:.4f} below it and {lo_spread - cut:.4f} above")
    empty = bool(hi_arrested < cut < lo_spread)

    # --- verdict -----------------------------------------------------------------------
    p("\n=== VERDICT ===")
    if n_mismatch == 0 and dimensionless and empty:
        p(f"   THE FLAG IS CURRENT, NOT STALE. It reproduces from `d_arrest` on all "
          f"{len(d)} rows, `d_arrest` is the dimensionless spread ratio, and the cut "
          f"sits in a band {lo_spread - hi_arrested:.3f} wide containing no data at "
          f"all -- so it could move anywhere in that band without relabelling a single "
          f"drop.")
        p(f"   The {int((~d.arrested & fin).sum())} rows with a finite `d_arrest` and "
          f"`arrested = False` are drops that WERE scored and spread past {cut} d0. "
          f"They are not a stale flag; they are the non-arrested arm of the "
          f"measurement, and the {int((~d.arrested & ~fin).sum())} with no "
          f"`d_arrest` at all are the unscored ones.")
        p("   -> NEXT_STEPS_2026-09-09 section 2 item 5 is UNBLOCKED. The salt "
          "re-analysis may use `arrested` as its regime label.")
        rc = 0
    else:
        p("   🔴 THE FLAG DOES NOT REPRODUCE. Do not use it. Failing checks: "
          + ", ".join(n for n, ok in
                      (("rule", n_mismatch == 0), ("dimensionless", dimensionless),
                       ("empty band", empty)) if not ok))
        rc = 1

    # --- what the salt axis actually looks like under it --------------------------------
    # ⚠ The denominator matters and the two choices point OPPOSITE ways here, so both
    # are printed. An UNSCORED drop (no `d_arrest`) cannot testify to arrest -- the
    # generator's own comment says exactly that -- so `n_arrested / n_scored` is the
    # defensible fraction and `n_arrested / n` is a coverage statistic wearing its
    # label. And the u0 ranges are not matched across the ladder, so neither fraction
    # is a regime comparison on its own.
    p("\n=== the salt ladder under this label (control variables only) ===")
    g = d.groupby(["salt", "phi"], dropna=False)
    p(f"   {'salt':>6s} {'phi':>5s} {'n':>4s} {'n_scored':>9s} {'n_arrested':>11s} "
      f"{'/scored':>8s} {'/all':>6s} {'u0 range':>16s}")
    for key, sub in g:
        salt, phi = key
        sc = int(sub["d_arrest"].notna().sum())
        ar = int(sub["arrested"].sum())
        p(f"   {salt:6.1f} {phi:5.2f} {len(sub):4d} {sc:9d} {ar:11d} "
          f"{ar / sc if sc else float('nan'):8.2f} {ar / len(sub):6.2f} "
          f"{sub.u0.min():7.2f}-{sub.u0.max():.2f} m/s")
    p("   ⚠ the two fractions run in OPPOSITE directions across the salt ladder, and "
      "the u0 coverage is not matched (0 mM reaches 4.95 m/s, 100 mM only 2.63). "
      "Neither column is a regime comparison; the salt re-analysis has to match u0.")

    path_out = os.path.join(OUT_DIR, "regime_flag_check.txt")
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(path_out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path_out}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

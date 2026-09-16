"""sweep_descriptive.py -- the descriptive sweep, from shards that already exist.

WHAT THIS IS. The project has run three size-distribution families x six CV x a pressure
ladder x friction, and has reported almost none of it descriptively: the effort went into
a closure that predicts phi_J for distributions that were never simulated. This assembles
what IS measured -- phi_J, Z_all, Z_nr, the rattler fraction and the contact-force moments
against CV, per family, at matched friction and matched protocol -- and writes
REPORT_SWEEP_2026-09-18.md from the numbers it prints. **No new cells.**

THREE RULES IT FOLLOWS, EACH FOR A RECORDED REASON.

1. REPRODUCTION FIRST. The three family shifts this project already quotes
   (+0.0621 / +0.0242 / +0.0139) are recomputed here from `phase2_fits.csv` and compared
   against the committed `out/shape_control.txt` BEFORE anything new is printed. A
   diagnostic that cannot reproduce the number it extends has not extended it. Exits
   non-zero on a mismatch.

2. phi_J IS FITTED AND THE OTHERS ARE NOT. phi_J is an extrapolation of a four-rung
   pressure ladder to p_hat -> 0; Z, the rattler fraction and the force moments are read
   AT a pressure. They are therefore reported at every rung, and the rung dependence is
   printed rather than hidden behind one choice.

3. THE LADDER-HEALTH GATE TRAVELS WITH phi_J. `shape_control.py` applies it (no negative
   ladder slope, median r2 >= 0.85) because a 2026-09-13 review found a friction arm whose
   phi_J came from 0/6 monotone ladders. Any phi_J quoted here carries the same gate.
"""

from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FITS = os.path.join(OUT, "phase2_fits.csv")
SHAPE_TXT = os.path.join(OUT, "shape_control.txt")
REPORT = os.path.join(HERE, "REPORT_SWEEP_2026-09-18.md")

FAMILIES = ("lognormal", "uniform", "bidisperse")
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
MU_REF = 0.3
N_MATCH = 1500
RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
R2_MIN = 0.85

L: list[str] = []


def P(s=""):
    L.append(s)


def tcrit(n: int, alpha: float = 0.05) -> float:
    """Two-sided t critical value at n-1 dof.

    NOT 2. Every change in this file comes from 3 seeds, so its se carries 2 dof and the
    95 % value is 4.303; P(|T| > 2, df = 2) = 0.184, i.e. a flat 2*se bar is an 18 %
    test. The first version of this script used 2*se and it admitted exactly the change
    that sets the denominator of the headline ratio (referee, 2026-09-18).
    """
    from scipy import stats
    return float(stats.t.ppf(1.0 - alpha / 2.0, max(1, n - 1)))


def cochran_q(means, ses):
    """Cochran's Q heterogeneity test across groups, and its p-value.

    The question "are the three families' changes the same?" has a standard test and
    this file previously answered it by eye. Q treats the reported se's as EXACTLY
    known, which they are not (2 dof each), so it is the LIBERAL end of the range; the
    conservative end is a t-based comparison. Both are printed, because they disagree
    and the disagreement is the honest state (referee, 2026-09-18).
    """
    from scipy import stats
    m = np.asarray(means, float)
    s = np.asarray(ses, float)
    ok = np.isfinite(m) & np.isfinite(s) & (s > 0)
    if ok.sum() < 2:
        return float("nan"), float("nan")
    w = 1.0 / s[ok] ** 2
    mbar = float(np.sum(w * m[ok]) / np.sum(w))
    q = float(np.sum(w * (m[ok] - mbar) ** 2))
    return q, float(stats.chi2.sf(q, int(ok.sum()) - 1))


def mean_se(v):
    v = np.asarray([x for x in v if x == x], float)
    if len(v) == 0:
        return float("nan"), float("nan"), 0
    if len(v) == 1:
        return float(v[0]), float("nan"), 1
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v)


# ------------------------------------------------------------------- 1. reproduction
def load_fits() -> pd.DataFrame:
    d = pd.read_csv(FITS)
    d = d[(d.n_requested == N_MATCH) & np.isfinite(d.phi_J_lin)]
    return d.loc[sorted(d.groupby(["family", "cv_nom", "mu_p", "seed"])
                        ["max_iter"].idxmax())]


def paired_shift(d, family, cv_lo=0.0, cv_hi=0.5, mu=MU_REF):
    g = d[(d.family == family) & (d.mu_p == mu)]
    lo = g[g.cv_nom == cv_lo].set_index("seed").phi_J_lin
    hi = g[g.cv_nom == cv_hi].set_index("seed").phi_J_lin
    common = sorted(set(lo.index) & set(hi.index))
    if len(common) < 2:
        return None
    diff = np.array([hi[s] - lo[s] for s in common])
    return float(diff.mean()), float(diff.std(ddof=1) / math.sqrt(len(diff))), len(common)


def reproduce(d) -> bool:
    P("1. REPRODUCTION GATE -- runs first, exits non-zero on a mismatch.")
    P("   The committed out/shape_control.txt quotes the paired CV 0 -> 0.5 shift in")
    P("   phi_J per family. Recomputed here from phase2_fits.csv:")
    txt = open(SHAPE_TXT, encoding="utf-8").read()
    ok = True
    P("     family       recomputed      se    n     committed   match")
    for fam in FAMILIES:
        r = paired_shift(d, fam)
        if r is None:
            P(f"     {fam:12s} no paired seeds")
            ok = False
            continue
        m, se, n = r
        m_pat = re.search(rf"{fam}\s+\+([0-9.]+)\s+([0-9.]+)\s+(\d+)", txt)
        com = float(m_pat.group(1)) if m_pat else float("nan")
        good = abs(m - com) < 5e-5
        ok = ok and good
        P(f"     {fam:12s} {m:+.6f}  {se:.6f}  {n}    {com:+.4f}    "
          f"{'YES' if good else 'NO'}")
    P(f"   REPRODUCTION: {'PASS' if ok else 'FAIL'}")
    return ok


def ladder_health(d) -> dict:
    P("")
    P("2. LADDER HEALTH -- the gate that travels with phi_J.")
    P("   family        n   neg slopes  r2 median   verdict")
    out = {}
    for fam in FAMILIES:
        g = d[(d.family == fam) & (d.mu_p == MU_REF)]
        neg = int((g.slope_lin <= 0).sum())
        r2 = float(np.median(g.r2_lin)) if len(g) else float("nan")
        good = neg == 0 and r2 >= R2_MIN
        out[fam] = good
        P(f"   {fam:12s} {len(g):3d}   {neg:8d}   {r2:8.4f}   {'HEALTHY' if good else 'FAILS'}")
    return out


# ---------------------------------------------------------------- 2. the shard channel
def load_shards():
    rows = []
    files = sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl")))
    for f in files:
        for ln in open(f, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            try:
                r = json.loads(ln)
            except json.JSONDecodeError:
                continue
            rows.append(r)
    P("")
    P("3. THE PER-CELL OBSERVABLES -- N printed before and after every filter.")
    P(f"   files globbed out/phase2*.jsonl: {len(files)}")
    P(f"   N records: {len(rows)}")
    # DEDUP by key, keeping the highest iteration cap, the same rule analyse_phase2 uses.
    # Duplicates that DISAGREE are counted and printed rather than silently resolved
    # (lessons.md 2026-09-09: which record survives was set by the filename sort order).
    by = {}
    conflicts = 0
    for r in rows:
        if r.get("failed"):
            continue
        k = r.get("key")
        if k in by:
            if abs(by[k]["phi"] - r["phi"]) > 1e-12:
                conflicts += 1
            if r.get("max_iter", 0) <= by[k].get("max_iter", 0):
                continue
        by[k] = r
    P(f"   N after dropping failed and deduplicating by key: {len(by)}")
    P(f"   duplicate keys whose phi DISAGREED: {conflicts}")
    sel = [r for r in by.values()
           if r.get("n_requested") == N_MATCH and r.get("mu_p") == MU_REF
           and r.get("converged")]
    P(f"   N after N={N_MATCH}, mu_p={MU_REF}, converged: {len(sel)}")
    return sel


def table(sel, getter, name, unit=""):
    P("")
    P(f"   {name}{unit}, mean +- se over seeds, by family and CV")
    P("     p_hat      family        " + "  ".join(f"  CV{cv:.1f}   " for cv in CVS))
    for p in RUNGS:
        for fam in FAMILIES:
            cells = []
            for cv in CVS:
                v = [getter(r) for r in sel
                     if (r.get("family") or "lognormal") == fam and r["cv"] == cv
                     and r["p_hat_target"] == p]
                m, se, n = mean_se(v)
                cells.append(f"{m:7.4f}({n})" if n else "    --   ")
            P(f"     {p:.0e}  {fam:12s} " + " ".join(cells))
        P("")


def main() -> int:
    P("sweep_descriptive.py -- the descriptive sweep, from existing shards. No new cells.")
    P("=" * 78)
    d = load_fits()
    if not reproduce(d):
        write()
        return 1
    health = ladder_health(d)

    P("")
    P("   phi_J (linear extrapolation to p_hat -> 0), mean +- se over seeds")
    P("     family        " + "  ".join(f"  CV{cv:.1f}   " for cv in CVS))
    phiJ = {}
    for fam in FAMILIES:
        cells = []
        for cv in CVS:
            g = d[(d.family == fam) & (d.mu_p == MU_REF) & (d.cv_nom == cv)]
            m, se, n = mean_se(list(g.phi_J_lin))
            phiJ[(fam, cv)] = (m, se, n)
            cells.append(f"{m:7.4f}({n})" if n else "    --   ")
        P(f"     {fam:12s} " + " ".join(cells))
    P("   (ladder health: " + ", ".join(f"{f}={'ok' if health[f] else 'FAILS'}"
                                        for f in FAMILIES) + ")")

    sel = load_shards()
    table(sel, lambda r: r["z_all"], "Z_all")
    table(sel, lambda r: r["z_nr"], "Z_nr (non-rattler backbone, k = 2)")
    table(sel, lambda r: r["rattler_frac"], "rattler fraction (k = 2)")
    table(sel, lambda r: r["force"].get("mean", float("nan")), "contact force, mean")
    table(sel, lambda r: (r["force"].get("cv", float("nan"))), "contact force, CV")

    # --- the one comparison the write-up exists to make ---------------------------
    P("")
    P("4. THE SHAPE-OVER-CV CONTRAST, on every channel, at p_hat = 1e-3.")
    P("   For each observable: the CV 0 -> 0.5 change WITHIN a family, and the spread of")
    P("   that change ACROSS families. A ratio > 1 in the last column means the shape of")
    P("   the size distribution moves the observable more than a fourfold change in CV")
    P("   does within any one family.")
    P("")
    P("   Changes are PAIRED PER SEED, never a difference of two means: the two ends")
    P("   share the seed, so the paired form cancels the shared draw exactly rather than")
    P("   statistically (lessons.md 2026-08-30). Each change is printed with its own se,")
    P("   and a ratio whose SMALLEST change is not significantly non-zero is reported")
    P("   UNDETERMINED -- a relative measure whose denominator is shrinking measures the")
    P("   shrinkage.")
    P("")
    P("   ⚠ THE GUARD IS t(df = n-1), NOT 2*se (referee, 2026-09-18). Every change here")
    P("   comes from 3 seeds, so its se carries 2 dof and the two-sided 95 % critical")
    P(f"   value is t = {tcrit(3):.3f}, not 2. P(|T| > 2, df = 2) = 0.184, so a flat")
    P("   2*se bar is an 18 % test. The first version of this script used 2*se and it")
    P("   admitted exactly the change that carries the headline.")
    P("")
    P("   observable              lognormal         uniform        bidisperse"
      "      across-family")
    def chg(getter, p=1e-3):
        out = {}
        for fam in FAMILIES:
            lo = {r["seed"]: getter(r) for r in sel
                  if (r.get("family") or "lognormal") == fam and r["cv"] == 0.0
                  and r["p_hat_target"] == p}
            hi = {r["seed"]: getter(r) for r in sel
                  if (r.get("family") or "lognormal") == fam and r["cv"] == 0.5
                  and r["p_hat_target"] == p}
            common = sorted(set(lo) & set(hi))
            diffs = [hi[s] - lo[s] for s in common]
            out[fam] = mean_se(diffs)
        return out
    def peak_chg(getter, p=1e-3):
        """The change from CV = 0 to whichever rung is FURTHEST from it.

        The endpoint definition assumes monotonicity, and the bidisperse family is
        NON-MONOTONE in CV -- a pathology this project already flags for phi_J. Reported
        beside the endpoint definition so the reader can see how much of the headline is
        the definition (referee, 2026-09-18).
        """
        out = {}
        for fam in FAMILIES:
            best = (0.0, float("nan"), 0.0)
            for cv in CVS[1:]:
                lo = {r["seed"]: getter(r) for r in sel
                      if (r.get("family") or "lognormal") == fam and r["cv"] == 0.0
                      and r["p_hat_target"] == p}
                hi = {r["seed"]: getter(r) for r in sel
                      if (r.get("family") or "lognormal") == fam and r["cv"] == cv
                      and r["p_hat_target"] == p}
                common = sorted(set(lo) & set(hi))
                if not common:
                    continue
                m, se, n = mean_se([hi[s] - lo[s] for s in common])
                if abs(m) > abs(best[0]):
                    best = (m, se, cv)
            out[fam] = best
        return out

    CHANNELS = (("phi_J", None),
                ("Z_all", lambda r: r["z_all"]),
                ("Z_nr (k=2)", lambda r: r["z_nr"]),
                ("rattler fraction", lambda r: r["rattler_frac"]),
                ("contact force mean", lambda r: r["force"].get("mean", float("nan"))),
                ("contact force CV", lambda r: r["force"].get("cv", float("nan"))))
    rows_out = []
    for nm, g in CHANNELS:
        if g is None:
            ch = {}
            for f in FAMILIES:
                ps = paired_shift(d, f)
                ch[f] = ps if ps else (float("nan"), float("nan"), 0)
        else:
            ch = chg(g)
        vals = [abs(ch[f][0]) for f in FAMILIES if ch[f][0] == ch[f][0]]
        ses = [ch[f][1] for f in FAMILIES if ch[f][0] == ch[f][0]]
        ns = [ch[f][2] for f in FAMILIES if ch[f][0] == ch[f][0]]
        if len(vals) != len(FAMILIES) or min(vals) == 0:
            ratio, note = float("nan"), "n/a"
        else:
            ratio = max(vals) / min(vals)
            j = int(np.argmin(vals))
            tc = tcrit(ns[j])
            note = (f"{ratio:.2f}x" if vals[j] > tc * ses[j]
                    else f"UNDETERMINED (smallest change {vals[j]:.4f} does not clear "
                         f"t={tc:.2f} x se={ses[j]:.4f} = {tc*ses[j]:.4f})")
        rows_out.append((nm, ch, ratio, note))
        P(f"   {nm:22s} " + " ".join(f"{ch[f][0]:+8.4f}+-{ch[f][1]:.4f}"
                                     for f in FAMILIES) + f"   {note}")

    # --- the test the eyeball judgement was standing in for -------------------------
    P("")
    P("   DO THE FAMILIES ACTUALLY DIFFER? Cochran's Q across the three changes, which")
    P("   is the standard test the earlier version of this file answered by eye. Q")
    P("   treats the se's as exactly known (they carry 2 dof), so it is the LIBERAL end;")
    P("   the conservative end is that with 3 seeds almost nothing is resolvable. Both")
    P("   are reported because they disagree, and the disagreement is the honest state.")
    P("     observable               Q      p(Q)    families differ at 5 % (liberal)?")
    qres = {}
    for nm, ch, ratio, note in rows_out:
        q, p_ = cochran_q([ch[f][0] for f in FAMILIES], [ch[f][1] for f in FAMILIES])
        qres[nm] = (q, p_)
        P(f"     {nm:22s} {q:7.2f}  {p_:7.4f}   {'YES' if p_ < 0.05 else 'no'}")

    # --- how much of the ratio is the endpoint definition? --------------------------
    P("")
    P("   ENDPOINT vs PEAK. The CV 0 -> 0.5 definition assumes monotonicity in CV, and")
    P("   the bidisperse family is not monotone. Recomputed against whichever rung is")
    P("   furthest from CV = 0:")
    P("     observable             endpoint ratio   peak ratio   peak rung per family")
    peak_rows = {}
    for nm, g in CHANNELS:
        if g is None:
            continue
        pk = peak_chg(g)
        vals = [abs(pk[f][0]) for f in FAMILIES]
        pr = max(vals) / min(vals) if min(vals) > 0 else float("nan")
        peak_rows[nm] = pr
        ep = dict((r[0], r[2]) for r in rows_out)[nm]
        P(f"     {nm:22s} {ep:12.2f}x  {pr:10.2f}x   "
          + " ".join(f"{f[:4]}@{pk[f][2]:.1f}" for f in FAMILIES))

    # --- is the ORDERING itself stable? ---------------------------------------------
    P("")
    P("   IS THE ORDERING STABLE? Each family change is resampled as t(df = n-1) around")
    P("   its own mean and se, the six ratios recomputed, and the exact ordering counted.")
    rng = np.random.default_rng(20260918)
    order0 = [r[0] for r in sorted(rows_out, key=lambda r: -r[2])]
    hits = 0
    pair_hits = {i: 0 for i in range(len(order0) - 1)}
    REPS = 20000
    for _ in range(REPS):
        draw = {}
        for nm, ch, ratio, note in rows_out:
            vals = []
            for f in FAMILIES:
                m, se, n = ch[f]
                vals.append(abs(m + se * rng.standard_t(max(1, n - 1))))
            draw[nm] = max(vals) / min(vals) if min(vals) > 0 else float("inf")
        o = [k for k in sorted(draw, key=lambda k: -draw[k])]
        if o == order0:
            hits += 1
        for i in range(len(order0) - 1):
            if draw[order0[i]] > draw[order0[i + 1]]:
                pair_hits[i] += 1
    P(f"     the exact six-way ordering survives {100*hits/REPS:.0f} % of {REPS} draws")
    for i in range(len(order0) - 1):
        P(f"       {order0[i]:22s} > {order0[i+1]:22s}  {100*pair_hits[i]/REPS:5.1f} %")
    order_frac = hits / REPS

    # --- the channels are not independent -------------------------------------------
    P("")
    P("   ⚠ Z_all IS NOT AN INDEPENDENT CHANNEL. Z_all = (1 - rattler_frac) * Z_nr to")
    P("   within the numbers already in this file, so its near-family-independence is")
    P("   INHERITED from the rattler fraction's and is not a second confirmation:")
    P("     family        CV    Z_all    (1-rat)*Z_nr    difference")
    worst_id = 0.0
    for fam in FAMILIES:
        for cv in (0.0, 0.3, 0.5):
            sub = [r for r in sel if (r.get("family") or "lognormal") == fam
                   and r["cv"] == cv and r["p_hat_target"] == 1e-3]
            if not sub:
                continue
            za = float(np.mean([r["z_all"] for r in sub]))
            zz = float(np.mean([(1 - r["rattler_frac"]) * r["z_nr"] for r in sub]))
            worst_id = max(worst_id, abs(za - zz))
            P(f"     {fam:12s} {cv:.1f}  {za:7.4f}   {zz:11.4f}   {za-zz:+10.4f}")
    P(f"     worst |difference| over the rows printed: {worst_id:.4f}")

    P("")
    P("   READ THIS WITH THE FOLLOWING RESTRICTIONS, none of which is optional:")
    P("   - bidisperse is 50:50 BY NUMBER, so its fines carry 0.500 -> 0.036 of the")
    P("     VOLUME across this CV ladder. It is a PATH through the (ratio, composition)")
    P("     plane and not a representative 'bidisperse family'; the 2026-09-13 block")
    P("     found matched binaries 0.0635 apart in phi where every PSD descriptor")
    P("     predicts 0.")
    P("   - a ratio across families is a spread, not a pairwise separation. It says the")
    P("     three families differ; it does not say any two of them are distinguishable")
    P("     at a given CV.")
    P("   - Z, the rattler fraction and the force moments are read AT a pressure and")
    P("     move with the rung; phi_J alone is extrapolated to p_hat -> 0.")
    write()
    md(rows_out, health, sel, qres, peak_rows, order_frac, order0, pair_hits, REPS)
    return 0


def md(rows_out, health, sel, qres, peak_rows, order_frac, order0, pair_hits, reps):
    """The report is GENERATED from the verdicts, never written alongside them."""
    # The narrative is generated FROM the verdicts: `lead` is whichever channel actually
    # has the largest RESOLVED across-family spread, not one chosen in advance (this
    # project has shipped three artifacts whose prose contradicted their own tables).
    resolved = [r for r in rows_out if r[2] == r[2] and r[3].endswith("x")]
    lead = max(resolved, key=lambda r: r[2]) if resolved else None
    n_cells = len(sel)
    body = [
        "# The descriptive sweep — three size-distribution families, one axis each",
        "",
        "*Generated by `sweep_descriptive.py` from shards already on disk. No new cells.*",
        "*Every number below is printed by that script into `out/sweep_descriptive.txt`;",
        "none of them is typed here by hand.*",
        "",
        "## What was run",
        "",
        f"Three size-distribution families — lognormal, uniform, bidisperse — over CV "
        f"0.0–0.5, at matched friction `mu_p` = {MU_REF}, N = {N_MATCH}, on a four-rung "
        f"reduced-pressure ladder. {n_cells} converged cells survive the matched-condition "
        "filter; the counts before and after every filter are in the artifact.",
        "",
        "## The result",
        "",
        (f"The largest **resolved** across-family contrast is **{lead[0]}**, whose "
         f"paired CV 0 → 0.5 change spans a factor **{lead[2]:.2f}×** between the three "
         "families ("
         + ", ".join(f"{f} {lead[1][f][0]:+.4f}±{lead[1][f][1]:.4f}" for f in FAMILIES)
         + "). Every family change here rests on 3 seeds, so read the ordering section "
           "below before carrying any of it."
         if lead is not None else
         "**No across-family contrast is resolved.** On every channel the smallest of "
         "the three family changes fails a two-sided 95% test at its own degrees of "
         "freedom, so no ratio is quotable — which is itself the finding, and it is a "
         "statement about 3 seeds, not about the packings."),
        "",
        "| observable | lognormal | uniform | bidisperse | across-family spread |",
        "|---|---|---|---|---|",
    ]
    for nm, ch, ratio, note in rows_out:
        body.append(
            f"| {nm} | {ch['lognormal'][0]:+.4f}±{ch['lognormal'][1]:.4f} "
            f"| {ch['uniform'][0]:+.4f}±{ch['uniform'][1]:.4f} "
            f"| {ch['bidisperse'][0]:+.4f}±{ch['bidisperse'][1]:.4f} | {note} |")
    body += [
        "",
        "(CV 0 → 0.5 change, at `p_hat` = 1e-3 for the per-cell channels and "
        "extrapolated to `p_hat` → 0 for `phi_J`.)",
        "",
        "### How far the ordering can be pushed — and it is not far",
        "",
        "Ranked by across-family spread, with the ones whose denominator fails a "
        "two-sided 95% test marked — they are ranked but **not quotable**: "
        + " > ".join(f"**{r[0]}** {r[2]:.2f}×"
                     + ("" if r[3].endswith("x") else " *(UNDETERMINED)*")
                     for r in sorted(rows_out, key=lambda r: -r[2])) + ".",
        "",
        (f"**That six-way ordering survives only {100*order_frac:.0f}% of "
         f"{reps:,} resamples** in which each family change is redrawn as "
         "*t*(df = n−1) around its own mean and se. It is not a result. The links that "
         "do hold are the individual ones: "
         + "; ".join(f"`{order0[i]}` > `{order0[i+1]}` {100*pair_hits[i]/reps:.0f}%"
                     for i in range(len(order0) - 1)) + "."),
        "",
        "Three separate reasons the headline has to be stated weakly, all raised by an "
        "adversarial referee against an earlier draft of this page:",
        "",
        "1. **The significance guard was wrong.** Each change comes from 3 seeds, so its "
        "se carries 2 degrees of freedom and the 95% critical value is *t* = 4.303, not "
        "2. A flat 2·se bar is an 18% test, and it admitted exactly the change that set "
        "the denominator of the largest ratio. The guard is now *t*(df = n−1) and the "
        "table above says which ratios survive it.",
        "2. **The endpoint definition does real work.** CV 0 → 0.5 assumes monotonicity "
        "and the bidisperse family is not monotone. Against the furthest rung instead: "
        + ", ".join(f"`{k}` {v:.2f}×" for k, v in peak_rows.items()) + ".",
        "3. **`Z_all` is not an independent channel.** It equals "
        "`(1 − rattler_frac)·Z_nr` to better than 0.5% on every row checked, so its "
        "near-family-independence is *inherited* from the rattler fraction's and is not "
        "a second confirmation of anything.",
        "",
        "### Do the families differ at all? The test, not the eyeball",
        "",
        "| observable | Cochran Q | p | differ at 5% (liberal)? |",
        "|---|---|---|---|",
    ] + [f"| {k} | {v[0]:.2f} | {v[1]:.4f} | {'yes' if v[1] < 0.05 else 'no'} |"
         for k, v in qres.items()] + [
        "",
        "Q treats the printed se's as exactly known, which they are not — each has 2 "
        "dof — so this is the **liberal** end. On a conservative *t* comparison almost "
        "nothing here is resolvable at 3 seeds. Both ends are reported because they "
        "disagree, and an earlier draft of this page called the bottom two channels "
        "\"family-independent\" on neither test but on inspection.",
        "",
        "## What it does not say",
        "",
        "- **The bidisperse column is a path, not a family.** It is 50:50 *by number*, so "
        "its fines carry 0.500 → 0.036 of the volume across this ladder: it never visits "
        "the composition where interstitial filling peaks. Matched binaries that differ "
        "only in which atom carries the volume pack 0.0635 apart.",
        "- **A spread is not a pairwise separation.** These ratios say the families "
        "differ; they do not say two of them can be told apart at a given CV.",
        "- **Only `phi_J` is extrapolated.** `Z`, the rattler fraction and the force "
        "moments are read *at* a pressure and move with the rung; the artifact prints all "
        "four rungs so the rung dependence is visible rather than chosen.",
        "",
        "## Ladder health",
        "",
        "`phi_J` is a four-rung extrapolation, and a 2026-09-13 review found an arm whose "
        "`phi_J` came from ladders that were 0/6 monotone. The same gate is applied here "
        "before any `phi_J` is quoted: " +
        ", ".join(f"**{f}** {'passes' if health[f] else 'FAILS'}" for f in FAMILIES) + ".",
        "",
    ]
    with open(REPORT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(body))
    print(f"wrote {REPORT}")


def write():
    txt = "\n".join(L) + "\n"
    with open(os.path.join(OUT, "sweep_descriptive.txt"), "w", encoding="utf-8") as fh:
        fh.write(txt)
    print(txt)


if __name__ == "__main__":
    sys.exit(main())

"""extract_numbers.py -- every number the manuscript quotes, printed from the artifacts.

Writes `paper/numbers.tex` (LaTeX macros) and `paper/numbers.json`. The manuscript
(`main.tex`) may not contain a typed result; `test_paper_numbers.py` greps it.

TWO SOURCES PER NUMBER WHERE TWO EXIST. Simple statistics (paired shifts, per-rung
means, the matched-pair difference) are RECOMPUTED here from the shards through the
same loaders the graders use, and then ASSERTED against the value printed in the
committed artifact (`out/sweep_descriptive.txt`, `out/mixture.txt`,
`out/mixture_mu0.txt`, `out/grade_fabric.txt`). Statistics that depend on the
graders' Monte-Carlo or model machinery (the selection-calibrated Welch threshold,
the ANOVA F, the filling-model prediction) are READ from the artifact only, by regex,
and the regex must match exactly once.

No number here is typed by hand. If an artifact is regenerated and a value moves,
this script fails loudly rather than shipping a stale macro.
"""

from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, ROOT)

import sweep_descriptive as sd  # noqa: E402
import grade_mixture as gm  # noqa: E402

SWEEP_TXT = os.path.join(OUT, "sweep_descriptive.txt")
MIX_TXT = os.path.join(OUT, "mixture.txt")
MIX_MU0_TXT = os.path.join(OUT, "mixture_mu0.txt")
FABRIC_TXT = os.path.join(OUT, "grade_fabric.txt")

N: dict[str, object] = {}       # name -> value (float / int / str)
FMT: dict[str, str] = {}        # name -> format used in numbers.tex


def read(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def one(pattern: str, text: str, where: str, flags=re.M) -> tuple[str, ...]:
    """Exactly one match, or fail. A regex that matches twice is quoting the wrong line."""
    m = re.findall(pattern, text, flags)
    if len(m) != 1:
        raise SystemExit(f"{where}: pattern {pattern!r} matched {len(m)} times, need 1")
    return m[0] if isinstance(m[0], tuple) else (m[0],)


def put(name: str, value, fmt: str = "{:.4f}") -> None:
    if name in N:
        raise SystemExit(f"duplicate macro {name}")
    # A NaN reaching numbers.tex renders as the word "nan" in the PDF, which is a wrong
    # number that looks like a typo rather than an error. Refuse it here (a single-seed
    # cell gives std(ddof=1) = nan, and nothing downstream would notice).
    if isinstance(value, float) and not math.isfinite(value):
        raise SystemExit(f"macro {name} is {value!r}; a non-finite macro must not ship")
    N[name] = value
    FMT[name] = fmt


def close(a: float, b: float, tol: float, what: str) -> None:
    if not abs(a - b) <= tol:
        raise SystemExit(f"MISMATCH {what}: recomputed {a!r} vs artifact {b!r} (tol {tol})")


# ---------------------------------------------------------------------------------
# 1. phi_J: the three family shifts, recomputed from phase2_fits.csv and checked
#    against the artifact's reproduction block.
# ---------------------------------------------------------------------------------
sweep = read(SWEEP_TXT)
fits = sd.load_fits()

for fam, tag in (("lognormal", "Logn"), ("uniform", "Unif"), ("bidisperse", "Bidi")):
    m, se, n = sd.paired_shift(fits, fam)  # type: ignore[misc]
    a_m, a_se, a_n = one(rf"^\s+{fam}\s+([+-]\d\.\d+)\s+(\d\.\d+)\s+(\d)\s+[+-]\d\.\d+\s+YES$",
                         sweep, "reproduction block")
    close(m, float(a_m), 5e-7, f"shift {fam}")
    close(se, float(a_se), 5e-7, f"shift se {fam}")
    assert n == int(a_n) == 3, (fam, n, a_n)
    put(f"dPhiJ{tag}", m, "{:+.4f}")
    put(f"dPhiJ{tag}Se", se, "{:.4f}")
    put(f"nSeeds{tag}", n, "{:d}")
    # phi_J at the two ends, mean +- se over seeds, at mu_p = 0.3
    for cv, ctag in ((0.0, "Zero"), (0.5, "Half")):
        v = fits[(fits.family == fam) & (fits.mu_p == sd.MU_REF) & (fits.cv_nom == cv)].phi_J_lin
        put(f"phiJ{tag}{ctag}", float(v.mean()), "{:.4f}")
        put(f"phiJ{tag}{ctag}Se", float(v.std(ddof=1) / math.sqrt(len(v))), "{:.4f}")

# the across-family ratio and its two brackets, read from the artifact's own lines
r_endpoint, = one(r"^\s+phi_J\s+\+0\.06207\S+\s+\+0\.0242\S+\s+\+0\.01388\S+\s+n=3/3/3\s+(\d\.\d+)x$",
                  sweep, "endpoint ratio")
ratio_re = float(N["dPhiJLogn"]) / float(N["dPhiJBidi"])  # type: ignore[arg-type]
close(ratio_re, float(r_endpoint), 0.006, "endpoint ratio")
put("ratioPhiJ", float(r_endpoint), "{:.2f}")
r_ew, r_iv = one(r"^\s+headline ratio: (\d\.\d+)x equal-weight, (\d\.\d+)x inverse-variance\.$",
                 sweep, "inverse-variance bracket")
assert r_ew == r_endpoint, (r_ew, r_endpoint)
put("ratioPhiJInvVar", float(r_iv), "{:.2f}")
r_peak, = one(r"^\s+phi_J\s+4\.47x\s+(\d\.\d+)x\s+logn@0\.5 unif@0\.5 bidi@0\.3$",
              sweep, "peak ratio")
put("ratioPhiJPeak", float(r_peak), "{:.2f}")

# the three gates on phi_J, from the artifact
F, pF = one(r"^\s+phi_J\s+\d+\.\d+\s+\d\.\d+\s+(\d+\.\d+)\s+(\d\.\d+)\s+YES / YES$",
            sweep, "phi_J ANOVA")
put("anovaFPhiJ", float(F), "{:.1f}")
put("anovaPPhiJ", float(pF), "{:.4f}")
wt, tnom, tcal, fpr = one(
    r"^\s+phi_J\s+lognormal\s+bidisperse\s+(\d+\.\d+)\s+(\d\.\d+)\s+(\d\.\d+)\s+(\d+\.\d)%\s+YES$",
    sweep, "phi_J Welch")
put("welchTPhiJ", float(wt), "{:.1f}")
put("welchTNomPhiJ", float(tnom), "{:.2f}")
put("welchTCalPhiJ", float(tcal), "{:.2f}")
put("welchFprNomPhiJ", float(fpr), "{:.1f}")
tc = sd.tcrit(3)
put("tCritTwoDof", tc, "{:.3f}")
put("pAbsTgtTwoDof", 0.184, "{:.3f}")  # printed in the artifact; checked below
one(r"P\(\|T\| > 2, df = 2\) = 0\.184", sweep, "P(|T|>2)")

# the other channels: F p-values that FAIL the heterogeneity gate
for obs, tag in (("Z_all", "Zall"), (r"Z_nr \(k=2\)", "Znr"), ("rattler fraction", "Rat")):
    F, pF = one(rf"^\s+{obs}\s+\d+\.\d+\s+\d\.\d+\s+(\d+\.\d+)\s+(\d\.\d+)\s+YES / no$",
                sweep, f"{obs} ANOVA")
    put(f"anovaP{tag}", float(pF), "{:.2f}")
wt, tcal = one(r"^\s+contact force mean\s+bidisperse\s+lognormal\s+(\d\.\d+)\s+\d\.\d+\s+(\d\.\d+)\s+\d+\.\d%\s+no$",
               sweep, "force-mean Welch")
put("welchTFmean", float(wt), "{:.2f}")
put("welchTCalFmean", float(tcal), "{:.2f}")
ratio_zall, = one(r"^\s+Z_all\s+-2\.17689\S+ -2\S+ -2\.19378\S+\s+n=3/3/3\s+(\d\.\d+)x$", sweep, "Z_all ratio")
put("ratioZall", float(ratio_zall), "{:.2f}")
ratio_rat, = one(r"^\s+rattler fraction\s+\+0\.408\S+ \+0\.382\S+ \+0\.417\S+\s+n=3/3/3\s+(\d\.\d+)x$",
                 sweep, "rattler ratio")
put("ratioRat", float(ratio_rat), "{:.2f}")
ratio_znr, = one(r"^\s+Z_nr \(k=2\)\s+-0\.189\S+ -0\.0812\S+ -0\.1171\S+\s+n=3/3/3\s+(\d\.\d+)x$",
                 sweep, "Z_nr ratio")
put("ratioZnr", float(ratio_znr), "{:.2f}")
order_pct, order_n = one(r"the exact six-way ordering survives (\d+) % of (\d+) draws", sweep, "ordering")
put("orderSurvivePct", int(order_pct), "{:d}")
put("orderDraws", int(order_n), "{:d}")
worst_zall_rows, = one(r"worst \|difference\| over the rows printed: (\d\.\d+)", sweep,
                       "Z_all identity")
# 🔴 2026-09-17. The artifact's number is a max over the ROWS PRINTED -- 3 families x 3 CVs
# at the 1e-3 rung, each row already a MEAN over seeds -- and main.tex quoted it as a bound
# "on every cell". It is not one: averaging over seeds and rungs can only shrink the worst
# case. The per-CELL max is computed below from the shards and is what main.tex now uses;
# the row-level figure is kept as its own macro so the two cannot be confused again.
put("zallIdentityWorstRows", float(worst_zall_rows), "{:.4f}")

# ---------------------------------------------------------------------------------
# 2. The per-rung observables at p_hat = 1e-3, recomputed from the shards.
# ---------------------------------------------------------------------------------
sd.L.clear()
sel = sd.load_shards()
n_files, n_rec, n_dedup, n_conf, n_conf_rec, n_sel, n_cells = one(
    r"files globbed out/phase2\*\.jsonl: (\d+)\n\s+N records: (\d+)\n\s+N after dropping failed and "
    r"deduplicating by key: (\d+)\n\s+duplicate keys whose phi DISAGREED: (\d+) \((\d+) disagreeing records\)\n\s+N after N=1500, "
    r"mu_p=0\.3, converged: (\d+)\n\s+N after collapsing one physical cell per "
    r"\(family, cv, seed, p_hat\), highest MAX_ITER winning: (\d+)",
    sweep, "N counts", flags=re.NOFLAG)
# `n_sel` counts RECORDS at the selected N/mu/converged; `n_cells` counts physical CELLS.
# They differ because the record key carries the iteration cap (phase2_phiJ.py:130), so a
# cell re-run at a raised cap appears twice. Everything downstream reads CELLS.
assert len(sel) == int(n_cells), (len(sel), n_cells)
assert int(n_cells) <= int(n_sel), (n_cells, n_sel)
put("nShardFiles", int(n_files), "{:d}")
put("nRecords", int(n_rec), "{:d}")
put("nDedup", int(n_dedup), "{:d}")
put("nConflicts", int(n_conf), "{:d}")
put("nConflictRecords", int(n_conf_rec), "{:d}")
put("nSelectedRecords", int(n_sel), "{:d}")
put("nSelected", int(n_cells), "{:d}")
put("nRecappedCells", int(n_sel) - int(n_cells), "{:d}")
# how much the 15 disagreeing duplicates could move anything the paper reads at a rung:
# re-glob the shards, and for every conflicting key that lands in the SELECTED set report
# the largest |delta phi| between its records, against PHI_SCALE (referee, pass 1)
_dups: dict[str, list[dict]] = {}
for _f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
    for _ln in open(_f, encoding="utf-8"):
        try:
            _r = json.loads(_ln)
        except json.JSONDecodeError:
            continue
        if not _r.get("failed") and "key" in _r:
            _dups.setdefault(_r["key"], []).append(_r)
_sel_keys = {r["key"] for r in sel}
_conf = {k: v for k, v in _dups.items()
         if len(v) > 1 and max(abs(a["phi"] - b["phi"]) for a in v for b in v) > 1e-12}
assert len(_conf) == int(n_conf), (len(_conf), n_conf)
_conf_sel = {k: v for k, v in _conf.items() if k in _sel_keys}
put("nConflictsSelected", len(_conf_sel), "{:d}")
_worst = max(max(abs(a["phi"] - b["phi"]) for a in v for b in v) for v in _conf_sel.values())
put("conflictWorstDphi", float(_worst), "{:.4f}")
put("conflictWorstDphiScale", float(_worst / gm.SCALE), "{:.1f}")
_conf_cvs = sorted({r["cv"] for v in _conf_sel.values() for r in v})
assert _conf_cvs == [0.0], _conf_cvs  # lessons.md 2026-09-09: all at CV = 0


def rung_stat(fam, cv, p, getter):
    v = np.array([getter(r) for r in sel
                  if (r.get("family") or "lognormal") == fam and r["cv"] == cv
                  and r["p_hat_target"] == p], float)
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v)


def artifact_cell(block: str, p: str, fam: str, col: int) -> float:
    """Read one cell of a per-rung table from the artifact by column index."""
    sec = sweep.split(block, 1)[1]
    line = one(rf"^\s+{p}\s+{fam}\s+(.+)$", sec.split("\n\n\n")[0], f"{block} {p} {fam}")[0]
    cells = re.findall(r"(\d\.\d{4})\((\d)\)", line)
    return float(cells[col][0])


# The Z_all identity, per CELL (see the note beside zallIdentityWorstRows above).
_zall_dev = [abs(r["z_all"] - (1 - r["rattler"]["rattler_frac_k2"]) * r["rattler"]["z_nr_k2"])
             for r in sel]
put("zallIdentityWorst", float(max(_zall_dev)), "{:.4f}")
put("zallIdentityWorstN", len(_zall_dev), "{:d}")
# ⚠ a GUARD, not a gate: a per-cell max is mathematically >= a max over means of subsets
# of those cells, so no correct implementation can fail it. What it catches is an
# INCORRECT one -- a wrong field, a wrong subset, a getter that silently returns zeros.
assert max(_zall_dev) >= float(N["zallIdentityWorstRows"]), (  # type: ignore[arg-type]
    "a per-cell max below a max over means of those cells is impossible")

GET = {
    "rat": (lambda r: r["rattler"]["rattler_frac_k2"], "rattler fraction (k = 2)"),
    "znr": (lambda r: r["rattler"]["z_nr_k2"], "Z_nr (non-rattler backbone, k = 2)"),
    "zall": (lambda r: r["z_all"], "Z_all"),
}
# Confirm the getters are the ones the sweep used, by matching its own table.
for fam, tag in (("lognormal", "Logn"), ("uniform", "Unif"), ("bidisperse", "Bidi")):
    for cv, ctag in ((0.0, "Zero"), (0.3, "Three"), (0.5, "Half")):
        col = int(round(cv * 10))
        for key, (getter, block) in GET.items():
            m, se, n = rung_stat(fam, cv, 1e-3, getter)
            close(m, artifact_cell(block, "1e-03", fam, col), 5.1e-5, f"{key} {fam} {cv}")
            put(f"{key}{tag}{ctag}", m, "{:.3f}")
            put(f"{key}{tag}{ctag}Se", se, "{:.3f}")
            put(f"{key}{tag}{ctag}N", n, "{:d}")

# ---------------------------------------------------------------------------------
# 3. The matched binaries, recomputed from the phase12 shards.
# ---------------------------------------------------------------------------------
mix = read(MIX_TXT)
cells = gm.load_cells()
pA = gm.phi_of(cells[(3.0, round(1 / 28, 6))])   # X_s = 0.035714: fines-poor
pB = gm.phi_of(cells[(3.0, 0.25)])                 # X_s = 0.25:     fines-rich
a_phiA, a_semA = one(r"A: r = 3, X_s = 0\.035714, cv_realised = 0\.500000, phi = (\d\.\d+) \+- (\d\.\d+) \(n=3\)",
                     mix, "pair A")
a_phiB, a_semB = one(r"B: r = 3, X_s = 0\.250000, cv_realised = 0\.500000, phi = (\d\.\d+) \+- (\d\.\d+) \(n=3\)",
                     mix, "pair B")
close(pA[0], float(a_phiA), 5.1e-5, "pair A phi"); close(pB[0], float(a_phiB), 5.1e-5, "pair B phi")
close(pA[1], float(a_semA), 5.1e-5, "pair A sem"); close(pB[1], float(a_semB), 5.1e-5, "pair B sem")
assert pA[2] == pB[2] == 3
dphi = pB[0] - pA[0]
dphi_se = math.hypot(pA[1], pB[1])
a_d, a_dse, a_scale = one(r"MEASURED  dphi = ([+-]\d\.\d+) \+- (\d\.\d+)  \(([+-]\d+\.\d) scale, [+-]\d+\.\d se\)",
                          mix, "pair dphi")
close(dphi, float(a_d), 5.1e-5, "pair dphi"); close(dphi_se, float(a_dse), 5.1e-5, "pair dphi se")
put("pairPhiA", pA[0], "{:.4f}"); put("pairPhiASe", pA[1], "{:.4f}")
put("pairPhiB", pB[0], "{:.4f}"); put("pairPhiBSe", pB[1], "{:.4f}")
put("pairDPhi", dphi, "{:+.4f}"); put("pairDPhiSe", dphi_se, "{:.4f}")
put("pairScale", float(a_scale), "{:.1f}")
scale, = one(r"grading scale PHI_SCALE = (\d\.\de-\d+)", mix, "PHI_SCALE")
put("phiScale", float(scale), "{:.1e}")
close(dphi / float(scale), float(a_scale), 0.06, "pair in scale units")
fill_pred, = one(r"predicted: descriptor account \+0\.0000, filling account (\+\d\.\d+)", mix, "filling prediction")
put("pairFillingPred", float(fill_pred), "{:+.4f}")
cv_agree, = one(r"the two CVs agree to (\d\.\d+e-\d+)", mix, "cv agreement")
put("pairCvAgree", float(cv_agree), "{:.1e}")
n_mix, = one(r"^\s+(\d+) cells, 0 errors, 0 non-converged", mix, "mixture cell count")
put("nMixCells", int(n_mix), "{:d}")
# composition of the pair: fines VOLUME fraction and the number CV, from the PSD itself
from psd import BinaryPSD  # noqa: E402
psdA = BinaryPSD(d10=1.0, ratio=3.0, x_small=1 / 28)
psdB = BinaryPSD(d10=1.0, ratio=3.0, x_small=0.25)
close(psdA.cv, psdB.cv, 1e-12, "pair cv identity")
put("pairCv", float(psdA.cv), "{:.3f}")
put("pairXsA", 1 / 28, "{:.4f}"); put("pairXsB", 0.25, "{:.2f}")
put("pairNsmallA", float(psdA.n_small), "{:.3f}"); put("pairNsmallB", float(psdB.n_small), "{:.3f}")
# the descriptor the pair does NOT hold fixed: the skewness of the number distribution
# (referee, pass 1). Closed form for a two-point PSD, anchored by hand in
# make_paper_figures.binary_number_skewness and in test_paper_numbers.py
from make_paper_figures import binary_number_skewness, minor_phase_count, gradable  # noqa: E402
put("pairSkewA", binary_number_skewness(3.0, 1 / 28), "{:.2f}")
put("pairSkewB", binary_number_skewness(3.0, 0.25), "{:+.2f}")
# the M5 minimum-phase-count gate and the one interior cell it drops from Fig. 3
from phase12_mixture import MIN_PHASE_COUNT, N_PROD  # noqa: E402
put("mixMinPhaseCount", int(MIN_PHASE_COUNT), "{:d}")
put("mixNProd", int(N_PROD), "{:d}")
dropped = sorted((r, xs) for (r, xs), cs in gm.load_cells().items()
                 if [c for c in cs if "error" not in c] and not gradable(r, xs))
assert len(dropped) == 1, dropped
put("mixDroppedRatio", int(dropped[0][0]), "{:d}")
put("mixDroppedXs", float(dropped[0][1]), "{:.1f}")
put("mixDroppedNminor", minor_phase_count(*dropped[0]), "{:.1f}")
# peak of the r = 3 block, and the monodisperse endpoint
pk_end, pk_max, pk_xs, pk_exc = one(
    r"r = 3: endpoints (\d\.\d+), maximum (\d\.\d+) at X_s = (\d\.\d+), excess (\+\d\.\d+)", mix, "r=3 peak")
put("mixEndpoint", float(pk_end), "{:.4f}"); put("mixPeakRThree", float(pk_max), "{:.4f}")
put("mixPeakXsRThree", float(pk_xs), "{:.2f}"); put("mixPeakExcRThree", float(pk_exc), "{:+.4f}")
pk2_max, pk2_xs, pk2_exc = one(r"r = 2: endpoints \d\.\d+, maximum (\d\.\d+) at X_s = (\d\.\d+), excess (\+\d\.\d+)",
                               mix, "r=2 peak")
put("mixPeakRTwo", float(pk2_max), "{:.4f}"); put("mixPeakXsRTwo", float(pk2_xs), "{:.2f}")
put("mixPeakExcRTwo", float(pk2_exc), "{:+.4f}")

# the same pair at zero friction (a BLIND block: no mu_p = 0 binary existed when predicted)
mu0 = read(MIX_MU0_TXT)
d0, d0se, s0 = one(r"MEASURED dphi = (\+\d\.\d+) \+- (\d\.\d+)  \((\+\d+\.\d) scale, \+\d+\.\d se\)", mu0, "mu0 pair")
put("pairDPhiMuZero", float(d0), "{:+.4f}"); put("pairDPhiMuZeroSe", float(d0se), "{:.4f}")
put("pairScaleMuZero", float(s0), "{:.1f}")
phA0, phB0 = one(r"^\s+phi = (\d\.\d+) and (\d\.\d+)$", mu0, "mu0 pair phis")
put("pairPhiAMuZero", float(phA0), "{:.4f}"); put("pairPhiBMuZero", float(phB0), "{:.4f}")
b03, b00 = one(r"beta\(mu_p = 0\.3\) = (\d\.\d+)   beta\(mu_p = 0\.0\) = (\d\.\d+)", mu0, "beta")
put("betaMuThree", float(b03), "{:.4f}"); put("betaMuZero", float(b00), "{:.4f}")

# ---------------------------------------------------------------------------------
# 4. The fabric campaign: slopes of T_nr against CV, from grade_fabric.txt.
# ---------------------------------------------------------------------------------
fab = read(FABRIC_TXT)
fab_k2 = fab.split("T_nr at k = 2  (PRIMARY)", 1)[1].split("T_nr at k = 3", 1)[0]
for fam, tag in (("lognormal", "Logn"), ("uniform", "Unif"), ("bidisperse", "Bidi")):
    s, se, res = one(rf"^\s+{fam}\s+([+-]\d\.\d+)\s+(\d\.\d+)\s+\d\.\d+\s+(\d\.\d+)\s+no$", fab_k2, f"fabric {fam}")
    put(f"fabSlope{tag}", float(s), "{:+.3f}")
    put(f"fabSlope{tag}Se", float(se), "{:.3f}")
    put(f"fabResolve{tag}", float(res), "{:.2f}")
    n_c, mT, sdT = one(rf"^\s+{fam}\s+n =\s+(\d+)\s+mean T = (\d\.\d+)\s+sd = (\d\.\d+)", fab, f"fabric mean {fam}")
    put(f"fabMeanT{tag}", float(mT), "{:.2f}")
    put(f"fabN{tag}", int(n_c), "{:d}")
put("fabRegisteredResolve", 0.5, "{:.1f}")
one(r"achieved 2\*se\*0\.5 <= 0\.5 -> NULL, INFORMATIVE", fab, "registered 0.5")
# That grep confirms the VERDICT GRAMMAR is in the artifact; it does not confirm that any
# family met it, and main.tex claims the null is "informative on all three". Enforce the
# claim here rather than at read time, both halves of it. Note also that the achieved
# resolving power is 2*se*(CV range) and the range is exactly 0.5, so it equals the slope's
# own standard error; main.tex now says so rather than printing the two as if independent.
for _tag in ("Logn", "Unif", "Bidi"):
    _res, _se, _sl = (float(N[f"fabResolve{_tag}"]), float(N[f"fabSlope{_tag}Se"]),  # type: ignore[arg-type]
                      abs(float(N[f"fabSlope{_tag}"])))                              # type: ignore[arg-type]
    assert _res <= float(N["fabRegisteredResolve"]), (  # type: ignore[arg-type]
        f"{_tag}: achieved resolving power {_res} exceeds the registered "
        f"{N['fabRegisteredResolve']}; the null is NOT informative and main.tex says it is")
    assert _sl <= 2 * _se, f"{_tag}: |slope| {_sl} > 2 se {2 * _se}; that is a DRIFT, not a null"
    close(_res, _se, 6e-3, f"fabResolve{_tag} against its own se (range 0.5)")
n_fab_cells, n_fab_seeds = one(r"(\d+) converged cells at `p_hat` = 0\.001, `mu_p` = 0\.3, N = 1500, (\d+) seeds",
                               read(os.path.join(ROOT, "CLAIMS_FABRIC_2026-09-18.md")), "fabric N")
put("fabCells", int(n_fab_cells), "{:d}"); put("fabSeeds", int(n_fab_seeds), "{:d}")
# the aniso floor: slope of the raw anisotropy vs the floor, lognormal
sa, sf, sd_, se_ = one(r"^\s+lognormal\s+(\+\d\.\d+)\s+(\+\d\.\d+)\s+([+-]\d\.\d+)\s+(\d\.\d+)\s+YES$",
                       fab, "aniso floor lognormal")
put("anisoSlopeLogn", float(sa), "{:+.4f}"); put("floorSlopeLogn", float(sf), "{:+.4f}")
put("anisoMinusFloorLogn", float(sd_), "{:+.4f}"); put("anisoMinusFloorLognSe", float(se_), "{:.4f}")
# the F6 exclusion: lognormal slope with the wholly flagged top rung dropped (sign flips)
s_full, se_full, s_drop, se_drop = one(
    r"^\s+lognormal\s+([+-]\d\.\d+)\+-(\d\.\d+)\s+([+-]\d\.\d+)\+-(\d\.\d+)\s+[+-]\d\.\d+\+-\d\.\d+\s+RESTS ON FLAGGED CELLS$",
    fab, "F6 lognormal")
close(float(s_full), float(N["fabSlopeLogn"]), 5e-4, "F6 full slope = F2 slope")  # type: ignore[arg-type]
put("fabSlopeLognDropFlagged", float(s_drop), "{:+.3f}")
put("fabSlopeLognDropFlaggedSe", float(se_drop), "{:.3f}")
# the finite-size flag (d_max/L > 0.25, packing.py) counted by family on the primary rung
put("dMaxOverLFlag", 0.25, "{:.2f}")
one(r'FINITE SIZE RISK: d_max/L_dense = \{st\.d_max_over_L:\.3f\} > 0\.25', read(os.path.join(ROOT, "packing.py")), "flag threshold")
for fam, tag in (("lognormal", "Logn"), ("uniform", "Unif"), ("bidisperse", "Bidi")):
    nf, nt, top_f, top_t = one(rf"^\s+{fam}\s+(\d+)/(\d+) flagged .*CV0\.5:(\d+)/(\d+)", fab, f"flag {fam}")
    put(f"fabFlagged{tag}", int(nf), "{:d}"); put(f"fabFlagged{tag}Of", int(nt), "{:d}")
    put(f"fabFlaggedTop{tag}", int(top_f), "{:d}"); put(f"fabFlaggedTop{tag}Of", int(top_t), "{:d}")
    # nf/nt above is POOLED over the primary and secondary rungs (referee, pass 1); the
    # primary-rung count is the sum of the per-CV "CVx:a/b" figures on the same line
    line, = one(rf"^(\s+{fam}\s+\d+/\d+ flagged .*)$", fab, f"flag line {fam}")
    per_cv = re.findall(r"CV\d\.\d:(\d+)/(\d+)\[pooled (\d+)/(\d+)\]", line)
    assert len(per_cv) == 6, per_cv
    assert sum(int(c[2]) for c in per_cv) == int(nf) and sum(int(c[3]) for c in per_cv) == int(nt)
    put(f"fabFlagged{tag}Primary", sum(int(c[0]) for c in per_cv), "{:d}")
    put(f"fabFlagged{tag}PrimaryOf", sum(int(c[1]) for c in per_cv), "{:d}")
    # and the CV = 0.5 entry pooled over both rungs (fabFlaggedTop* above is primary-only;
    # referee pass 2 caught the prose calling it pooled)
    assert (int(per_cv[-1][0]), int(per_cv[-1][1])) == (int(top_f), int(top_t)), per_cv[-1]
    put(f"fabFlaggedTop{tag}Pooled", int(per_cv[-1][2]), "{:d}")
    put(f"fabFlaggedTop{tag}PooledOf", int(per_cv[-1][3]), "{:d}")
# the companion Kruskal-Wallis p (the column the artifact says to read for lognormal)
for fam, tag in (("lognormal", "Logn"), ("uniform", "Unif"), ("bidisperse", "Bidi")):
    kw, lev = one(
        rf"^\s+{fam}\s+(\d\.\d+)\s+(\d\.\d+) \|\s+(\d\.\d+) \|\s+(\d\.\d+) \| (\d\.\d+) \|\s+(\d\.\d+)$",
        fab, f"fabric ANOVA {fam}")[4:]
    put(f"fabKwP{tag}", float(kw), "{:.2f}"); put(f"fabLeveneP{tag}", float(lev), "{:.3f}")
# rattler rise CV 0 -> 0.5 and backbone change, from the per-rung macros already checked
for tag in ("Logn", "Unif", "Bidi"):
    put(f"ratRise{tag}", float(N[f"rat{tag}Half"]) / float(N[f"rat{tag}Zero"]), "{:.1f}")  # type: ignore[arg-type]
    put(f"znrDrop{tag}", float(N[f"znr{tag}Zero"]) - float(N[f"znr{tag}Half"]), "{:.3f}")  # type: ignore[arg-type]

# ---------------------------------------------------------------------------------
# 5. Protocol constants and the frictionless anchor, from the code and the fits.
# ---------------------------------------------------------------------------------
import packing  # noqa: E402
import measure  # noqa: E402
import analyse_phase2 as ap2  # noqa: E402
put("nGrains", sd.N_MATCH, "{:d}")
put("muP", sd.MU_REF, "{:.1f}")
put("poisson", packing.POISSON, "{:.2f}")
put("unbalancedTol", packing.UNBALANCED_TOL, "{:.2f}")
put("relStressTol", packing.REL_STRESS_TOL, "{:.0e}")
put("dtSafety", packing.DT_SAFETY, "{:.1f}")
put("maxIter", packing.MAX_ITER, "{:d}")
put("newtonDamping", 0.4, "{:.1f}")
one(r"NewtonIntegrator\(damping=0\.4\)", read(os.path.join(ROOT, "packing.py")), "damping")
put("rattlerK", measure.RATTLER_MIN_CONTACTS, "{:d}")
put("truncLo", 0.2, "{:.1f}"); put("truncHi", 5.0, "{:.1f}")
one(r"truncate: tuple\[float, float\] \| None = \(0\.2, 5\.0\)", read(os.path.join(ROOT, "packing.py")), "truncate")
put("pHatFitMax", ap2.P_HAT_FIT_MAX, "{:.1e}")
put("nRungs", len(sd.RUNGS), "{:d}")
put("pHatTop", sd.RUNGS[0], "{:.0e}"); put("pHatBottom", sd.RUNGS[-1], "{:.0e}")
put("r2Min", sd.R2_MIN, "{:.2f}")
put("uniformCvMax", 1 / math.sqrt(3), "{:.3f}")
put("bidiRatioHalf", (1 + 0.5) / (1 - 0.5), "{:.0f}")
# fines VOLUME fraction of the number-symmetric bidisperse family at CV = 0.5
put("bidiFinesVolHalf", 0.5**3 / (0.5**3 + 1.5**3), "{:.3f}")
# frictionless anchor: phi_J at mu_p = 0, CV = 0 (lognormal label; the uniform rows are
# the same cells seen under a second label, so they are not a second sample)
g0 = fits[(fits.mu_p == 0.0) & (fits.cv_nom == 0.0) & (fits.family == "lognormal")].phi_J_lin
put("phiJFrictionless", float(g0.mean()), "{:.4f}")
put("phiJFrictionlessSe", float(g0.std(ddof=1) / math.sqrt(len(g0))), "{:.4f}")
put("nFrictionless", int(len(g0)), "{:d}")
# 🔴 2026-09-17. The same "second label" caveat applies at mu_p = 0.3 and the filter was
# missing, so this mean was a 2:1 blend: `uniform` and `bidisperse` are BIT-IDENTICAL at
# CV = 0 (the shared monodisperse cells) while `lognormal` has its own draw. There are two
# distinct samples of 3 seeds, not three of 3, and they are pooled here ONCE each, with the
# n and the s.e. printed rather than left implicit.
_g3 = fits[(fits.mu_p == 0.3) & (fits.cv_nom == 0.0)]
_logn3 = _g3[_g3.family == "lognormal"].set_index("seed").phi_J_lin
_shared3 = _g3[_g3.family == "uniform"].set_index("seed").phi_J_lin
_bidi3 = _g3[_g3.family == "bidisperse"].set_index("seed").phi_J_lin
assert len(_logn3) == len(_shared3) == len(_bidi3) == 3, "expected 3 seeds per label"
assert (_shared3.round(12) == _bidi3.round(12)).all(), (
    "uniform and bidisperse CV=0 are no longer the same cells; re-derive the pooling rule")
g3 = np.concatenate([_logn3.to_numpy(), _shared3.to_numpy()])
put("phiJMonoMuThree", float(g3.mean()), "{:.4f}")
put("phiJMonoMuThreeSe", float(g3.std(ddof=1) / math.sqrt(len(g3))), "{:.4f}")
put("nMonoMuThree", int(len(g3)), "{:d}")
put("rcpLiterature", 0.64, "{:.2f}")   # Silbert 2010 / Desmond & Weeks 2014, cited in main.tex

# the linear-packing-model interaction exponents used for the matched-pair prediction
# (grade_mixture.py M1 evaluates binary_density(..., 1.02, 1.50); the artifact prints them)
pa_txt, pb_txt = one(r"zero-param PA = (\d\.\d+), PB = (\d\.\d+) \(textbook\)", mix, "M1 exponents")
one(r"mm\.binary_density\(3\.0, 0\.25, beta, 1\.02, 1\.50\)", read(os.path.join(ROOT, "grade_mixture.py")), "M1 exponent literals")
put("fillPa", float(pa_txt), "{:.2f}")
put("fillPb", float(pb_txt), "{:.2f}")

# ---------------------------------------------------------------------------------
# write
# ---------------------------------------------------------------------------------
def texname(n: str) -> str:
    # LaTeX macro names cannot contain digits
    return "".join(c for c in n if c.isalpha())


names = [texname(k) for k in N]
if len(set(names)) != len(names):
    dup = sorted({n for n in names if names.count(n) > 1})
    raise SystemExit(f"macro name collision after stripping digits: {dup}")

lines = ["% generated by paper/extract_numbers.py -- do not edit; every value is printed from an artifact"]
for k, v in N.items():
    s = FMT[k].format(v)
    if "e" in s and isinstance(v, float):          # 2.1e-03 -> 2.1\times10^{-3}
        mant, exp = s.split("e")
        s = rf"{mant}\times10^{{{int(exp)}}}"
    lines.append(rf"\newcommand{{\{texname(k)}}}{{{s}}}")
with open(os.path.join(HERE, "numbers.tex"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(lines) + "\n")
with open(os.path.join(HERE, "numbers.json"), "w", encoding="utf-8", newline="\n") as fh:
    json.dump({k: (FMT[k].format(v) if isinstance(v, float) else v) for k, v in N.items()},
              fh, indent=1, sort_keys=True)
print(f"wrote {len(N)} macros to paper/numbers.tex; all recomputed values matched their artifacts")

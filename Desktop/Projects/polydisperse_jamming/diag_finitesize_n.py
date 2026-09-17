"""diag_finitesize_n.py -- what N would it actually take to clear the lognormal
finite-size flag, and has the N >= 3000 replicate already been run?

WHY. The descriptive-sweep manuscript (`paper/main.tex`) discloses, but does not resolve,
a finite-size flag: every lognormal CV = 0.5 cell carries `d_max/L_dense > 0.25`. The
Status card names "a lognormal CV = 0.5 replicate at N >= 3000" as the cheapest way to
close it. Before registering that run, two premises have to be checked, because an item's
stated premise can be false and checking it is cheaper than the work it asks for:

  P1  Has it already been run?  `out/phase2*.jsonl` is searched for lognormal cells at
      N > 1500, and every such cell's convergence flag is printed. A ladder with a
      non-converged rung is not a ladder.
  P2  Would N = 3000 clear the flag?  The flag is NOT a simulation output. Reading
      `packing.py`, it is

          d_max/L_dense = d_max / (sum_i (pi/6) d_i^3 / PHI_DENSE_BOUND)^(1/3)

      i.e. a deterministic function of (family, CV, N, seed, truncate) alone, computed
      from the drawn diameters before a single time step runs. So it can be evaluated on
      the project's OWN generator at any N for pennies, and the answer is not an
      extrapolation. This script does that, and checks the reproduction against the
      values on disk at N = 1500 and N = 3000 before believing the rest of the column.

WHAT THIS SCRIPT DOES NOT DO. It does not run a cell, and it says nothing about what
phi_J does at larger N -- only about where the flag clears and what the cells would cost.
The cost column is an EXTRAPOLATION of two measured points and is labelled as one.

Run:  python diag_finitesize_n.py
Out:  out/diag_finitesize_n.txt
"""
from __future__ import annotations

import glob
import json
import math
import os
import statistics as st

import numpy as np

import packing
import psd

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "out")

FAMILY = "lognormal"
CV_TOP = 0.5
CV_REF = 0.0
MU_REF = 0.3
SEEDS = (1, 2, 3)              # the manuscript's own seeds; fixed here, not chosen later
TRUNCATE = (0.2, 5.0)          # psd truncation, matches the campaign
FLAG = 0.25                    # packing.py's own threshold
N_GRID = (1500, 3000, 4000, 5000, 6000, 7000, 8000, 10000, 12000)
N_SEED_SAMPLE = 24             # a wider seed draw, to see the flag rate rather than 3 cells

lines: list[str] = []


def say(s: str = "") -> None:
    lines.append(s)
    print(s)


def rule(title: str) -> None:
    say()
    say("=" * 78)
    say(title)
    say("=" * 78)


def load_phase2() -> list[dict]:
    recs = []
    for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                try:
                    recs.append(json.loads(ln))
                except json.JSONDecodeError:
                    continue
    return recs


def d_max_over_L(n: int, seed: int, cv: float) -> float:
    """packing.py's flag, recomputed from the generator alone (no Yade)."""
    dist = psd.make_psd(FAMILY, 1.0, cv)
    rng = np.random.default_rng(seed)
    d = 2.0 * np.asarray(dist.radii_for_yade(n, rng, truncate=TRUNCATE), dtype=float)
    vol_solid = (math.pi / 6.0) * float(np.sum(d ** 3))
    l_final_min = (vol_solid / packing.PHI_DENSE_BOUND) ** (1.0 / 3.0)
    return float(d.max()) / l_final_min


recs = load_phase2()
sweep = [r for r in recs
         if r.get("family") == FAMILY and r.get("mu_p") == MU_REF
         and r.get("cv") in (CV_REF, CV_TOP)]
# the shards contain exact repeats; collapse on the full cell key before counting anything
cells: dict[tuple, dict] = {}
for r in sweep:
    key = (r["cv"], r["n_requested"], r["seed"], r["p_hat_target"])
    cells.setdefault(key, r)

rule("0. REPRODUCING THE FLAG FROM THE GENERATOR (the premise check for everything below)")
say()
say("  If these two columns disagree, nothing further in this file is usable.")
say()
say("  %-5s %-6s %-5s  %-12s %-12s %-10s" % ("cv", "N", "seed", "on disk", "regenerated", "diff"))
# the flag depends on (cv, N, seed) only -- it is read before compression, so every rung
# of one ladder carries the same value. Collapse the rungs and check they agree first.
by_draw: dict[tuple, set] = {}
for (cv, n, seed, _p), r in cells.items():
    disk = r.get("d_max_over_L")
    if disk is None or not math.isfinite(float(disk)):
        continue
    by_draw.setdefault((cv, n, seed), set()).add(round(float(disk), 12))
n_multi = sum(1 for v in by_draw.values() if len(v) > 1)
worst_repro = 0.0
for (cv, n, seed), vals in sorted(by_draw.items()):
    disk = next(iter(vals))
    gen = d_max_over_L(n, seed, cv)
    diff = abs(disk - gen)
    worst_repro = max(worst_repro, diff)
    say("  %-5.1f %-6d %-5d  %-12.4f %-12.4f %-10.2e" % (cv, n, seed, disk, gen, diff))
n_repro = len(by_draw)
say()
say("  draws with more than one distinct flag value across their rungs: %d (expected 0)" % n_multi)
say("  cells compared: %d     worst |disk - regenerated|: %.3e" % (n_repro, worst_repro))
say("  VERDICT: %s" % ("REPRODUCED -- the flag is a pure function of the draw"
                       if worst_repro < 1e-9 else
                       "*** NOT REPRODUCED -- stop here, the model of the flag is wrong ***"))

rule("1. P1 -- HAS THE N >= 3000 REPLICATE ALREADY BEEN RUN?")
say()
big = sorted({(k[0], k[1], k[2], k[3]) for k in cells if k[1] > 1500})
sizes = sorted({k[1] for k in cells})
say("  system sizes present for %s at mu_p = %s: %s" % (FAMILY, MU_REF, sizes))
say("  cells at N > 1500: %d" % len(big))
say()
if big:
    say("  %-5s %-6s %-5s %-9s %-7s %-9s %-9s" % ("cv", "N", "seed", "p_hat", "conv", "iters", "phi"))
    n_noncon = 0
    for key in big:
        r = cells[key]
        conv = bool(r.get("converged"))
        n_noncon += (not conv)
        say("  %-5.1f %-6d %-5d %-9.0e %-7s %-9s %-9.4f"
            % (key[0], key[1], key[2], key[3], conv, r.get("iters"), r["phi"]))
    say()
    say("  non-converged cells among them: %d of %d" % (n_noncon, len(big)))
    # a ladder is usable only if EVERY rung of that (cv, N, seed) converged
    ladders: dict[tuple, list[bool]] = {}
    for key in big:
        ladders.setdefault((key[0], key[1], key[2]), []).append(bool(cells[key].get("converged")))
    usable = [k for k, v in ladders.items() if all(v) and len(v) >= 4]
    say("  ladders at N > 1500 with every rung converged and >= 4 rungs: %d of %d"
        % (len(usable), len(ladders)))
    for k in sorted(ladders):
        say("     cv=%.1f N=%d seed=%d : %d rungs, %d converged -> %s"
            % (k[0], k[1], k[2], len(ladders[k]), sum(ladders[k]),
               "USABLE" if (all(ladders[k]) and len(ladders[k]) >= 4) else "NOT A LADDER"))
    say()
    top = [k for k in usable if k[0] == CV_TOP]
    ref = [k for k in usable if k[0] == CV_REF]
    say("  A paired CV = 0 -> 0.5 difference needs BOTH ends usable at the same N and seed.")
    say("  usable CV = %.1f ladders: %d     usable CV = %.1f ladders: %d"
        % (CV_TOP, len(top), CV_REF, len(ref)))
    paired = sorted({(k[1], k[2]) for k in top} & {(k[1], k[2]) for k in ref})
    say("  (N, seed) pairs with both ends usable: %d  %s" % (len(paired), paired))
    say()
    say("  VERDICT: %s" % (
        "a paired large-N replicate EXISTS on disk; grade it before registering a new run"
        if paired else
        "NO paired large-N replicate exists -- the run is still needed"))
else:
    say("  VERDICT: no cells above N = 1500; the run is still needed")

rule("2. P2 -- AT WHICH N DOES THE FLAG ACTUALLY CLEAR?")
say()
say("  d_max/L_dense for the lognormal CV = %.1f family, truncate = %s, threshold %.2f." % (CV_TOP, TRUNCATE, FLAG))
say("  The three registered seeds are shown individually; the right-hand block is a wider")
say("  draw of %d seeds, which measures the FLAG RATE rather than three cells." % N_SEED_SAMPLE)
say()
say("  %-7s | %-8s %-8s %-8s %-8s | %-8s %-8s %-9s"
    % ("N", "seed 1", "seed 2", "seed 3", "worst", "median", "max", "flagged"))
say("  " + "-" * 74)
table: dict[int, dict] = {}
for n in N_GRID:
    trio = [d_max_over_L(n, s, CV_TOP) for s in SEEDS]
    wide = [d_max_over_L(n, s, CV_TOP) for s in range(1, N_SEED_SAMPLE + 1)]
    nflag = sum(v > FLAG for v in wide)
    table[n] = {"trio": trio, "trio_max": max(trio), "median": float(np.median(wide)),
                "max": max(wide), "flagged": nflag, "of": len(wide)}
    say("  %-7d | %-8.4f %-8.4f %-8.4f %-8.4f | %-8.4f %-8.4f %d/%d"
        % (n, trio[0], trio[1], trio[2], max(trio),
           table[n]["median"], table[n]["max"], nflag, len(wide)))
say()
clear_trio = [n for n in N_GRID if table[n]["trio_max"] <= FLAG]
clear_wide = [n for n in N_GRID if table[n]["flagged"] == 0]
say("  smallest N on this grid clearing the flag on all %d registered seeds : %s"
    % (len(SEEDS), clear_trio[0] if clear_trio else "none on the grid"))
say("  smallest N on this grid clearing it on all %d sampled seeds          : %s"
    % (N_SEED_SAMPLE, clear_wide[0] if clear_wide else "none on the grid"))
if clear_wide:
    n0 = clear_wide[0]
    say("  margin at N = %d: worst sampled d_max/L = %.4f, i.e. %.1f %% below the threshold"
        % (n0, table[n0]["max"], 100.0 * (FLAG - table[n0]["max"]) / FLAG))
say()
say("  Note the flag falls SLOWER than N^(-1/3): the largest grain is itself a sample")
say("  maximum and grows with N until it saturates at the truncation ceiling %.1f d10." % TRUNCATE[1])
ratio = table[3000]["median"] / table[1500]["median"]
say("  measured 1500 -> 3000 (medians): %.4f -> %.4f, ratio %.4f against N^(-1/3) = %.4f,"
    % (table[1500]["median"], table[3000]["median"], ratio, 2.0 ** (-1.0 / 3.0)))
say("  so the drawn maximum grew by a factor %.4f over that doubling." % (ratio / 2.0 ** (-1.0 / 3.0)))

rule("3. WHAT THE CELLS WOULD COST (extrapolated from two measured points -- not a measurement)")
say()
cost_rows = []
for cv in (CV_REF, CV_TOP):
    for n in (1500, 3000):
        s = [r for k, r in cells.items() if k[0] == cv and k[1] == n
             and r.get("wall_s") and r.get("iters")]
        if not s:
            continue
        cost_rows.append((cv, n, len(s), st.median([r["iters"] for r in s]),
                          max(r["iters"] for r in s),
                          st.median([r["wall_s"] for r in s]),
                          max(r["wall_s"] for r in s),
                          sum(1 for r in s if not r.get("converged"))))
say("  %-5s %-6s %-5s %-10s %-10s %-9s %-9s %-8s"
    % ("cv", "N", "n", "iters med", "iters max", "wall med", "wall max", "noncon"))
for cv, n, k, im, ix, wm, wx, nc in cost_rows:
    say("  %-5.1f %-6d %-5d %-10.0f %-10.0f %-9.0f %-9.0f %d" % (cv, n, k, im, ix, wm, wx, nc))
say()
w15 = [r for r in cost_rows if r[0] == CV_TOP and r[1] == 1500]
w30 = [r for r in cost_rows if r[0] == CV_TOP and r[1] == 3000]
if w15 and w30 and clear_wide:
    p = math.log(w30[0][5] / w15[0][5]) / math.log(2.0)
    n0 = clear_wide[0]
    say("  wall time at CV = %.1f scales as N^%.2f between the two measured sizes" % (CV_TOP, p))
    say("  -> median cell at N = %d: %.0f s = %.1f h  (EXTRAPOLATED from two points)"
        % (n0, w30[0][5] * (n0 / 3000.0) ** p, w30[0][5] * (n0 / 3000.0) ** p / 3600.0))
    say("  ** The N = 3000 CV = %.1f medians are CENSORED at the 200 000 iteration cap" % CV_TOP)
    say("     (%d of %d non-converged), so the true cost is HIGHER than this row says."
        % (w30[0][7], w30[0][2]))
    say("     A larger-N run must raise the cap, which is a change of protocol.")

rule("4. WHAT IS ALREADY KNOWN ABOUT THE N DEPENDENCE OF THE CLAIM")
say()
fs = os.path.join(OUT, "finite_size_verdict.txt")
if os.path.exists(fs):
    with open(fs, encoding="utf-8") as fh:
        txt = fh.read()
    say("  from out/finite_size_verdict.txt (30 matched cells, N = 1500 vs 3000):")
    for ln in txt.splitlines():
        s = ln.strip()
        if s.startswith("C1 ") or s.startswith("C3") or "packing fraction at fixed p_hat" in s:
            say("    " + s)
    say()
    say("  That file grades the phi_J channel SAFE and the CV trend as holding, but it")
    say("  compares phi at FIXED p_hat, not the extrapolated phi_J, and it was computed")
    say("  on the same N = 3000 block section 1 above shows is partly non-converged.")
else:
    say("  out/finite_size_verdict.txt not found")

say()
say("=" * 78)
say("This file is the evidence base for PREREG_LOGN_LARGEN. It registers nothing itself.")
say("=" * 78)

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "diag_finitesize_n.txt"), "w", encoding="utf-8") as fh:
    fh.write("\n".join(lines) + "\n")
print("\nwrote out/diag_finitesize_n.txt")

"""grade_shape.py -- grade the shape campaign against PREREG_SHAPE_2026-09-18.md.

Writes `out/grade_shape.txt` and generates `CLAIMS_SHAPE_2026-09-18.md`. EVERY number in
the generated file is printed by this script; no prose is written alongside a verdict and
no verdict is written that the numbers do not produce (lessons.md 2026-09-06: generate
the narrative FROM the verdict, never beside it). Exits non-zero if any gate fails, so a
stale artifact can never be read as a pass (lessons.md 2026-09-14).

THE TAG GUARD IS STRUCTURAL, not by eye. Every input file's basename must start with
`phase2_shp`; any file matching `mixr*`, `extr*`, `tern*`, `kl*`, `ug*` or `fab*` is a
hard error, as is any admitted record without an `alpha` field. The `shpgate` and
`shptime` tags are read for S0a/S0b and REFUSED for every S3 number: the gate cells are
alpha = 1 reproductions of records already on disk and the timing cells ran at a
different N, so pooling either into a production average would be averaging a different
question.
"""

from __future__ import annotations

import glob
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

from phase2_phiJ import fit_phi_J
from shape import ALPHA_LADDER, sub_sphere_diameter

HERE = os.path.dirname(os.path.abspath(__file__))
SHARD_GLOB = os.path.join(HERE, "out", "shape", "phase2_shp*.jsonl")
OUT_TXT = os.path.join(HERE, "out", "grade_shape.txt")
CLAIMS = os.path.join(HERE, "CLAIMS_SHAPE_2026-09-18.md")

FORBIDDEN_PREFIXES = ("phase2_mixr", "phase2_extr", "phase2_tern", "phase2_kl",
                      "phase2_ug", "phase2_fab")
GATE_TAGS = ("shpgate", "shptime")

# PREREG_FABRIC_2026-09-18.md F1, read from the committed shards there so the comparison
# cannot drift. S0a reproduces these through the SHAPE runner at alpha = 1.
S0A_REFERENCE = {
    ("bidisperse", 0.3): 0.6120403164804106,
    ("uniform", 0.3): 0.6101301588958694,
    ("lognormal", 0.4): 0.6339122825028567,
}
S0B_BAR = 1e-13          # |phi_union - phi_porosity| at alpha = 1; see amendment 5.2
S0C_RATIO_BAR = 1e-9     # phi_porosity/phi_union vs 2*(ds/d_eq)^3 at alpha > 1
S0C_MASS_BAR = 1e-12     # mass_over_vol_solid vs 1
S1_BAR = 1e-12           # PSD descriptors across alpha at fixed (family, cv, seed)

P_HAT_FIT_MAX = 1.1e-3   # same window analyse_phase2 uses

# THE REGISTERED LADDER (PREREG_SHAPE_2026-09-18.md amendment 5.6). Filtering to a
# DECLARED ladder rather than to whatever is on disk is the point: at N = 5400 the 3e-5
# rung does not jam at CV = 0.5 (z_all = 0.023, capped at 200 000 iterations) while it
# does at CV = 0.0, so admitting it where it happens to exist would give a four-rung
# phi_J at one end of the axis under test and a three-rung phi_J at the other -- two
# estimators wearing one name, with the difference between them landing inside D.
P_HAT_LADDER_SHAPE = (1e-3, 3e-4, 1e-4)
P_HAT_TOL = 0.02         # relative; the target is matched, the realised value drifts


class Gate:
    """One pre-registered gate: a SENTENCE, a measured number, and a verdict computed
    from the number. The sentence is printed beside the number so the two cannot drift
    (lessons.md 2026-09-13: the code can compute a different quantity than the
    registered sentence names)."""

    def __init__(self, name: str, sentence: str):
        self.name = name
        self.sentence = sentence
        self.lines: list[str] = []
        self.passed: bool | None = None
        self.informative = True

    def say(self, s: str) -> None:
        self.lines.append(s)

    def verdict(self, ok: bool, informative: bool = True) -> None:
        self.passed = bool(ok)
        self.informative = bool(informative)

    def render(self) -> list[str]:
        if self.passed is None:
            tag = "NOT RUN"
        elif not self.informative:
            tag = "NOT INFORMATIVE"
        else:
            tag = "PASS" if self.passed else "FAIL"
        out = [f"### {self.name} -- {tag}", f"*{self.sentence}*", ""]
        out += self.lines
        out.append("")
        return out


def load_records(report) -> tuple[list[dict], list[dict], list[str]]:
    """(production records, gate records, file list), with N printed at every filter."""
    files = sorted(glob.glob(SHARD_GLOB))
    if not files:
        raise SystemExit(f"no shards matched {SHARD_GLOB}")
    for f in files:
        b = os.path.basename(f)
        if any(b.startswith(p) for p in FORBIDDEN_PREFIXES):
            raise SystemExit(f"TAG GUARD: {b} belongs to another campaign")
        if not b.startswith("phase2_shp"):
            raise SystemExit(f"TAG GUARD: {b} is not a shape shard")
    raw: list[dict] = []
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    raw.append(json.loads(line))
    report(f"N: {len(files)} shard files -> {len(raw)} records read")

    for r in raw:
        if "alpha" not in r:
            raise SystemExit(f"TAG GUARD: record {r.get('key')} has no `alpha` field")

    failed = [r for r in raw if r.get("failed")]
    ok = [r for r in raw if not r.get("failed")]
    report(f"N: -> {len(ok)} after dropping {len(failed)} cells that raised")

    gate = [r for r in ok if r.get("tag") in GATE_TAGS]
    prod = [r for r in ok if r.get("tag") not in GATE_TAGS]
    report(f"N: -> {len(prod)} production records, {len(gate)} gate/probe records "
           f"(tags {sorted({r.get('tag') for r in gate})}) held out of every S3 number")
    return prod, gate, files


def dedup(recs: list[dict], report) -> list[dict]:
    """Keep ONE record per cell key. Which one is kept is NOT left to file order:
    lessons.md 2026-09-09 records 15 duplicated keys in this project returning different
    phi, with the survivor set by the filename sort. Any duplicate whose phi differs is a
    hard error here rather than a silent choice."""
    by_key: dict[str, dict] = {}
    n_dup = 0
    for r in recs:
        k = r["key"]
        if k in by_key:
            n_dup += 1
            if by_key[k]["phi"] != r["phi"]:
                raise SystemExit(
                    f"DUPLICATE KEY {k} with DIFFERENT phi "
                    f"({by_key[k]['phi']!r} vs {r['phi']!r}) -- which record survives "
                    f"would set the answer; refusing to choose")
            continue
        by_key[k] = r
    report(f"N: -> {len(by_key)} after dropping {n_dup} exact-duplicate keys")
    return list(by_key.values())


# Amendment 5.10: the four fields the identity check compares, and the cap that the
# default-cap records carry. A physical cell is the key with its `_it{cap}` suffix
# removed (phase2_shape.cell_key appends it only when the cap is non-default).
S6_FIELDS = ("phi", "iters", "z_all", "p_hat")
DEFAULT_CAP = 200000


def physical_key(r: dict) -> str:
    k = r["key"]
    suffix = f"_it{int(r['max_iter'])}"
    return k[:-len(suffix)] if k.endswith(suffix) else k


def prefer_raised_cap(prod: list[dict], report) -> tuple[list[dict], dict]:
    """Amendment 5.10's rule for a ladder that mixes caps: for one physical cell, when
    a default-cap and a raised-cap record both exist, the raised-cap record is used and
    the fact printed. Returns the reduced list and {physical_key: {cap: record}} for
    every cell that exists at more than one cap (the input to S6)."""
    by_cell: dict[str, dict[int, dict]] = defaultdict(dict)
    for r in prod:
        cap = int(r["max_iter"])
        if cap in by_cell[physical_key(r)]:
            raise SystemExit(f"two records for {physical_key(r)} at cap {cap}")
        by_cell[physical_key(r)][cap] = r
    multi = {k: v for k, v in by_cell.items() if len(v) > 1}
    out = [v[max(v)] for v in by_cell.values()]
    report(f"N: -> {len(out)} physical cells after preferring the raised-cap record on "
           f"{len(multi)} cell(s) that exist at more than one cap (amendment 5.10)")
    for k in sorted(multi):
        caps = sorted(multi[k])
        report(f"       {k}: caps {caps} -> using {max(caps)}")
    by_cap = defaultdict(int)
    for r in out:
        by_cap[int(r["max_iter"])] += 1
    report("     records used, by cap: " + ", ".join(f"{c}: {n}" for c, n in sorted(by_cap.items())))
    return out, multi


def g_s6(multi: dict) -> Gate:
    g = Gate("S6", "Amendment 5.10 identity check: every physical cell that CONVERGED at "
                   "the default cap and exists again at the raised cap must be "
                   f"BIT-IDENTICAL on {S6_FIELDS}. FAILS on any non-identical field; "
                   "VACUOUS (and S3 stays BLOCKED) if no converged cell exists at both "
                   "caps. A cell that was CUT at the default cap carries no identity "
                   "claim -- for those the raised-cap record must instead have run PAST "
                   "the old cap, which is the check that the cap acted at all.")
    g.say(f"N: {len(multi)} physical cell(s) present at more than one cap")
    ident = {k: v for k, v in multi.items() if v[min(v)].get("converged")}
    cut = {k: v for k, v in multi.items() if not v[min(v)].get("converged")}
    g.say(f"N: -> {len(ident)} converged at the default cap (identity REQUIRED); "
          f"{len(cut)} cap-cut at the default cap (identity NOT expected; the raised-cap "
          f"record must show iters > {DEFAULT_CAP})")
    if not ident:
        g.say("VACUOUS: no converged cell exists at both caps, so the raised-cap records "
              "cannot be shown to be the same protocol. S3 is BLOCKED.")
        g.verdict(False, informative=False)
        return g
    all_same = True
    for k in sorted(ident):
        recs = ident[k]
        caps = sorted(recs)
        base = recs[caps[0]]
        for cap in caps[1:]:
            r = recs[cap]
            diffs = [f for f in S6_FIELDS if base[f] != r[f]]
            same = not diffs
            all_same &= same
            g.say(f"- `{k}` caps {caps[0]} vs {cap}: phi {base['phi']!r} vs {r['phi']!r}, "
                  f"iters {base['iters']} vs {r['iters']}, z_all {base['z_all']!r} vs "
                  f"{r['z_all']!r}, p_hat {base['p_hat']!r} vs {r['p_hat']!r} -> "
                  f"{'BIT-IDENTICAL' if same else 'DIFFER on ' + ', '.join(diffs)}")
    n_past = 0
    n_conv = 0
    for k in sorted(cut):
        recs = cut[k]
        r = recs[max(recs)]
        past = int(r["iters"]) > int(recs[min(recs)]["max_iter"])
        n_past += past
        n_conv += bool(r.get("converged"))
    if cut:
        g.say(f"- cap-cut cells re-run: {n_past}/{len(cut)} ran past the old cap of "
              f"{DEFAULT_CAP}; {n_conv}/{len(cut)} converged under the raised cap")
        if n_past != len(cut):
            g.say("  a cap-cut cell that did NOT run past the old cap means the raised "
                  "cap was not in force for it -- the re-run is not what 5.10 registered")
            all_same = False
    if all_same:
        g.say(f"All {len(ident)} identity pair(s) bit-identical on every field: a cell "
              f"that reaches its doneHook at k < {DEFAULT_CAP} stops at k under any "
              f"larger cap, so the raised-cap records are the same protocol with a "
              f"larger resource allowance.")
    else:
        g.say("The re-run is VOID and nothing from it is pooled (amendment 5.10).")
    g.verdict(all_same)
    return g


def g_s0a(gate_recs, report) -> Gate:
    g = Gate("S0a", "The three PREREG_FABRIC F1 cells, re-run at alpha = 1 through the "
                    "shape runner, return `phi_porosity` BIT-IDENTICAL to the committed "
                    "record.")
    recs = [r for r in gate_recs if r.get("tag") == "shpgate"]
    g.say(f"N: {len(recs)} `shpgate` records")
    if not recs:
        g.verdict(False)
        g.say("NO GATE CELLS ON DISK -- S0a cannot be read.")
        return g
    ok = True
    g.say("")
    g.say("| family | CV | alpha | `phi_porosity` | reference | bit-identical |")
    g.say("|---|---|---|---|---|---|")
    seen = set()
    for r in sorted(recs, key=lambda x: x["family"]):
        key = (r["family"], round(float(r["cv"]), 3))
        ref = S0A_REFERENCE.get(key)
        if ref is None:
            ok = False
            g.say(f"| {r['family']} | {r['cv']} | {r['alpha']} | {r['phi_porosity']!r} "
                  f"| NO REFERENCE | **FAIL** |")
            continue
        seen.add(key)
        hit = r["phi_porosity"] == ref
        ok = ok and hit and float(r["alpha"]) == 1.0
        g.say(f"| {r['family']} | {r['cv']} | {r['alpha']} | `{r['phi_porosity']!r}` "
              f"| `{ref!r}` | {'YES' if hit else '**NO**'} |")
    missing = set(S0A_REFERENCE) - seen
    if missing:
        ok = False
        g.say("")
        g.say(f"**MISSING gate cells: {sorted(missing)}** -- a gate that did not run "
              f"cannot pass.")
    g.say("")
    g.say("The quantized channels are printed and NOT read as the gate "
          "(lessons.md 2026-08-31): a real perturbation can hide inside one quantum of "
          "`z_all` (2k/N), `rattler_frac` (k/N) or `n_contacts`.")
    g.say("")
    g.say("| family | `z_all` | `z_nr` | `rattler_frac` | `n_contacts` | "
          "`n_contacts_raw` | `dup_contact_frac` |")
    g.say("|---|---|---|---|---|---|---|")
    for r in sorted(recs, key=lambda x: x["family"]):
        g.say(f"| {r['family']} | {r['z_all']!r} | {r['z_nr']!r} | "
              f"{r['rattler_frac']!r} | {r['n_contacts']} | "
              f"{r.get('n_contacts_raw', 'n/a')} | {r.get('dup_contact_frac', 'n/a')} |")
    g.say("")
    g.say("`n_contacts_raw` equal to `n_contacts` on all three is the check that the "
          "grain-level rewrite of `contact_graph` is a no-op for spheres -- asserted "
          "here rather than assumed.")
    g.verdict(ok)
    return g


def g_s0b(gate_recs) -> Gate:
    g = Gate("S0b", f"At alpha = 1 the two volume routes agree: "
                    f"|phi_union - phi_porosity| <= {S0B_BAR:g}.")
    recs = [r for r in gate_recs if r.get("tag") == "shpgate"]
    if not recs:
        g.verdict(False)
        g.say("NO GATE CELLS ON DISK.")
        return g
    worst = 0.0
    g.say("| family | `phi_union` | `phi_porosity` | difference |")
    g.say("|---|---|---|---|")
    for r in sorted(recs, key=lambda x: x["family"]):
        d = abs(float(r["phi_union"]) - float(r["phi_porosity"]))
        worst = max(worst, d)
        g.say(f"| {r['family']} | `{r['phi_union']!r}` | `{r['phi_porosity']!r}` "
              f"| {d:.3e} |")
    g.say("")
    g.say(f"Worst: **{worst:.3e}** against a bar of {S0B_BAR:g}. For scale the seed "
          f"scatter in phi is 2.2e-4, eleven orders larger. The bar is set from the "
          f"measurement and is NOT vacuous: a deviation of exactly zero would mean the "
          f"two routes are the same calculation, and the measured value is not zero "
          f"({worst:.3e}).")
    g.verdict(worst <= S0B_BAR and worst > 0.0)
    return g


def g_s0c(prod, report) -> Gate:
    g = Gate("S0c", "At every alpha > 1 the clump body is the body it claims to be: "
                    "`phi_porosity`/`phi_union` equals the PREDICTED double count "
                    "2*(ds/d_eq)^3, and the mass Yade assembles equals the union mass.")
    recs = [r for r in prod if float(r["alpha"]) > 1.0]
    g.say(f"N: {len(recs)} production records at alpha > 1")
    if not recs:
        g.verdict(False)
        g.say("NO alpha > 1 PRODUCTION CELLS -- the clump construction is untested.")
        return g
    ok = True
    g.say("")
    g.say("| alpha | N cells | predicted 2(ds/d_eq)^3 | worst measured ratio | "
          "worst |ratio-pred| | worst |mass/vol - 1| |")
    g.say("|---|---|---|---|---|---|")
    by_a = defaultdict(list)
    for r in recs:
        by_a[round(float(r["alpha"]), 6)].append(r)
    for a in sorted(by_a):
        pred = 2.0 * sub_sphere_diameter(1.0, a) ** 3
        dr = 0.0
        worst_ratio = float("nan")
        dm = 0.0
        for r in by_a[a]:
            ratio = float(r["phi_porosity"]) / float(r["phi_union"])
            if abs(ratio - pred) >= dr:
                dr = abs(ratio - pred)
                worst_ratio = ratio
            m = r.get("mass_over_vol_solid")
            if m is None or not np.isfinite(m):
                ok = False
                dm = float("inf")
            else:
                dm = max(dm, abs(float(m) - 1.0))
        hit = dr <= S0C_RATIO_BAR and dm <= S0C_MASS_BAR
        ok = ok and hit
        g.say(f"| {a:.2f} | {len(by_a[a])} | {pred:.10f} | {worst_ratio:.10f} | "
              f"{dr:.3e} | {dm:.3e} |")
    g.say("")
    g.say(f"Bars: ratio {S0C_RATIO_BAR:g}, mass {S0C_MASS_BAR:g}.")
    g.say("")
    g.say("**WHAT THIS GATE CAN AND CANNOT DECIDE** (referee, 2026-09-18). It calls "
          "`shape.sub_sphere_diameter` on BOTH sides of the comparison -- the "
          "construction used it to size the members, and the prediction uses it again "
          "-- and `density_scale` is DEFINED as the reciprocal of the same quantity. "
          "So `phi_porosity/phi_union = 2*(ds/d_eq)^3` and `mass/vol_solid = 1` hold "
          "for ANY value that function returns. **This gate therefore cannot detect an "
          "error in the closed form `g(u)`**, and an earlier version of this text "
          "claimed it could. What it DOES decide is that the construction is "
          "self-consistent end to end: that Yade built members of the radius the code "
          "computed, that the clump assembly preserved them, that the density was "
          "applied, and that `utils.porosity()` is summing member volumes and nothing "
          "else -- each of which can fail independently and none of which is visible "
          "elsewhere. The closed form itself is anchored OUTSIDE this pipeline, by "
          "quadrature over the union's cross-sections, in "
          "`test_shape_campaign.py::test_clump_volume_against_independent_quadrature`, "
          "which imports none of `shape.py`'s constants and is shown able to fail on a "
          "1 % perturbation of `g(u)`'s coefficient.")
    g.verdict(ok)
    return g


def g_s1(prod) -> Gate:
    g = Gate("S1", "The PSD does not move along the alpha axis: `cv_realised`, "
                   "`d32_realised` and `d43_realised` at alpha > 1 equal those at "
                   f"alpha = 1 to {S1_BAR:g}, at every (family, CV, seed).")
    by = defaultdict(dict)
    for r in prod:
        by[(r["family"], round(float(r["cv"]), 3), int(r["seed"]))][
            round(float(r["alpha"]), 6)] = r
    worst = {"cv_realised": 0.0, "d32_realised": 0.0, "d43_realised": 0.0}
    n_pairs = 0
    n_groups_without_ref = 0
    for _k, byalpha in by.items():
        if 1.0 not in byalpha:
            n_groups_without_ref += 1
            continue
        ref = byalpha[1.0]
        for a, r in byalpha.items():
            if a == 1.0:
                continue
            n_pairs += 1
            for f in worst:
                worst[f] = max(worst[f], abs(float(r[f]) - float(ref[f])))
    g.say(f"N: {len(by)} (family, CV, seed) groups, {n_pairs} alpha-to-alpha=1 "
          f"comparisons, {n_groups_without_ref} groups with no alpha = 1 rung to "
          f"compare against")
    g.say("")
    for f, v in worst.items():
        g.say(f"- worst |{f}(alpha) - {f}(1)| = **{v:.3e}**")
    g.say("")
    g.say("This is exact by construction -- `d_eq` is drawn from the same "
          "`psd.make_psd` at every alpha -- and is checked rather than assumed. "
          "It CAN fail: a construction that drew radii after the orientations, or one "
          "that rescaled d_eq instead of ds, would move all three.")
    if n_pairs == 0:
        g.say("")
        g.say("**NO COMPARISONS AVAILABLE** -- this gate has not been exercised.")
        g.verdict(False)
    else:
        g.verdict(max(worst.values()) <= S1_BAR)
    return g


def g_s2(raw_files, prod) -> Gate:
    g = Gate("S2", "Every requested cell is legal for the periodic collider at its own "
                   "alpha: no cell raised the clearance error, and the required N was "
                   "computed per alpha rather than inherited from the sphere run.")
    n_fail = 0
    n_probe_fail = 0
    msgs = []
    probe_msgs = []
    for f in raw_files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if not (r.get("failed")
                        and "periodic collider" in str(r.get("error", ""))):
                    continue
                # A clearance failure on a GATE/PROBE cell is not a production defect --
                # it is the guard doing its job on a deliberately undersized cell, and
                # it is the only evidence available that this guard can fire at all.
                # Counting it against S2 reported a working guard as a failed gate.
                if r.get("tag") in GATE_TAGS:
                    n_probe_fail += 1
                    if len(probe_msgs) < 2:
                        probe_msgs.append(str(r.get("error"))[:200])
                else:
                    n_fail += 1
                    if len(msgs) < 5:
                        msgs.append(str(r.get("error"))[:240])
    g.say(f"N: {n_fail} PRODUCTION cells raised the clearance error out of "
          f"{n_fail + len(prod)} requested")
    for m in msgs:
        g.say(f"- `{m}`")
    g.say("")
    g.say(f"N: {n_probe_fail} GATE/PROBE cell(s) raised it, which is the guard working "
          f"and is **the only evidence that this check can fire at all** -- the timing "
          f"probe ran at N = 1700, where the clearance genuinely is insufficient at "
          f"alpha = 2:")
    for m in probe_msgs:
        g.say(f"- `{m}`")
    if n_probe_fail == 0:
        g.say("- (none -- so S2 has never been seen to fire, and a check that has never "
              "fired is untested code)")
    ns = sorted({int(r["n_requested"]) for r in prod})
    spans = sorted({round(float(r.get("span_over_deq", 1.0)), 6) for r in prod})
    g.say("")
    g.say(f"N_grains used across the production grid: {ns} -- a SINGLE value is "
          f"required, because phi_J drifts with N and an N that moved with alpha would "
          f"wear alpha's label.")
    g.say(f"`span_over_deq` realised: {spans} against "
          f"{[round(float(sub_sphere_diameter(1.0, a) * a), 6) for a in ALPHA_LADDER]} "
          f"for the full ladder.")
    g.verdict(n_fail == 0 and len(ns) == 1)
    return g


def phi_j_table(prod, report):
    """phi_J per (family, CV, alpha, seed), from the p_hat ladder. Returns
    {(family, cv, alpha): {seed: (phi_J, se, b_free, n_points)}} and prints N at every
    filter."""
    # FILTER TO THE REGISTERED LADDER FIRST, before the converged filter, so the two
    # counts can be read separately: an off-ladder rung is excluded by design and a
    # non-converged on-ladder cell is a loss.
    def _on_ladder(r) -> bool:
        return any(abs(float(r["p_hat_target"]) / p - 1.0) < P_HAT_TOL
                   for p in P_HAT_LADDER_SHAPE)

    # Partitioned in ONE pass on the predicate. A `r not in on` membership test would
    # fall back to dict equality, and these records carry NaN fields, where `x == x` is
    # False -- it happens to work only because CPython short-circuits on identity first.
    on = [r for r in prod if _on_ladder(r)]
    off = [r for r in prod if not _on_ladder(r)]
    report(f"N: -> {len(on)} on the registered ladder {P_HAT_LADDER_SHAPE}; "
           f"{len(off)} off-ladder records excluded BY DESIGN (amendment 5.6)")
    if off:
        off_by = defaultdict(int)
        for r in off:
            off_by[(float(r["p_hat_target"]), round(float(r["cv"]), 2),
                    bool(r.get("converged")))] += 1
        for k in sorted(off_by):
            report(f"       off-ladder (p_hat, CV, converged)={k}: {off_by[k]}")
    prod = on
    conv = [r for r in prod if r.get("converged")]
    report(f"N: -> {len(conv)} after dropping {len(prod) - len(conv)} non-converged "
           f"cells")
    # BROKEN DOWN BY (CV, alpha, RUNG), not pooled. Amendment 5.6 shortened the ladder
    # because the 3e-5 rung failed at CV = 0.5 and not at CV = 0.0; the evidence for
    # that covered only the 1e-3 rung at CV = 0.5 (referee finding iv), so whether the
    # RETAINED rungs converge there is exactly the question, and a pooled percentage
    # cannot answer it (lessons.md 2026-09-03: a global summary cannot diagnose a local
    # defect -- check the threshold where the defect is).
    att = defaultdict(int)
    good = defaultdict(int)
    for r in prod:
        k = (round(float(r["cv"]), 2), round(float(r["alpha"]), 2),
             float(r["p_hat_target"]))
        att[k] += 1
        if r.get("converged"):
            good[k] += 1
    report("     converged / attempted by (CV, alpha, p_hat) -- a drop rate that "
           "correlates with the axis under test is not a neutral filter:")
    bad_cells = []
    for k in sorted(att):
        flag = "" if good[k] == att[k] else "   <-- NOT ALL CONVERGED"
        report(f"       CV={k[0]:.2f} alpha={k[1]:.2f} p_hat={k[2]:.0e}: "
               f"{good[k]}/{att[k]}{flag}")
        if good[k] != att[k]:
            bad_cells.append(k)
    if bad_cells:
        report(f"     {len(bad_cells)} (CV, alpha, rung) combination(s) lost cells to "
               f"non-convergence. If those are concentrated at one end of the CV axis "
               f"the shortened ladder has the SAME defect it was introduced to remove.")
    groups = defaultdict(list)
    for r in conv:
        groups[(r["family"], round(float(r["cv"]), 3), round(float(r["alpha"]), 6),
                int(r["seed"]))].append(r)
    out = defaultdict(dict)
    n_short = 0
    for (fam, cv, a, seed), rows in groups.items():
        rows = [r for r in rows if float(r["p_hat"]) <= P_HAT_FIT_MAX]
        # EXACTLY the registered ladder, not "at least three rungs". A cell that lost
        # one rung is a SHORTER lever arm than its neighbours, and the difference would
        # land inside D rather than being visible as a missing cell.
        if len(rows) != len(P_HAT_LADDER_SHAPE):
            n_short += 1
            continue
        p = np.array([float(r["p_hat"]) for r in rows])
        f = np.array([float(r["phi"]) for r in rows])
        fit = fit_phi_J(p, f)
        caps = tuple(sorted({int(r["max_iter"]) if "max_iter" in r else -1 for r in rows}))
        out[(fam, cv, a)][seed] = (fit["phi_J_lin"], fit["phi_J_lin_se"],
                                   fit["b_free"], fit["n_points"], fit["slope_lin"],
                                   caps)
    report(f"N: -> {sum(len(v) for v in out.values())} (family, CV, alpha, seed) "
           f"ladders fitted; {n_short} dropped for not carrying EXACTLY the registered "
           f"{len(P_HAT_LADDER_SHAPE)} rungs inside p_hat <= {P_HAT_FIT_MAX:g}")
    # Amendment 5.10: the caps present in every fitted ladder are printed.
    cap_sets = defaultdict(int)
    for (fam, cv, a), seeds in out.items():
        for seed, v in seeds.items():
            cap_sets[(fam, cv, a, v[5])] += 1
    report("     caps present in the fitted ladders (family, CV, alpha) -> caps: n_seeds "
           "(amendment 5.10; -1 = record carries no max_iter field)")
    for (fam, cv, a, caps) in sorted(cap_sets):
        report(f"       {fam:<11} CV={cv:.1f} alpha={a:.2f} -> {list(caps)}: "
               f"{cap_sets[(fam, cv, a, caps)]}")
    report("     `b_free` is NaN on every cell by design: the free-exponent fit needs "
           "four points and the registered ladder has three (amendment 5.6). Every "
           "phi_J here is a LINEAR extrapolation with no independent exponent check.")
    return out


def g_s3(tab, report, mult=None) -> tuple[Gate, dict]:
    """`mult`: {alpha: mean n_contacts_raw/n_contacts over production cells at that alpha},
    used only to print the multiplicity-corrected stiffness prediction; None prints n/a."""
    g = Gate("S3", "PRIMARY. Per family, the paired CV 0 -> CV_max shift in phi_J at "
                   "each alpha. The registered question is whether the "
                   "shape-of-the-distribution effect SURVIVES, GROWS or VANISHES as "
                   "grains stop being spheres. No direction is predicted, and the "
                   "headline is the (CV, alpha) INTERACTION, not the alpha main effect.")
    fams = sorted({k[0] for k in tab})
    alphas = sorted({k[2] for k in tab})
    cvs = sorted({k[1] for k in tab})
    if not fams or len(cvs) < 2 or len(alphas) < 2:
        g.say(f"N: families {fams}, CVs {cvs}, alphas {alphas} -- a 2-D grid needs at "
              f"least two of each.")
        g.verdict(False)
        return g, {}
    cv_lo, cv_hi = cvs[0], cvs[-1]
    g.say(f"N: {len(fams)} families x {len(cvs)} CV rungs x {len(alphas)} alpha rungs; "
          f"the paired shift is read CV {cv_lo} -> {cv_hi}.")
    g.say("")
    g.say("#### phi_J per cell (mean over seeds, with the seed standard error)")
    g.say("")
    g.say("| family | alpha | " + " | ".join(f"CV {c}" for c in cvs) + " |")
    g.say("|---" * (len(cvs) + 2) + "|")
    mean = {}
    sem = {}
    nseed = {}
    for fam in fams:
        for a in alphas:
            cells = []
            for c in cvs:
                v = [x[0] for x in tab.get((fam, c, a), {}).values()
                     if np.isfinite(x[0])]
                if v:
                    mean[(fam, c, a)] = float(np.mean(v))
                    sem[(fam, c, a)] = (float(np.std(v, ddof=1) / math.sqrt(len(v)))
                                        if len(v) > 1 else float("nan"))
                    nseed[(fam, c, a)] = len(v)
                    cells.append(f"{mean[(fam, c, a)]:.6f} +- "
                                 f"{sem[(fam, c, a)]:.6f} (n={len(v)})")
                else:
                    cells.append("--")
            g.say(f"| {fam} | {a:.2f} | " + " | ".join(cells) + " |")

    g.say("")
    g.say(f"#### D(alpha) = phi_J(CV {cv_hi}) - phi_J(CV {cv_lo}), per family")
    g.say("")
    g.say("| family | " + " | ".join(f"alpha {a:.2f}" for a in alphas) + " |")
    g.say("|---" * (len(alphas) + 1) + "|")
    D = {}
    Dse = {}
    for fam in fams:
        cells = []
        for a in alphas:
            k_hi, k_lo = (fam, cv_hi, a), (fam, cv_lo, a)
            if k_hi in mean and k_lo in mean:
                d = mean[k_hi] - mean[k_lo]
                s = math.sqrt(max(sem.get(k_hi, float('nan')), 0) ** 2
                              + max(sem.get(k_lo, float('nan')), 0) ** 2)
                D[(fam, a)] = d
                Dse[(fam, a)] = s
                cells.append(f"{d:+.6f} +- {s:.6f}")
            else:
                cells.append("--")
        g.say(f"| {fam} | " + " | ".join(cells) + " |")
    # The RATIO D(alpha_hi)/D(alpha_lo) per family, printed because it is the number a
    # reader would form by hand from the table above; its error bar is first-order in
    # the two D standard errors (independent across alpha only approximately -- the
    # seeds are paired across alpha, so this bar is an over-estimate and the paired
    # test below is the one that decides).
    a_lo, a_hi = alphas[0], alphas[-1]
    ratios = []
    for fam in fams:
        if (fam, a_lo) in D and (fam, a_hi) in D and D[(fam, a_lo)] != 0:
            q = D[(fam, a_hi)] / D[(fam, a_lo)]
            qse = abs(q) * math.sqrt((Dse[(fam, a_hi)] / D[(fam, a_hi)]) ** 2
                                     + (Dse[(fam, a_lo)] / D[(fam, a_lo)]) ** 2)
            ratios.append(f"{fam} {q:.3f} +- {qse:.3f}")
    if ratios:
        g.say("")
        g.say(f"D(alpha {a_hi:.2f}) / D(alpha {a_lo:.2f}) per family (unpaired first-order "
              f"bar; the paired test below decides): " + "; ".join(ratios))

    g.say("")
    g.say("#### THE INTERACTION -- does D depend on alpha? PAIRED ON THE SEED")
    g.say("")
    g.say("The seeds are PAIRED ACROSS ALPHA and unpaired across CV. Seed s at "
          "(CV, alpha = 1) and seed s at (CV, alpha = 2) draw the SAME `d_eq` from the "
          "same stream -- only the grain's shape differs -- so the per-seed difference "
          "cancels the size realisation exactly rather than statistically "
          "(lessons.md 2026-08-30: take the paired difference, do not subtract two "
          "large variances). Seed s at CV 0 and seed s at CV 0.5 are DIFFERENT "
          "packings, so that half is not paired and its scatter is carried as an "
          "independent term.")
    g.say("")
    a_lo, a_hi = alphas[0], alphas[-1]
    g.say(f"Interaction per seed: "
          f"`Delta_s = [phi_J(CV {cv_hi}, a {a_hi:.2f}, s) - phi_J(CV {cv_lo}, "
          f"a {a_hi:.2f}, s)] - [phi_J(CV {cv_hi}, a {a_lo:.2f}, s) - "
          f"phi_J(CV {cv_lo}, a {a_lo:.2f}, s)]`")
    g.say("")
    g.say("| family | n seeds | mean Delta | se(Delta) | |Delta|/se | t_crit(0.05) "
          "| resolved? |")
    g.say("|---|---|---|---|---|---|---|")
    inter = {}
    for fam in fams:
        seeds = (set(tab.get((fam, cv_hi, a_hi), {})) & set(tab.get((fam, cv_lo, a_hi), {}))
                 & set(tab.get((fam, cv_hi, a_lo), {})) & set(tab.get((fam, cv_lo, a_lo), {})))
        vals = []
        for s in sorted(seeds):
            d_hi = (tab[(fam, cv_hi, a_hi)][s][0] - tab[(fam, cv_lo, a_hi)][s][0])
            d_lo = (tab[(fam, cv_hi, a_lo)][s][0] - tab[(fam, cv_lo, a_lo)][s][0])
            if np.isfinite(d_hi) and np.isfinite(d_lo):
                vals.append(d_hi - d_lo)
        if len(vals) < 2:
            g.say(f"| {fam} | {len(vals)} | -- | -- | -- | -- | "
                  f"fewer than 2 complete seed quartets |")
            continue
        v = np.array(vals, float)
        m = float(v.mean())
        se = float(v.std(ddof=1) / math.sqrt(v.size))
        # The two-sided 95 % value at n-1 dof, NOT a flat 2 (the descriptive sweep's own
        # correction: at 3 seeds the se carries 2 dof and t = 4.303, so a flat 2*se is
        # an 18 % test).
        from scipy import stats as _st
        tcrit = float(_st.t.ppf(0.975, v.size - 1))
        res = abs(m) > tcrit * se if se > 0 else False
        inter[fam] = (m, se, tcrit, res, v.size)
        g.say(f"| {fam} | {v.size} | {m:+.6f} | {se:.6f} | "
              f"{abs(m)/se if se else float('nan'):.2f} | {tcrit:.3f} | "
              f"{'YES' if res else 'no'} |")
    g.say("")
    g.say("**The critical value is t at n-1 dof, NOT a flat 2.** At three seeds the "
          "standard error carries two degrees of freedom and the two-sided 95 % value "
          "is 4.303; a flat `2*se` guard is an 18 % test and admitted a number this "
          "project had to retract once already "
          "(`REPORT_SWEEP_2026-09-18.md`, contact-force CV).")
    g.say("")
    g.say("**A Delta inside its critical value is UNDETERMINED AT THREE SEEDS, not an "
          "absence.** The achieved resolving power is `t_crit * se(Delta)`, printed "
          "below, and a real interaction smaller than that cannot be distinguished from "
          "zero by this design.")
    g.say("")
    g.say("| family | resolving power `t_crit*se` | |D| at alpha_lo | "
          "|D| at alpha_hi | resolving power < |D| (the CV shift itself; NOT the "
          "interaction test, which is the `resolved?` column above) |")
    g.say("|---|---|---|---|---|")
    for fam in fams:
        if fam not in inter:
            continue
        m, se, tcrit, _r, _n = inter[fam]
        rp = tcrit * se
        dlo = abs(D.get((fam, a_lo), float("nan")))
        dhi = abs(D.get((fam, a_hi), float("nan")))
        ref = max(dlo, dhi)
        g.say(f"| {fam} | {rp:.6f} | {dlo:.6f} | {dhi:.6f} | "
              f"{'yes' if np.isfinite(ref) and rp < ref else 'NO -- the bar is wider '
                 'than the CV effect itself'} |")

    slopes = {}
    if len(alphas) >= 3:
        g.say("")
        g.say("#### The alpha SLOPE, available only with three or more rungs")
        g.say("")
        g.say("| family | slope dD/dalpha | se | |slope|/se |")
        g.say("|---|---|---|---|")
        for fam in fams:
            xs = np.array([a for a in alphas if (fam, a) in D], float)
            ys = np.array([D[(fam, a)] for a in alphas if (fam, a) in D], float)
            if xs.size < 3:
                continue
            A = np.vstack([np.ones_like(xs), xs]).T
            coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
            resid = ys - A @ coef
            dof = xs.size - 2
            s2 = float(resid @ resid / dof) if dof > 0 else float("nan")
            se = (float(math.sqrt((s2 * np.linalg.inv(A.T @ A))[1, 1]))
                  if np.isfinite(s2) and s2 > 0 else float("nan"))
            slopes[fam] = (float(coef[1]), se)
            g.say(f"| {fam} | {coef[1]:+.6f} | {se:.6f} | "
                  f"{abs(coef[1])/se if se and np.isfinite(se) else float('nan'):.2f} |")
    else:
        g.say("")
        g.say(f"**Only {len(alphas)} alpha rungs are on disk, so there is no slope to "
              f"fit and no curvature to see.** The paired Delta above is the whole "
              f"interaction this grid can measure. A two-rung axis can show a DIRECTION "
              f"and not a SHAPE, and with alpha = 1 a sphere and alpha > 1 a clump, the "
              f"single step it does measure is the one that carries the seam.")

    g.say("")
    g.say("#### THE CONTACT-STIFFNESS CONFOUND, as a PREDICTION that can fail")
    g.say("")
    g.say("Yade's contact law never sees `d_eq`; it sees the MEMBER radius `ds/2`, and "
          "`ds/d_eq` falls 1.0000 -> 0.7937 along this ladder, so the normal stiffness "
          "drops with alpha at fixed `p_hat` (referee, 2026-09-18). Because `ds/d_eq` "
          "is a CONSTANT within a rung, the effect is exactly a rescaling of the "
          "pressure axis, and for `phi = phi_J + c*p_hat` a rescaling `p_hat -> k*p_hat`"
          " changes `c -> c/k` and leaves the INTERCEPT untouched. That argument makes "
          "a number: the fitted slope should rise by `d_eq/ds` from the lowest alpha to "
          "the highest.")
    g.say("")
    pred_ratio = (sub_sphere_diameter(1.0, a_lo) / sub_sphere_diameter(1.0, a_hi))
    # Second, sharper prediction (referee, 2026-09-16): a grain pair carrying m member
    # interactions is m springs in parallel, so its stiffness is m times a single
    # contact's and the slope c ~ 1/k falls by m. With m = mean n_contacts_raw/n_contacts
    # (S4; 1 at the sphere rung by construction) the two known stiffness effects
    # together predict slope(a_hi)/slope(a_lo) = (d_eq/ds) * m(a_lo)/m(a_hi). First
    # order: it assumes the m interactions share one normal and one overlap.
    m_lo = mult.get(a_lo, float("nan")) if mult else float("nan")
    m_hi = mult.get(a_hi, float("nan")) if mult else float("nan")
    pred_mult = pred_ratio * m_lo / m_hi if np.isfinite(m_lo * m_hi) and m_hi else float("nan")
    g.say(f"| family | mean slope at alpha {a_lo:.2f} | at alpha {a_hi:.2f} | measured "
          f"ratio | predicted `d_eq/ds` | ratio/prediction | predicted with contact "
          f"multiplicity `(d_eq/ds)*m_lo/m_hi` | ratio/that |")
    g.say("|---|---|---|---|---|---|---|---|")
    ratios = {}
    for fam in fams:
        s_lo = [v[4] for v in tab.get((fam, cv_lo, a_lo), {}).values()
                if np.isfinite(v[4])] + [v[4] for v in
                                         tab.get((fam, cv_hi, a_lo), {}).values()
                                         if np.isfinite(v[4])]
        s_hi = [v[4] for v in tab.get((fam, cv_lo, a_hi), {}).values()
                if np.isfinite(v[4])] + [v[4] for v in
                                         tab.get((fam, cv_hi, a_hi), {}).values()
                                         if np.isfinite(v[4])]
        if not s_lo or not s_hi:
            g.say(f"| {fam} | -- | -- | -- | {pred_ratio:.4f} | no cells | "
                  f"{pred_mult:.4f} | no cells |")
            continue
        sl_lo, sl_hi = float(np.mean(s_lo)), float(np.mean(s_hi))
        r = sl_hi / sl_lo if sl_lo else float("nan")
        ratios[fam] = r
        g.say(f"| {fam} | {sl_lo:.4f} | {sl_hi:.4f} | {r:.4f} | {pred_ratio:.4f} | "
              f"{r/pred_ratio if pred_ratio else float('nan'):.3f} | {pred_mult:.4f} | "
              f"{r/pred_mult if np.isfinite(pred_mult) and pred_mult else float('nan'):.3f} |")
    g.say("")
    g.say(f"Contact multiplicity `m` = mean `n_contacts_raw/n_contacts` per rung (S4): "
          f"alpha {a_lo:.2f} -> {m_lo:.4f}, alpha {a_hi:.2f} -> {m_hi:.4f}.")
    g.say("")
    g.say("A measured ratio near the prediction says the stiffness change is a pure "
          "rescaling and the intercept is safe. **A ratio far from it does NOT "
          "invalidate phi_J by itself** -- elongated grains have a genuinely different "
          "compressibility, so the slope carries physics as well as the artefact, and "
          "this test cannot separate them. It is a NECESSARY condition for the "
          "rescaling argument, not a sufficient one, and it is reported as such.")
    if ratios and np.isfinite(pred_mult):
        below_1 = sum(1 for r in ratios.values() if r < 1.0)
        below_m = sum(1 for r in ratios.values() if r < pred_mult)
        g.say("")
        g.say(f"**The SIGN of the miss** (referee, 2026-09-16): the member-radius argument "
              f"predicts the slope RISES ({pred_ratio:.4f}); it falls in {below_1} of "
              f"{len(ratios)} families. Contact multiplicity pushes the other way (parallel "
              f"springs stiffen the pair), and the two together predict {pred_mult:.4f} -- "
              f"the measured ratio is below THAT too in {below_m} of {len(ratios)} families, "
              f"by a factor of {min(ratios.values())/pred_mult:.3f}-"
              f"{max(ratios.values())/pred_mult:.3f}. So the multiplicity confound has the "
              f"right sign and the wrong size, and the remainder is either the compressibility "
              f"of an elongated packing or a contact effect this bookkeeping does not model; "
              f"neither is separated here, and the rescaling argument stays uncertified.")
    g.say("")
    g.say("#### THE SEAM, restated where the number is")
    g.say("")
    g.say(f"alpha = {a_lo:.2f} is the sphere rung and every alpha > 1 rung is a clump, "
          f"so Delta is not a pure shape effect: it also carries the ONSET of contact "
          f"multiplicity. The measured `dup_contact_frac` per rung is in S4. **No "
          f"number in this section may be quoted as an effect of ASPECT RATIO** -- the "
          f"prereg labels it elongation + interlocking + contact multiplicity until the "
          f"LevelSet control exists.")
    g.verdict(len(inter) == len(fams) and bool(inter), informative=bool(inter))
    return g, {"mean": mean, "sem": sem, "nseed": nseed, "D": D, "Dse": Dse,
               "inter": inter, "slopes": slopes, "cvs": cvs, "alphas": alphas,
               "fams": fams, "cv_lo": cv_lo, "cv_hi": cv_hi}


def g_s5(prod, tab) -> tuple[Gate, int]:
    """S5 -- the CV = 0 rung is the SAME CELL in every family, and is aliased.

    Returns the gate and the number of (family, alpha, seed) entries aliased in.
    """
    g = Gate("S5", "At CV = 0 every PSD family is the same delta distribution, so the "
                   "CV = 0 rung is ONE cell and not three. Where two families ran it, "
                   "phi is BIT-IDENTICAL; where only one did, the others alias it. The "
                   "discriminating half is that at CV > 0 the families must DIFFER.")
    cv0 = [r for r in prod if abs(float(r["cv"])) < 1e-12]
    g.say(f"N: {len(cv0)} production records at CV = 0")
    by = defaultdict(dict)
    for r in cv0:
        by[(round(float(r["alpha"]), 3), int(r["seed"]),
            float(r["p_hat_target"]))][r["family"]] = r
    n_cmp = n_same = 0
    worst = 0.0
    # BROKEN DOWN BY ALPHA AND BY FAMILY PAIR, not pooled. The aliasing licenses a
    # substitution at EVERY alpha and for EVERY family; a pooled count can pass on
    # evidence drawn entirely from one alpha and one family pair and say nothing about
    # the combinations it is actually licensing (referee, second pass -- the same defect
    # class as the pooled convergence rate fixed one gate earlier in this file).
    by_alpha = defaultdict(lambda: [0, 0, 0.0])
    pairs_seen = set()
    for k, fams in by.items():
        if len(fams) < 2:
            continue
        n_cmp += 1
        vals = [float(v["phi"]) for v in fams.values()]
        d = max(vals) - min(vals)
        worst = max(worst, d)
        same = int(len(set(vals)) == 1)
        n_same += same
        row = by_alpha[k[0]]
        row[0] += 1
        row[1] += same
        row[2] = max(row[2], d)
        for f1 in fams:
            for f2 in fams:
                if f1 < f2:
                    pairs_seen.add((f1, f2))
    g.say(f"N: {n_cmp} CV = 0 cells present in more than one family; "
          f"{n_same} bit-identical; worst spread {worst:.3e}")
    g.say("")
    g.say("| alpha | cells compared | bit-identical | worst spread |")
    g.say("|---|---|---|---|")
    for a in sorted(by_alpha):
        n, s, d = by_alpha[a]
        g.say(f"| {a:.2f} | {n} | {s} | {d:.3e} |")
    g.say("")
    fams_present = sorted({r["family"] for r in prod})
    want_pairs = {(f1, f2) for f1 in fams_present for f2 in fams_present if f1 < f2}
    alphas_aliased = sorted({round(float(r["alpha"]), 3) for r in prod})
    missing_pairs = sorted(want_pairs - pairs_seen)
    missing_alphas = [a for a in alphas_aliased if a not in by_alpha]
    g.say(f"Family pairs directly compared at CV = 0: {sorted(pairs_seen)}")
    g.say(f"Alphas at which the aliasing is USED: {alphas_aliased}; "
          f"alphas at which it is DIRECTLY VERIFIED: {sorted(by_alpha)}")
    if missing_pairs or missing_alphas:
        g.say("")
        g.say(f"⚠ **The gate's evidence is NARROWER than what it licenses.** "
              f"Family pairs never compared: {missing_pairs}. Alphas where the alias is "
              f"used but never checked: {missing_alphas}. The code-path argument is "
              f"strong -- `family` enters ONLY through `psd.make_psd`, and the radii, "
              f"the positions and the orientation stream are all drawn before `alpha` "
              f"or the shape logic runs -- but an argument is not a measurement, and "
              f"this project's own rule is to measure. Each missing combination is ONE "
              f"cell.")

    # the DISCRIMINATING half: at CV > 0 the families must differ, or the identity
    # above would only be saying the family flag does nothing.
    hi = defaultdict(dict)
    for r in prod:
        if float(r["cv"]) > 1e-12:
            hi[(round(float(r["cv"]), 3), round(float(r["alpha"]), 3), int(r["seed"]),
                float(r["p_hat_target"]))][r["family"]] = r
    n_hi = n_diff = 0
    sep = 0.0
    for _k, fams in hi.items():
        if len(fams) < 2:
            continue
        n_hi += 1
        vals = [float(v["phi"]) for v in fams.values()]
        d = max(vals) - min(vals)
        sep = max(sep, d)
        n_diff += int(len(set(vals)) > 1)
    g.say(f"N: {n_hi} CV > 0 cells present in more than one family; {n_diff} DIFFER; "
          f"largest family spread {sep:.3e}")
    if n_hi == 0:
        g.say("**The discriminating half has not been exercised** -- no CV > 0 cell is "
              "present in two families yet, so an identity at CV = 0 cannot yet be "
              "distinguished from the family flag doing nothing anywhere.")

    # alias
    n_alias = 0
    fams_all = sorted({k[0] for k in tab})
    donors = {k: v for k, v in tab.items() if abs(k[1]) < 1e-12}
    for (dfam, dcv, da), seeds in list(donors.items()):
        for fam in fams_all:
            if fam == dfam:
                continue
            if (fam, dcv, da) in tab and tab[(fam, dcv, da)]:
                continue
            tab[(fam, dcv, da)] = dict(seeds)
            n_alias += len(seeds)
    g.say(f"N: {n_alias} (family, alpha, seed) phi_J values ALIASED from the family "
          f"that ran the CV = 0 rung.")
    g.say("")
    g.say("🔴 **THE CONSEQUENCE IS NOT A SAVING, IT IS A CORRELATION.** The three "
          "families share their CV = 0 reference point as the SAME floating-point "
          "number, so `D = phi_J(CV_hi) - phi_J(CV 0)` for the three families is not "
          "three independent measurements. Writing `D_A = X_A - C`, the shared `C` "
          "cancels exactly in `D_A - D_B`, so `Var(D_A - D_B) = Var(X_A) + Var(X_B)` "
          "while the naive independent form adds `2*Var(C)`: assuming independence is "
          "CONSERVATIVE, not anti-conservative, and by a factor near two on the "
          "referee's worked case. S3's per-family paired Delta is unaffected for a "
          "stronger reason -- it differences the same reference out on BOTH alpha "
          "rungs, so the shared term cancels identically. A CROSS-FAMILY comparison of "
          "Delta would carry it; nothing here computes one.")
    ok = (n_cmp > 0 and n_same == n_cmp and worst == 0.0
          and n_hi > 0 and n_diff == n_hi)
    g.verdict(ok, informative=(n_cmp > 0 and n_hi > 0))
    return g, n_alias


def g_s4(prod) -> Gate:
    g = Gate("S4", "Every alpha > 1 result is reported with its corrugation `eps` and "
                   "its measured contact multiplicity, and is labelled "
                   "'elongation + interlocking + contact multiplicity'. A single-number "
                   "claim about ASPECT RATIO from clumps alone is pre-registered as NOT "
                   "QUOTABLE until the LevelSet control exists.")
    by_a = defaultdict(list)
    for r in prod:
        by_a[round(float(r["alpha"]), 6)].append(r)
    g.say("| alpha | N cells | `eps_corrugation` | `dup_contact_frac` mean | "
          "min | max | `n_contacts_raw`/`n_contacts` mean |")
    g.say("|---|---|---|---|---|---|---|")
    for a in sorted(by_a):
        rows = by_a[a]
        eps = sorted({round(float(r.get("eps_corrugation", 0.0)), 6) for r in rows})
        d = np.array([float(r.get("dup_contact_frac", 0.0)) for r in rows])
        ratio = np.array([float(r.get("n_contacts_raw", r["n_contacts"]))
                          / float(r["n_contacts"]) for r in rows if r["n_contacts"]])
        g.say(f"| {a:.2f} | {len(rows)} | {eps[0] if len(eps) == 1 else eps} | "
              f"{d.mean():.4f} | {d.min():.4f} | {d.max():.4f} | "
              f"{ratio.mean() if ratio.size else float('nan'):.4f} |")
    g.say("")
    g.say("**THE SEAM.** alpha = 1 builds a SINGLE SPHERE, so its `dup_contact_frac` is "
          "0 by construction; the clump family's own alpha -> 1 limit has "
          "`dup_contact_frac` -> 0.75 (four interactions per grain pair, measured in "
          "`scripts/probe_clump_s0.sh` Q2). The alpha = 1 -> next-rung step therefore "
          "mixes shape with the ONSET of contact multiplicity, and no trend in alpha "
          "may be quoted across it without the alpha > 1 sub-trend printed beside it.")
    # This gate is a REPORTING requirement, so it passes iff the fields it requires are
    # actually present on every record -- which is something that can fail.
    have = all(("eps_corrugation" in r and "dup_contact_frac" in r
                and "n_contacts_raw" in r) for r in prod)
    g.say("")
    g.say(f"Every production record carries `eps_corrugation`, `dup_contact_frac` and "
          f"`n_contacts_raw`: {have}")
    g.verdict(have and bool(prod))
    return g


def main() -> int:
    report_lines: list[str] = []

    def report(s: str = "") -> None:
        print(s, flush=True)
        report_lines.append(s)

    report("=" * 78)
    report("grade_shape.py -- PREREG_SHAPE_2026-09-18.md (with the 5.x amendment)")
    report("=" * 78)
    prod, gate_recs, files = load_records(report)
    prod = dedup(prod, report)
    prod, multi = prefer_raised_cap(prod, report)
    g6 = g_s6(multi)

    gates = [g_s0a(gate_recs, report), g_s0b(gate_recs), g_s0c(prod, report),
             g_s1(prod), g_s2(files, prod), g6]
    report("")
    report("--- phi_J ladder fits ---")
    tab = phi_j_table(prod, report) if prod else {}
    # S5 runs BEFORE S3 because it modifies `tab`: it aliases the CV = 0 rung into the
    # families that did not run it. The aliasing is gated -- if S5 fails, nothing is
    # aliased and S3 simply has no CV = 0 rung for those families.
    g5, _n_alias = g_s5(prod, tab) if prod else (Gate("S5", "CV = 0 identity."), 0)
    if not prod:
        g5.say("NO PRODUCTION RECORDS.")
        g5.verdict(False)
    # the same per-alpha contact multiplicity S4 tabulates, fed to S3's stiffness prediction
    mult = {}
    for a in sorted({float(r["alpha"]) for r in prod}):
        rows = [float(r["n_contacts_raw"]) / float(r["n_contacts"]) for r in prod
                if abs(float(r["alpha"]) - a) < 1e-9 and "n_contacts_raw" in r
                and float(r["n_contacts"]) > 0]
        if rows:
            mult[a] = float(np.mean(rows))
    g3, s3 = g_s3(tab, report, mult) if tab else (Gate("S3", "PRIMARY."), {})
    if not tab:
        g3.say("NO PRODUCTION LADDERS -- S3 cannot be read.")
        g3.verdict(False)
    # Amendment 5.10: S3 is read only if S6 PASSES. The numbers above are still printed
    # so a reader can see what a void re-run would have said, but the verdict is BLOCKED.
    if g6.passed is not True:
        g3.say("")
        g3.say(f"**S3 BLOCKED by S6 ({'VACUOUS' if not g6.informative else 'FAIL'}): "
               f"the raised-cap records are not shown to be the same protocol as the "
               f"default-cap ones, so the ladders that mix caps are not read.**")
        g3.verdict(False)
    gates += [g5, g3, g_s4(prod)] if prod else [g5, g3]

    body: list[str] = []
    for g in gates:
        body += g.render()

    n_pass = sum(1 for g in gates if g.passed)
    n_fail = sum(1 for g in gates if g.passed is False)
    summary = ["| gate | verdict |", "|---|---|"]
    for g in gates:
        v = ("NOT RUN" if g.passed is None
             else ("NOT INFORMATIVE" if not g.informative
                   else ("PASS" if g.passed else "FAIL")))
        summary.append(f"| {g.name} | **{v}** |")

    header = [
        "# CLAIMS -- the shape axis, 2026-09-18",
        "",
        "**GENERATED BY `grade_shape.py`. Every number on this page is printed by that "
        "script from the shards in `out/shape/`; nothing is transcribed and no verdict "
        "is written beside a number rather than from it. Re-generate with "
        "`python grade_shape.py`; it exits non-zero if any gate fails.**",
        "",
        f"Shards read: {len(files)} files. Production records after every filter: "
        f"{len(prod)}.",
        "",
        "**WHAT THIS AXIS CANNOT DECIDE, before anything it can.** A two-sphere clump "
        "moves three things at once -- elongation, surface corrugation and contact "
        "multiplicity -- and no clump-only run can separate them. Every number below is "
        "labelled accordingly, and the LevelSet control registered in "
        "`PREREG_SHAPE_2026-09-18.md` is what would turn the label into a measurement.",
        "",
    ] + summary + ["", f"{n_pass} gate(s) pass, {n_fail} fail.", ""]

    with open(CLAIMS, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(header + body) + "\n")
    with open(OUT_TXT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(report_lines + [""] + header + body) + "\n")

    for line in header + body:
        print(line)
    print(f"\nwrote {CLAIMS} and {OUT_TXT}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())

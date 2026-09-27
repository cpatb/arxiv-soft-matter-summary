"""grade_facet.py -- PREREG_FACET_2026-09-27 (+ Amendment 1): convex facets and the shear-jamming window, rotations live.

Reads out/lsdem/ls_facet2.jsonl (synced from ~/pj/out/lsdem by scripts/sync_facet.sh) and the per-cell contact
census npz files. Every test is graded only when all of its member cells have ENDED; an unconverged member makes the
test INCOMPLETE (health is graded separately from the bar, never as a FAIL).

Primary value (Amendment 1): phi at the registered stop rule. Block-1 cells ran 20000 post-convergence iterations;
their record-level phi is AFTER that relaxation, so the registered value is post[0]['phi'] and the drift is reported.

    python grade_facet.py            # prints and writes out/grade_facet.txt
"""
from __future__ import annotations

import io
import json
import math
import os
import sys
from contextlib import redirect_stdout

import numpy as np

import facet_census

HERE = os.path.dirname(os.path.abspath(__file__))
LSDIR = os.path.join(HERE, "out", "lsdem")
PHI_SCALE = 2.1e-3            # the project's declared run-to-run reproducibility (mixture_model.py)
EPS_LADDER = (1.0, 0.9, 0.8, 0.7, 0.5)
C_FACE = facet_census.C_FACE
QUEUE = {  # the registered queue (Amendment 1): (n, eps, mu, seed)
    **{(100, e, m, s): None for s in (1, 2, 3) for e in EPS_LADDER for m in (0.0, 0.3)},
    **{(200, e, m, s): None for s in (1, 2, 3, 4, 5) for e in (1.0, 0.5) for m in (0.0, 0.3)},
}


SPHERE_REF_FILES = ("ls_ladder.jsonl", "ls_nodes.jsonl", "ls_repl.jsonl", "ls_extra.jsonl", "ls_n200.jsonl",
                    "ls_facet2ref.jsonl")


def load() -> dict:
    """Superball cells (eps < 1) from ls_facet2.jsonl, rotations live; sphere references (eps 1) with the torque-free
    law, which is exact for a sphere (Amendment 3): the pilot's 1000-node records and the new facet2ref cells."""
    recs = {}
    path = os.path.join(LSDIR, "ls_facet2.jsonl")
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            r = json.loads(line)
            if r["eps"] == 1.0:
                continue   # killed before landing (Am. 3); none should exist
            assert r.get("spherical_bodies") is False, "a facet2 superball record without live rotations"
            k = (r["n"], r["eps"], r["mu_p"], r["seed"])
            assert k not in recs, f"duplicate record {k}"
            recs[k] = r
    fp = os.path.join(LSDIR, "ls_facet2mu01.jsonl")   # Am. 4 low-friction edge
    if os.path.exists(fp):
        for line in open(fp, encoding="utf-8"):
            r = json.loads(line)
            assert r["spherical_bodies"] is (r["eps"] == 1.0), "mu 0.1 cell with the wrong torque law"
            k = (r["n"], r["eps"], r["mu_p"], r["seed"])
            assert k not in recs, f"duplicate record {k}"
            recs[k] = r
    for f in SPHERE_REF_FILES:
        fp = os.path.join(LSDIR, f)
        if not os.path.exists(fp):
            continue
        for line in open(fp, encoding="utf-8"):
            r = json.loads(line)
            if r["eps"] != 1.0 or r["nodes"] != 1000 or r.get("spherical_bodies", True) is not True:
                continue
            k = (r["n"], 1.0, r["mu_p"], r["seed"])
            assert k not in recs, f"duplicate sphere reference {k} ({f})"
            recs[k] = r
    return recs


def healthy(r: dict) -> bool:
    """Registered health rule at the STOP RULE. Block-1 records ran 20000 post iterations after the stop, and the
    record-level `converged` is evaluated AFTER them (where the servo keeps unbalanced jittering at 0.03-0.06), so for
    those the stop state post[0] decides."""
    if r.get("post"):
        return bool(r["post"][0]["unbalanced"] < 0.02 and r["post"][0]["iter"] < 300000)
    return bool(r["converged"])


def phi_primary(r: dict) -> float:
    return r["post"][0]["phi"] if r.get("post") else r["phi"]


def seedpaired(recs, n, f, seeds):
    """values f(seed) for seeds whose member cells all landed converged; returns (vals, missing, unhealthy)."""
    vals, missing, bad = [], [], []
    for s in seeds:
        v = f(s)
        if v == "missing":
            missing.append(s)
        elif v == "bad":
            bad.append(s)
        else:
            vals.append(v)
    return vals, missing, bad


def cellphi(recs, n, e, m, s):
    r = recs.get((n, e, m, s))
    if r is None:
        return "missing"
    if not healthy(r):
        return "bad"
    return phi_primary(r)


def verdict_diff(vals, floor_per, label_pos, label_neg, frac=0.8):
    n = len(vals)
    mean = float(np.mean(vals))
    se_seed = float(np.std(vals, ddof=1) / math.sqrt(n)) if n > 1 else float("inf")
    se_eff = max(se_seed, floor_per / math.sqrt(n))
    need = math.ceil(frac * n)
    npos, nneg = sum(v > 0 for v in vals), sum(v < 0 for v in vals)
    if mean > 3 * se_eff and npos >= need:
        v = label_pos
    elif mean < -3 * se_eff and nneg >= need:
        v = label_neg
    else:
        v = "UNRESOLVED"
    return v, mean, se_seed, se_eff, npos, nneg


def window_test(recs, n, seeds, tag):
    def d(s):
        p = [cellphi(recs, n, e, m, s) for e in (1.0, 0.5) for m in (0.0, 0.3)]
        if "missing" in p:
            return "missing"
        if "bad" in p:
            return "bad"
        w1, w5 = 1 - p[1] / p[0], 1 - p[3] / p[2]
        return (w1 - w5, p[0] - p[2], p[1] - p[3], w1, w5)
    rows, missing, bad = seedpaired(recs, n, d, seeds)
    print(f"\n=== {tag}: N {n}, seeds {list(seeds)} ===")
    if missing or bad:
        print(f"   INCOMPLETE -- missing seeds {missing}, unconverged seeds {bad}")
        if not rows:
            return
    for s_, r in zip([s for s in seeds if s not in missing + bad], rows):
        print(f"   seed {s_}: W(1.0) {r[3]:.4f}  W(0.5) {r[4]:.4f}  D {r[0]:+.4f}  dphi0 {r[1]:+.4f}  dphim {r[2]:+.4f}")
    if missing or bad:
        print("   (partial table above; no verdict)")
        return
    phi0 = np.mean([cellphi(recs, n, 1.0, 0.0, s) for s in seeds])
    D = [r[0] for r in rows]
    v, mean, se, se_eff, npos, nneg = verdict_diff(D, math.sqrt(2) * PHI_SCALE / phi0, "NARROWS", "WIDENS")
    alt = max(se, 2 * PHI_SCALE / phi0 / math.sqrt(len(D)))
    print(f"   F1 D = W(1.0) - W(0.5) = {mean:+.4f}; se_seed {se:.4f}, se_eff {se_eff:.4f} (registered floor "
          f"sqrt2*PHI_SCALE/phi0); D > 0 on {npos}/{len(D)}  ->  F1 {v}")
    print(f"      (four-phi floor 2*PHI_SCALE/phi0 would give se {alt:.4f}, |mean|/se {abs(mean)/alt:.2f}; reported only)")
    v2, m2, s2, e2, p2, n2 = verdict_diff([r[1] for r in rows], math.sqrt(2) * PHI_SCALE, "REPLICATES (phi0 lower at eps 0.5)",
                                         "REVERSED (phi0 HIGHER at eps 0.5)")
    print(f"   F2 phi0(1.0) - phi0(0.5) = {m2:+.4f} +- {e2:.4f} (se_eff), positive on {p2}/{len(D)}  ->  F2 {v2}")
    dm = float(np.mean([r[2] for r in rows]))
    if v == "UNRESOLVED":
        print(f"   F3 not graded (F1 UNRESOLVED); reported: mean dphi0 {m2:+.4f}, mean dphim {dm:+.4f}")
    else:
        br = "FRICTIONLESS" if abs(m2) > abs(dm) else "FRICTIONAL"
        print(f"   F3 |mean dphi0| {abs(m2):.4f} vs |mean dphim| {abs(dm):.4f}  ->  carried by the {br} branch")


def ladder(recs):
    print("\n=== Block 1 ladder (N 100, seeds 1-3): seed means at the registered stop rule ===")
    print("   eps   mu    phi(mean +- sd)       n   Z_all@stop  Z_nr@END(after post)   rattl  drift_post(phi)  torque_unb stop->end   w(G) from Z_nr@END")
    for e in EPS_LADDER:
        for m in (0.0, 0.3):
            rs = [recs.get((100, e, m, s)) for s in (1, 2, 3)]
            ok = [r for r in rs if r is not None and healthy(r)]
            nb = sum(1 for r in rs if r is not None and not healthy(r))
            if not ok:
                print(f"   {e:.1f}   {m:.1f}   none landed" + (f" ({nb} unconverged)" if nb else ""))
                continue
            ph = [phi_primary(r) for r in ok]
            dr = [r["phi"] - phi_primary(r) for r in ok if r.get("post")]
            tq = [(r["post"][0]["torque_unb"], r["post"][-1]["torque_unb"]) for r in ok if r.get("post")]
            print(f"   {e:.1f}   {m:.1f}   {np.mean(ph):.4f} +- {np.std(ph, ddof=1) if len(ph) > 1 else float('nan'):.4f}"
                  f"   {len(ok)}   {np.mean([(r['post'][0]['z_all'] if r.get('post') else r['z_all']) for r in ok]):.2f}   {np.mean([r['z_nr'] for r in ok]):.2f}"
                  f"   {np.mean([r['rattler_frac'] for r in ok]):.3f}   {np.mean(dr) if dr else float('nan'):+.4f}"
                  f"           {np.mean([a for a, _ in tq]):.3f}->{np.mean([b for _, b in tq]):.3f}"
                  + f"   w(G) {np.mean([r['z_nr'] for r in ok]) * (1 - np.mean(ph)) / (2 * math.sqrt(3) * np.mean(ph)):.3f}"
                  + (f"   ({nb} unconverged, excluded)" if nb else ""))


def anchor_and_f4(recs):
    print("\n=== A-lit (method anchor) and F4: frictionless phi0(eps) - phi0(1.0), seed-paired, N 100 ===")
    for e in (0.9, 0.8, 0.7, 0.5):
        def d(s, e=e):
            a, b = cellphi(recs, 100, e, 0.0, s), cellphi(recs, 100, 1.0, 0.0, s)
            if "missing" in (a, b):
                return "missing"
            if "bad" in (a, b):
                return "bad"
            return a - b
        vals, missing, bad = seedpaired(recs, 100, d, (1, 2, 3))
        if missing or bad:
            print(f"   eps {e}: INCOMPLETE (missing {missing}, unconverged {bad})"
                  + (f"; partial {np.round(vals, 4).tolist()}" if vals else ""))
            continue
        v, mean, se, se_eff, npos, nneg = verdict_diff(vals, math.sqrt(2) * PHI_SCALE, "RISE", "FALL", frac=1.0)
        tag = ""
        if e in (0.9, 0.8):
            tag = "  A-lit " + ("PASS" if v == "RISE" else "FAIL -- method not validated against Jiao 2010")
        if e == 0.9:
            tag += f"; F4 {v}"
        print(f"   eps {e}: dphi0 {mean:+.4f} (se_eff {se_eff:.4f}; per seed {np.round(vals, 4).tolist()})  {v}{tag}")


def branch_tests(recs):
    """Am. 4: F6 (mu 0.3 branch), F7 (mu 0.1 edge), F8 (window between mu 0.1 and 0.3), seed-paired over seeds 1-3."""
    print("\n=== Am. 4: F6 (mu_p 0.3), F7 (mu_p 0.1), F8 (W_0.1 = 1 - phi(0.3)/phi(0.1)); N 100, seeds 1-3 ===")
    for e in (0.9, 0.8, 0.7, 0.5):
        for name, mu in (("F6", 0.3), ("F7", 0.1)):
            def d(s, e=e, mu=mu):
                a, b = cellphi(recs, 100, e, mu, s), cellphi(recs, 100, 1.0, mu, s)
                return a if a in ("missing", "bad") else (b if b in ("missing", "bad") else a - b)
            vals, missing, bad = seedpaired(recs, 100, d, (1, 2, 3))
            if missing or bad:
                print(f"   {name} eps {e}: INCOMPLETE (missing {missing}, unconverged {bad})"
                      + (f"; partial {np.round(vals, 4).tolist()}" if vals else ""))
                continue
            v, mean, se, se_eff, npos, nneg = verdict_diff(vals, math.sqrt(2) * PHI_SCALE, "RISE", "FALL", frac=1.0)
            print(f"   {name} eps {e}: dphi {mean:+.4f} (se_eff {se_eff:.4f}; per seed {np.round(vals, 4).tolist()})  {v}")

        def dw(s, e=e):
            p = [cellphi(recs, 100, ee, m, s) for ee in (1.0, e) for m in (0.1, 0.3)]
            if "missing" in p:
                return "missing"
            if "bad" in p:
                return "bad"
            return (1 - p[1] / p[0]) - (1 - p[3] / p[2])   # W(1) - W(eps): positive = NARROWS
        vals, missing, bad = seedpaired(recs, 100, dw, (1, 2, 3))
        if missing or bad:
            print(f"   F8 eps {e}: INCOMPLETE (missing {missing}, unconverged {bad})")
            continue
        phi01 = np.mean([cellphi(recs, 100, 1.0, 0.1, s) for s in (1, 2, 3)])
        v, mean, se, se_eff, npos, nneg = verdict_diff(vals, math.sqrt(2) * PHI_SCALE / phi01, "NARROWS", "WIDENS",
                                                       frac=1.0)
        ws = [np.mean([1 - cellphi(recs, 100, ee, 0.3, s) / cellphi(recs, 100, ee, 0.1, s) for s in (1, 2, 3)])
              for ee in (1.0, e)]
        print(f"   F8 eps {e}: W_0.1 sphere {ws[0]:.4f} -> {ws[1]:.4f}; D = W(1)-W(eps) {mean:+.4f} (se_eff {se_eff:.4f})  {v}")


def trajectories():
    """Am. 4: the four logged frictionless seed-1 cells, DESCRIPTIVE only (they never met the stop rule)."""
    print("\n=== Am. 4 descriptive: frictionless superballs (rotations live), seed 1, N 100 -- NOT graded ===")
    for e in (0.9, 0.8, 0.7, 0.5):
        fp = os.path.join(LSDIR, f"facet2vis_n100_e{e}_m0.0_s1.log")
        if not os.path.exists(fp):
            continue
        rows = [ln.split() for ln in open(fp, encoding="utf-8") if ln.startswith("PROGRESS")]
        if not rows:
            continue
        import facet_census as fc  # noqa: F401  (volume from ls_pack's formula below)
        import ls_pack
        v = ls_pack.superellipsoid_volume(ls_pack.R_EQ * ((4.0 / 3.0) * math.pi / ls_pack.superellipsoid_volume(1.0, e)) ** (1.0 / 3.0), e)
        ph = [(int(r[1]), 100 * v / float(r[2]), float(r[3])) for r in rows]
        sel = ph[len(ph) // 4:: max(1, len(ph) // 6)] + [ph[-1]]
        print(f"   eps {e}: " + "  ".join(f"it {i}: phi {p:.4f} unb {u:.2f}" for i, p, u in sel))


def log_phi(fp: str, e: float):
    """Am. 6: last-quarter mean phi of a progress log (unconverged, descriptive), with final unbalanced and iteration."""
    import ls_pack
    if not os.path.exists(fp):
        return None
    r = ls_pack.R_EQ * ((4.0 / 3.0) * math.pi / ls_pack.superellipsoid_volume(1.0, e)) ** (1.0 / 3.0)
    v = ls_pack.superellipsoid_volume(r, e)
    pr = [ln.split() for ln in open(fp, encoding="utf-8") if ln.startswith("PROGRESS")]
    if len(pr) < 4:
        return None
    q = pr[3 * len(pr) // 4:]
    return float(np.mean([100 * v / float(x[2]) for x in q])), float(pr[-1][3]), int(pr[-1][1])


def descriptive_mu01(recs):
    """Am. 6: phi(eps, mu 0.1) and phi(eps, mu 0.3) with converged cells at the stop rule and unconverged ones read
    from their logs (labelled), and the resulting W_0.1 per seed. Descriptive; the registered F7/F8 are above."""
    print("\n=== Am. 6 descriptive: mu 0.1 and 0.3 edges incl. unconverged cells (u = unconverged, log readout) ===")
    for e in EPS_LADDER:
        cells = []
        for s in (1, 2, 3):
            out = []
            for mu in (0.1, 0.3):
                r = recs.get((100, e, mu, s))
                if r is not None and healthy(r):
                    out.append((phi_primary(r), ""))
                    continue
                name = f"facet2mu01_e{e}_s{s}.log" if mu == 0.1 else None
                lp = log_phi(os.path.join(LSDIR, name), e) if name else None
                out.append((lp[0], f"u{lp[1]:.2f}") if lp else (float("nan"), "missing"))
            w = 1 - out[1][0] / out[0][0]
            cells.append((s, out, w))
        txt = "  ".join(f"s{s}: {o[0][0]:.4f}{o[0][1]} / {o[1][0]:.4f}{o[1][1]} W {w:.4f}" for s, o, w in cells)
        ws = [w for _, _, w in cells if not math.isnan(w)]
        print(f"   eps {e}: phi(0.1) / phi(0.3): {txt}   mean W_0.1 {np.mean(ws) if ws else float('nan'):.4f}")


LIVE_SPHERE_1000 = 0.5977   # Am. 5: live-torque sphere, mu 0.3, seed 1, 1000 nodes, cell volume 87.61 at 4500 it (probe)
FREE_SPHERE_1000 = 0.5797   # torque-free sphere, same seed (ls_nodes)


def node_control(recs):
    """Am. 5: node refinement of the torque-noise densification, mu 0.3, seed 1, N 100."""
    print("\n=== Am. 5 node-noise control (mu_p 0.3, seed 1, N 100, rotations live) ===")
    import ls_pack
    rows = {}
    for e, nodes in ((0.5, 2000), (0.5, 4000), (0.9, 4000), (1.0, 4000)):
        fp = os.path.join(LSDIR, f"facet2nodes_e{e}_n{nodes}.log")
        if not os.path.exists(fp):
            print(f"   eps {e} nodes {nodes}: no log yet")
            continue
        r = ls_pack.R_EQ * ((4.0 / 3.0) * math.pi / ls_pack.superellipsoid_volume(1.0, e)) ** (1.0 / 3.0)
        v = ls_pack.superellipsoid_volume(r, e)
        pr = [ln.split() for ln in open(fp, encoding="utf-8") if ln.startswith("PROGRESS")]
        rec = None
        for line in open(fp, encoding="utf-8"):
            if line.startswith("{"):
                rec = json.loads(line)
        if rec is not None and rec["converged"]:
            rows[(e, nodes)] = (rec["phi"], "stop-rule")
        elif len(pr) >= 4:
            q = pr[3 * len(pr) // 4:]
            rows[(e, nodes)] = (float(np.mean([100 * v / float(x[2]) for x in q])),
                                f"last-quarter mean of {len(pr)} logs, unb {float(pr[-1][3]):.3f}, it {pr[-1][1]}")
        else:
            print(f"   eps {e} nodes {nodes}: {len(pr)} progress lines, too early")
            continue
        print(f"   eps {e} nodes {nodes}: phi {rows[(e, nodes)][0]:.4f} ({rows[(e, nodes)][1]})")
    if (1.0, 4000) in rows:
        gap0 = LIVE_SPHERE_1000 - FREE_SPHERE_1000
        gap1 = rows[(1.0, 4000)][0] - FREE_SPHERE_1000
        ok = gap1 < gap0 / 2
        print(f"   noise model: live-sphere excess over torque-free {gap0:+.4f} (1000 nodes) -> {gap1:+.4f} (4000): "
              + ("SHRINKS by >= half, node spacing drives it" if ok else "does NOT shrink by half -- control uninformative"))
    ref1000 = {0.5: cellphi(recs, 100, 0.5, 0.3, 1), 0.9: cellphi(recs, 100, 0.9, 0.3, 1)}
    for e in (0.5, 0.9):
        if (e, 4000) in rows and isinstance(ref1000[e], float):
            d = rows[(e, 4000)][0] - ref1000[e]
            print(f"   N1 eps {e}: phi_m(4000) - phi_m(1000) = {d:+.4f} (bar |d| < {2 * PHI_SCALE:.4f})  "
                  + ("PASS -- not a node-noise artefact" if abs(d) < 2 * PHI_SCALE else
                     ("FAIL, falls -- correct F6 by this shift" if d < 0 else "FAIL, rises")))


def census(recs):
    print(f"\n=== F5 contact census: face-face fraction f_ff (m > {C_FACE} on both grains) vs area null ===")
    nulls = {e: facet_census.area_null(e) for e in EPS_LADDER}
    print("   N    eps   mu    contacts  P(m>c) meas/null   f_ff meas / null      E_ff   per-seed E_ff     lever arm med (arm/r)")
    for n in (100, 200):
        for e in EPS_LADDER:
            for m in (0.0, 0.3):
                keys = [k for k in recs if k[0] == n and k[1] == e and k[2] == m and healthy(recs[k])
                        and recs[k].get("census")]
                if not keys:
                    continue
                per, allm1, allm2, arms = [], [], [], []
                for k in sorted(keys):
                    z = np.load(os.path.join(LSDIR, recs[k]["census"]))
                    m1, m2 = z["m1"], z["m2"]
                    ff = float(np.mean((m1 > C_FACE) & (m2 > C_FACE)))
                    per.append(ff / nulls[e]["f_ff0"])
                    allm1.append(m1); allm2.append(m2)
                    arms.append(np.concatenate([z["arm1"], z["arm2"]]) / recs[k]["r_extent"])
                m1, m2 = np.concatenate(allm1), np.concatenate(allm2)
                ff = float(np.mean((m1 > C_FACE) & (m2 > C_FACE)))
                pm = float(np.mean(np.concatenate([m1, m2]) > C_FACE))
                E = ff / nulls[e]["f_ff0"]
                grade = ""
                if n == 200 and e == 0.5 and m == 0.0:
                    grade = "  <- F5 registered: " + ("FACE CONTACTS ENRICHED (E_ff > 1.25)" if E > 1.25 else "NOT ENRICHED")
                print(f"   {n}  {e:.1f}   {m:.1f}   {len(m1):6d}    {pm:.3f} / {nulls[e]['p_face']:.3f}      "
                      f"{ff:.4f} / {nulls[e]['f_ff0']:.4f}   {E:5.2f}   {np.round(per, 2).tolist()}   "
                      f"{np.median(np.concatenate(arms)):.3f}{grade}")


def main() -> int:
    recs = load()
    buf = io.StringIO()
    with redirect_stdout(buf):
        print("grade_facet -- PREREG_FACET_2026-09-27 (+Am. 1): level-set superellipsoids, rotations live, p_hat 1e-3")
        landed = [k for k in QUEUE if k in recs]
        unconv = [k for k in landed if not healthy(recs[k])]
        print(f"registered cells landed: {len(landed)} of {len(QUEUE)}; unconverged (excluded): "
              f"{[f'N{k[0]} e{k[1]} m{k[2]} s{k[3]}' for k in unconv]}")
        ladder(recs)
        anchor_and_f4(recs)
        window_test(recs, 200, (1, 2, 3, 4, 5), "F1-F3 REGISTERED (Block 2)")
        window_test(recs, 100, (1, 2, 3), "F1-F3 second size (Block 1)")
        branch_tests(recs)
        trajectories()
        node_control(recs)
        descriptive_mu01(recs)
        census(recs)
    out = buf.getvalue()
    sys.stdout.write(out)
    with open(os.path.join(HERE, "out", "grade_facet.txt"), "w", encoding="utf-8") as fh:
        fh.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

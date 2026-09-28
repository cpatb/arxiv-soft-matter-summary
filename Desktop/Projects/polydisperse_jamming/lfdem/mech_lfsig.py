"""mech_lfsig.py -- POST HOC mechanism readout (not in PREREG_LFSIG2): why does the dry-DEM chain's geometric
factor survive lubrication while the activation stress does not move?

    python lfdem/mech_lfsig.py [lfsig|lfsig2]   ->  out/mech_<campaign>.txt

Per cell and per interaction snapshot with strain in [2, 4] (int_*.dat, written every 1 strain), in LF_DEM
repulsion units (F_R0 = 1; the repulsion on a pair at contact is the reduced radius R* = a_i a_j/(a_i + a_j),
RepulsiveForce.cpp:32, so the activation THRESHOLD of a pair is R*):

  touching contacts   : state (col 3) >= 1, i.e. h < 0.  Normal load F_n = spring (col 12) + repulsion (col 16).
  g = F_n / R*        : load in units of the pair's own threshold (>= 1 by construction for a touching pair).
  N_c/V, Z_c = 2 N_c/N, <a_i a_j>_c                    -- the chain's geometric factor G = (N_c/V) <a_i a_j>_c
  sigma_c (identity)  : (1/V) sum F_n |l| n_x n_z over touching contacts, |l| = a_i + a_j + h (h < 0),
                        the normal part of the Love-Weber contact shear stress; compared with LF_DEM's own
                        Viscosity(contact) x shear rate from data_*.dat at the same strain (the identity check).
  near pairs          : state 0 with gap h < lambda = 0.02 (in units of R*... LF_DEM's col 7 is dimensionless
                        gap = s - 2, s = 2r/(a_i + a_j), i.e. h/(a_i + a_j)/2... see NOTE), their normal load
                        (lubrication col 8 + repulsion col 16) over the threshold R*: how close the not-yet-
                        touching population sits to activation.

NOTE on the gap: LF_DEM's dimensionless gap is 2 r/(a_i+a_j) - 2, so the physical gap is h = gap (a_i+a_j)/2.

Every number here is computed from the dumped files; nothing is typed. The chain's prediction for the RATIO
CV 0.5 / CV 0.2 of G at matched sigma is what the Yade sheared packs gave (0.768 for sigma*); this script
prints the LF_DEM ratios of G, <g>, sigma_c and the contact share so the discrepancy can be located.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 2.0, 4.0
LAMBDA = 0.02   # repulsive_length, params_base.txt (in units of the number-mean radius a = 1)
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def snapshots(path: str):
    strain, rate, rows = None, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# cumulated strain"):
                if strain is not None and rows:
                    yield strain, rate, np.array(rows)
                strain, rows = float(line.split(":")[1]), []
            elif line.startswith("# shear rate"):
                rate = float(line.split(":")[1])
            elif line.startswith("#") or not line.strip():
                continue
            elif strain is not None:
                rows.append([float(x) for x in line.split()[:17]])
    if strain is not None and rows:
        yield strain, rate, np.array(rows)


def cell(d: str) -> dict | None:
    ints = glob.glob(os.path.join(d, "int_*.dat"))
    data = glob.glob(os.path.join(d, "data_*.dat"))
    if not ints or not data or not os.path.exists(os.path.join(d, "conf.dat")):
        return None
    kv = dict(re.findall(r"(\w+)=(\S+)", open(os.path.join(d, "META"), encoding="utf-8").read()))
    conf = np.loadtxt(os.path.join(d, "conf.dat"), comments="#", ndmin=2)
    r = conf[:, 3]
    hdr = [l for l in open(os.path.join(d, "conf.dat"), encoding="utf-8") if l.startswith("#")][1].split()
    L = float(hdr[4]); V = L**3; N = r.size
    a = np.loadtxt(data[0], comments="#", ndmin=2)
    rows = []
    for strain, rate, s in snapshots(ints[0]):
        if not (W0 <= strain <= W1):
            continue
        i, j = s[:, 0].astype(int), s[:, 1].astype(int)
        ai, aj = r[i], r[j]
        rstar = ai * aj / (ai + aj)
        touching = s[:, 2] >= 1
        gap = s[:, 6]                      # dimensionless: 2 rr/(ai+aj) - 2
        h = gap * (ai + aj) / 2.0          # physical gap (negative for touching)
        fn = s[:, 11] + s[:, 15]           # spring + repulsion at the pair's separation
        t = touching
        nc = int(t.sum())
        if nc == 0:
            continue
        g = fn[t] / rstar[t]
        ell = (ai + aj + h)[t]
        # Love-Weber: sigma_ab = (1/V) sum l_a F_b with l = ell n (i -> j) and F the force on j (normal F_n n plus the
        # tangential vector, cols 13-15). LF_DEM's printed stress is in units F_R/(6 pi a^2) (SimulationInit.cpp), so
        # the dump's F_R-unit sum is multiplied by 6 pi; the sign convention is fixed by the component that matches.
        F = fn[t, None] * s[t, 3:6] + s[t, 12:15]
        lvec = ell[:, None] * s[t, 3:6]
        T = 6.0 * np.pi * np.einsum("ka,kb->ab", lvec, F) / V
        Tn = 6.0 * np.pi * np.einsum("ka,kb->ab", lvec, fn[t, None] * s[t, 3:6]) / V
        Tf = 6.0 * np.pi * np.einsum("ka,kb->ab", lvec, fn[t, None] * s[t, 3:6] - s[t, 12:15]) / V   # tangential sign flipped
        sig_c = -float(T[0, 2]); sig_c_n = -float(Tn[0, 2]); sig_c_zx = -float(Tf[0, 2])
        # LF_DEM's own contact stress at this strain: Viscosity(contact) * shear rate, nearest data row
        k = int(np.argmin(np.abs(a[:, 1] - strain)))
        sig_c_lf = float(a[k, 4] * a[k, 2]); sig_tot_lf = float(a[k, 3] * a[k, 2])
        near = (~t) & (h < LAMBDA) & (h >= 0)
        gnear = (s[near, 7] + s[near, 15]) / rstar[near] if near.any() else np.array([np.nan])
        rows.append({"strain": strain, "rate": rate, "nc": nc, "zc": 2.0 * nc / N, "ncV": nc / V,
                     "aa": float(np.mean(ai[t] * aj[t])), "G": nc / V * float(np.mean(ai[t] * aj[t])),
                     "g_mean": float(g.mean()), "g_med": float(np.median(g)), "g_q90": float(np.quantile(g, 0.9)),
                     "g_max": float(g.max()), "sig_c": sig_c, "sig_c_n": sig_c_n, "sig_c_zx": sig_c_zx, "sig_c_lf": sig_c_lf, "sig_tot_lf": sig_tot_lf,
                     "n_near": int(near.sum()), "gnear_mean": float(np.nanmean(gnear)), "gnear_q90": float(np.nanquantile(gnear, 0.9)),
                     "a3": float(np.mean(r**3)), "phi": float(np.sum(4 / 3 * np.pi * r**3) / V)})
    if not rows:
        return None
    m = {k: float(np.mean([x[k] for x in rows])) for k in rows[0] if k != "strain"}
    m.update({"label": kv["LABEL"], "cv": float(kv["CV"]), "sig": float(kv["SIG"]), "seed": int(kv.get("SEED", 1)),
              "nsnap": len(rows), "N": N})
    return m


def grade_m(cells: list) -> None:
    """Amendment 1, M1-M3 (descriptive, registered before any lfsig2 cell was read). Needs every arm complete at
    sigma_0 3, 10, 15, 30, 100 (M1, M2) and the grader's sigma_half (M3, imported from grade_lfsig2)."""
    import grade_lfsig2 as G
    say("")
    say("==== Amendment 1: M1-M3 on the lfsig2 dumps (3 seeds) ====")
    by = {(c["cv"], c["seed"], c["sig"]): c for c in cells}
    seeds = (1, 2, 3)
    # M1: <aa>_c(3)/<aa>_c(100) > 1 at CV 0.5 on 3/3; excess larger at CV 0.5 than at CV 0.2 on 3/3
    m1a = m1b = n1 = 0
    for sd in seeds:
        if all((cv, sd, s) in by for cv in (0.2, 0.5) for s in (3.0, 100.0)):
            n1 += 1
            e5 = by[(0.5, sd, 3.0)]["aa"] / by[(0.5, sd, 100.0)]["aa"]
            e2 = by[(0.2, sd, 3.0)]["aa"] / by[(0.2, sd, 100.0)]["aa"]
            m1a += int(e5 > 1); m1b += int(e5 > e2)
            say(f"  M1 seed {sd}: <aa>_c(3)/<aa>_c(100) CV 0.2 {e2:.3f}, CV 0.5 {e5:.3f}")
    say(f"  M1: {'PASS' if n1 == 3 and m1a == 3 and m1b == 3 else 'NOT GRADED' if n1 < 3 else 'FAIL'} (excess > 1 at CV 0.5: {m1a}/{n1}; larger than CV 0.2: {m1b}/{n1})")
    # M2: G ratio > 1 at sigma 3; < 0.90 at sigma >= 10, on 3/3 seeds each
    say("  M2: G(0.5)/G(0.2) per seed at each sigma")
    ok_lo = ok_hi = n2 = 0
    for sd in seeds:
        sigs = sorted({s for (cv, d, s) in by if d == sd and cv in (0.2, 0.5)})
        row = []
        for s in sigs:
            if (0.2, sd, s) in by and (0.5, sd, s) in by:
                row.append((s, by[(0.5, sd, s)]["G"] / by[(0.2, sd, s)]["G"]))
        if not row:
            continue
        say(f"    seed {sd}: " + " ".join(f"{s:g}:{r:.3f}" for s, r in row))
        d = dict(row)
        if 3.0 in d and all(s in d for s in (10.0, 15.0, 30.0, 100.0)):
            n2 += 1; ok_lo += int(d[3.0] > 1.0); ok_hi += int(all(d[s] < 0.90 for s in (10.0, 15.0, 30.0, 100.0)))
    say(f"  M2: {'PASS' if n2 == 3 and ok_lo == 3 and ok_hi == 3 else 'NOT GRADED' if n2 < 3 else 'FAIL'} (> 1 at 3: {ok_lo}/{n2}; < 0.90 at every sigma >= 10: {ok_hi}/{n2})")
    # M3: G ratio interpolated (ln G vs ln sigma) at each arm's own sigma_half, within +-0.15 of 1
    import glob as _g
    dirs = [x for x in sorted(_g.glob(os.path.join(OUT, "lfsig2", "c*"))) if ".killed_" not in x]
    gc = {(c["cv"], c["seed"], c["sig"]): c for c in (G.load_cell(d, False) for d in dirs) if c and c["healthy"]}
    vals = []
    for sd in seeds:
        r2 = G.arm_readouts(gc, 0.2, sd, G.SIGMAS); r5 = G.arm_readouts(gc, 0.5, sd, G.SIGMAS)
        if not (r2["complete"] and r5["complete"] and np.isfinite(r2["shalf"]) and np.isfinite(r5["shalf"])):
            continue
        def g_at(cv, s0):
            pts = sorted((s, by[(cv, sd, s)]["G"]) for (c, d, s) in by if c == cv and d == sd)
            xs = np.log([p[0] for p in pts]); ys = np.log([p[1] for p in pts])
            return float(np.exp(np.interp(np.log(s0), xs, ys)))
        ratio = g_at(0.5, r5["shalf"]) / g_at(0.2, r2["shalf"])
        vals.append(ratio)
        say(f"  M3 seed {sd}: sigma_half 0.2 {r2['shalf']:.3f}, 0.5 {r5['shalf']:.3f}; G(0.5 at its half)/G(0.2 at its half) = {ratio:.3f}")
    if len(vals) == 3:
        v = np.array(vals); m, se = v.mean(), v.std(ddof=1) / np.sqrt(3)
        say(f"  M3: mean {m:.3f} se {se:.3f} band [{m - se:.3f}, {m + se:.3f}]; within +-0.15 of 1: {'PASS' if abs(m - 1) <= 0.15 else 'FAIL'}")
    else:
        say(f"  M3: NOT GRADED ({len(vals)}/3 seeds)")


def main() -> int:
    camp = sys.argv[1] if len(sys.argv) > 1 else "lfsig"
    say(f"mech_lfsig -- POST HOC, not registered: the chain's geometry vs the activation stress in LF_DEM [{camp}], snapshots strain [{W0}, {W1}]")
    say("  identity check: sigma_c from the dump (normal + tangential, 6 pi) vs LF_DEM's Viscosity(contact) x rate at the nearest data row")
    say("  columns: sig_c = tangential as written; (norm) = normal part only; (t+) = tangential read as the force on particle 1 -- the")
    say("  (t+) column is the one that matches (0.88-1.02 on lfsig; the residual is the ~0.025-strain offset between dump and data row)")
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, camp, "c*"))) if ".killed_" not in x and "_x95_" not in x]
    cells = [c for c in (cell(d) for d in dirs) if c]
    say(f"  {len(cells)} cells with interaction dumps in the window")
    say(f"  {'label':18s} {'sig':>5s} sn {'Zc':>6s} {'Nc/V':>7s} {'<aa>_c':>7s} {'G':>7s} {'<g>':>6s} {'g_med':>6s} {'g_q90':>6s} {'sig_c':>7s} {'(norm)':>7s} {'(t+)':>7s} {'sig_cLF':>7s} {'ratio':>6s} {'fcont':>6s} {'n_near':>6s} {'gnear':>6s} {'gn_q90':>6s}")
    for c in sorted(cells, key=lambda c: (c["cv"], c["seed"], c["sig"])):
        say(f"  {c['label']:18s} {c['sig']:5.1f} {c['nsnap']:2d} {c['zc']:6.3f} {c['ncV']:7.4f} {c['aa']:7.4f} {c['G']:7.4f} "
            f"{c['g_mean']:6.3f} {c['g_med']:6.3f} {c['g_q90']:6.2f} {c['sig_c']:7.3f} {c['sig_c_n']:7.3f} {c['sig_c_zx']:7.3f} {c['sig_c_lf']:7.3f} {c['sig_c_zx'] / c['sig_c_lf']:6.3f} "
            f"{c['sig_c_lf'] / c['sig_tot_lf']:6.3f} {c['n_near']:6.0f} {c['gnear_mean']:6.3f} {c['gnear_q90']:6.3f}")
    # ratios CV 0.5 / CV 0.2 at matched (sigma, seed)
    say("\n  CV 0.5 / CV 0.2 at matched (sigma, seed): the chain's factor G = (Nc/V)<aa>_c, the mean load per threshold <g>,")
    say("  the contact stress sigma_c and the contact share fcont.  (chain, dry Yade: sigma* ratio 0.768 <=> G ratio at fixed C/q)")
    say(f"  {'sig':>5s} {'seed':>4s} {'Zc':>6s} {'Nc/V':>6s} {'<aa>_c':>6s} {'G':>6s} {'<g>':>6s} {'g_med':>6s} {'sig_c':>6s} {'fcont':>6s} {'G*<g>':>6s}")
    by = {(c["cv"], c["seed"], c["sig"]): c for c in cells}
    for (cv, sd, sig), c5 in sorted(by.items()):
        if cv != 0.5 or (0.2, sd, sig) not in by:
            continue
        c2 = by[(0.2, sd, sig)]
        rr = lambda k: c5[k] / c2[k]  # noqa: E731
        say(f"  {sig:5.1f} {sd:4d} {rr('zc'):6.3f} {rr('ncV'):6.3f} {rr('aa'):6.3f} {rr('G'):6.3f} {rr('g_mean'):6.3f} {rr('g_med'):6.3f} "
            f"{rr('sig_c'):6.3f} {c5['sig_c_lf'] / c5['sig_tot_lf'] / (c2['sig_c_lf'] / c2['sig_tot_lf']):6.3f} {rr('G') * rr('g_mean'):6.3f}")
    if camp == "lfsig2":
        grade_m(cells)
    outp = os.path.join(OUT, f"mech_{camp}.txt")
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

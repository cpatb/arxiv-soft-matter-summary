"""grade_lspilot.py -- PREREG_LSPILOT_2026-09-26.md (+ Amendments 1-2). Writes out/grade_lspilot.txt.

Inputs: out/lsdem/ls_{pilot,nodes,ladder}.jsonl (synced by scripts/ls_*.sh), out/lsdem/phase2_lsref_lognormal.jsonl
(FrictMat N 100 references). Written before the 1000-node ladder landed. A PILOT: one seed, N 100, p_hat 1e-3.
"""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "out", "lsdem")
SQ3 = 2 * math.sqrt(3)
L = []


def say(s=""):
    print(s)
    L.append(s)


def load(name):
    p = os.path.join(D, name)
    return [json.loads(x) for x in open(p, encoding="utf-8") if x.strip()] if os.path.exists(p) else []


def main():
    say("grade_lspilot -- PREREG_LSPILOT_2026-09-26.md (pilot: seed 1, N 100, p_hat 1e-3)")
    ref = {r["mu_p"]: r for r in load("phase2_lsref_lognormal.jsonl") if r.get("converged")}
    say(f"FrictMat N100 references: " + ", ".join(f"mu {m}: phi {r['phi']:.4f} Z {r['z_nr']:.2f}" for m, r in sorted(ref.items())))
    cells = {}
    for r in load("ls_nodes.jsonl") + load("ls_ladder.jsonl"):
        if r["nodes"] == 1000:
            cells[(r["eps"], r["mu_p"])] = r
    say(f"1000-node cells landed: {len(cells)} of 8 registered")
    for k, r in sorted(cells.items()):
        if not r["converged"]:
            say(f"   NOT CONVERGED eps {k[0]} mu {k[1]}: unbalanced {r['unbalanced']:.3f} iters {r['iters']} -- reported, not graded")
    ok = {k: r for k, r in cells.items() if r["converged"]}
    say("\n=== P0 at 1000 nodes: level-set sphere vs FrictMat (bar 0.01) ===")
    p0 = []
    for m in (0.0, 0.3):
        r, f = ok.get((1.0, m)), ref.get(m)
        if r and f:
            d = r["phi"] - f["phi"]
            p0.append(abs(d) < 0.01)
            say(f"   mu {m}: LS {r['phi']:.4f} vs FrictMat {f['phi']:.4f}  d {d:+.4f}  {'ok' if p0[-1] else 'OUT'}")
        else:
            p0.append(None)
            say(f"   mu {m}: INCOMPLETE")
    say(f"P0: {'INCOMPLETE' if None in p0 else ('PASS' if all(p0) else 'FAIL')}")
    say("\n=== the ladder: phi_0, phi_m, window W = 1 - phi_m/phi_0, and the (G) decomposition ===")
    say("   eps    phi_0   Z0    w0     | phi_m   Zm    wm     | Zm/Z0   W")
    W = {}
    for e in (1.0, 0.7, 0.5, 0.35):
        a, b = ok.get((e, 0.0)), ok.get((e, 0.3))
        if not (a and b):
            say(f"   {e:<5}  INCOMPLETE ({'mu0 ' if not a else ''}{'mu.3' if not b else ''} missing)")
            continue
        w0 = a["z_nr"] * (1 - a["phi"]) / (SQ3 * a["phi"])
        wm = b["z_nr"] * (1 - b["phi"]) / (SQ3 * b["phi"])
        W[e] = 1 - b["phi"] / a["phi"]
        say(f"   {e:<5}  {a['phi']:.4f}  {a['z_nr']:.2f}  {w0:.3f}  | {b['phi']:.4f}  {b['z_nr']:.2f}  {wm:.3f}  | "
            f"{b['z_nr']/a['z_nr']:.3f}   {W[e]:.4f}")
    say("\n=== P1: W(eps < 1) > W(1) ===")
    if 1.0 in W:
        for e in (0.7, 0.5, 0.35):
            if e in W:
                say(f"   eps {e}: W {W[e]:.4f} vs sphere {W[1.0]:.4f}  ({'wider' if W[e] > W[1.0] else 'NOT wider'}; x{W[e]/W[1.0]:.2f})")
        rest = [W[e] > W[1.0] for e in (0.7, 0.5, 0.35) if e in W]
        say(f"P1: {'INCOMPLETE' if len(rest) < 3 else ('PASS' if all(rest) else 'FAIL')} ({sum(rest)}/{len(rest)}) -- pilot, one seed")
    else:
        say("P1: INCOMPLETE (no eps 1 pair)")
    with open(os.path.join(HERE, "out", "grade_lspilot.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

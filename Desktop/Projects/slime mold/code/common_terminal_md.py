#!/usr/bin/env python3
"""
common_terminal_md.py -- rescore Table 1's MD column on a common terminal set.

Why this exists
---------------
`network_metrics()` averages the mean inter-terminal distance MD over `tin`, the
terminals that fall inside the network's largest connected component. That set is
not the same for every row of Table 1:

    MST, Steiner, Physarum   reach all 38 terminals
    hub-and-spoke            reaches 28, because baseline_hub_spoke() routes only
                             to the 22 DEMAND_NODES

and the 10 it misses are the line ends (O'Hare, Midway, 95th/Dan Ryan, Howard,
Linden, Dempster-Skokie, Forest Park, Harlem/Lake, Cermak-54th, Clark/Lake). Mean
inter-terminal distance is exactly the quantity that shrinks when you delete the
most peripheral stops, so the tabulated MD flatters hub-and-spoke and the column
is not like-for-like.

This script recomputes MD for every network over the intersection of the terminal
sets, which is the comparison the text quotes. The headline Physarum network is
rebuilt from the saved conductivity field `data/_headline_D.npy` rather than
re-solved, so this is a pure post-processing step with no RNG in it -- it asserts
that the rebuild reproduces `data/_metrics.json` before reporting anything.

Run:  python common_terminal_md.py        (from code/)
"""
from __future__ import annotations

import json
import os

import networkx as nx
import numpy as np

import physarum as P
from run_experiments import GAMMA_HEAD, SPACING, repair

DAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
RTOL = 1e-9


def md_over(G: nx.Graph, termset) -> tuple[float, int]:
    """Mean shortest-path length (km) between terminals in `termset` that G reaches."""
    comp = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(comp)
    tin = [t for t in termset if t in Gc]
    lengths = dict(nx.all_pairs_dijkstra_path_length(Gc, weight="length"))
    total, count = 0.0, 0
    for i, a in enumerate(tin):
        for b in tin[i + 1:]:
            if b in lengths.get(a, {}):
                total += lengths[a][b]
                count += 1
    if count == 0:
        return float("nan"), len(tin)
    return total / count, len(tin)


def main() -> None:
    cg = P.CandidateGraph(spacing_km=SPACING, weight_coeffs=(0.4, 0.3, 0.3))
    D = np.load(os.path.join(DAT, "_headline_D.npy"))
    head = repair(cg, P.physarum_network(cg, D, frac=0.10)[0])

    nets = {
        "MST": P.baseline_mst(cg),
        "Steiner": P.baseline_steiner(cg),
        "Hub-and-spoke": P.baseline_hub_spoke(cg),
        f"Physarum g={GAMMA_HEAD}": head,
    }

    # Guard: the rebuilt headline must reproduce the stored run exactly, otherwise
    # the numbers below describe a different network than the paper's.
    stored = json.load(open(os.path.join(DAT, "_metrics.json")))
    got = P.network_metrics(head, cg)
    want = stored[f"Physarum g={GAMMA_HEAD}"]
    for key in ("TL", "MD", "FT"):
        np.testing.assert_allclose(
            got[key], want[key], rtol=RTOL,
            err_msg=f"rebuilt headline {key} does not match data/_metrics.json")
    print(f"rebuild reproduces data/_metrics.json to rtol={RTOL:g}\n")

    common: set = set(cg.terminals)
    for G in nets.values():
        comp = max(nx.connected_components(G), key=len)
        common &= {t for t in cg.terminals if t in comp}

    inv = {v: k for k, v in cg.terminal_index.items()}
    missed = sorted(set(cg.terminals) - common)
    print(f"terminals: {len(cg.terminals)} total, {len(common)} reached by every network")
    print(f"missed by hub-and-spoke: {[inv.get(t, t) for t in missed]}\n")

    print(f"{'network':24s} {'TL':>7s} {'MD_tab':>8s} {'n_tab':>6s} {'MD_common':>10s} {'n':>4s}")
    md_common = {}
    for name, G in nets.items():
        m = P.network_metrics(G, cg)
        md, n = md_over(G, common)
        md_common[name] = md
        print(f"{name:24s} {m['TL']:7.2f} {m['MD']:8.2f} {m['n_term']:6d} {md:10.2f} {n:4d}")

    hub = md_common["Hub-and-spoke"]
    print(f"\nMST - hub MD gap:  tabulated {P.network_metrics(nets['MST'], cg)['MD'] - hub:5.2f} km"
          f"   like-for-like {md_common['MST'] - hub:5.2f} km")
    phys = md_common[f"Physarum g={GAMMA_HEAD}"]
    print(f"Physarum advantage over hub-and-spoke: "
          f"tabulated {hub - P.network_metrics(nets[f'Physarum g={GAMMA_HEAD}'], cg)['MD']:5.2f} km"
          f"   like-for-like {hub - phys:5.2f} km")


if __name__ == "__main__":
    main()

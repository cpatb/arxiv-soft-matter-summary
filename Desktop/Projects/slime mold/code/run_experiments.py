"""
run_experiments.py
==================
Drives the full Chicago 'L' slime-mold study and writes all figures + results.

Stages (each fits the sandbox time budget):
  python run_experiments.py 1   -> headline run + baselines + figs 1-4 + csv
  python run_experiments.py 2   -> multi-seed stats + gamma sweep + figs 5-6 + json

Outputs (../figures, ../data):
  fig1_candidate_graph.png   constraints: lattice, lake, existing 'L', demand nodes
  fig2_evolution.png         Physarum conductivity field at 4 timesteps
  fig3_grown_network.png     headline grown network (gamma=1.3) over Chicago
  fig4_baselines.png         MST / Steiner / hub-spoke / Physarum side by side
  fig5_pareto.png            Tero cost-benefit: TL ratio vs FT and vs MD ratio
  fig6_gamma_sweep.png       TL, MD, FT versus adaptation exponent gamma
  results.json / results.csv numerical metrics (incl. multi-seed statistics)
"""
import json, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import networkx as nx

import data_chicago as dc
import physarum as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DAT = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(FIG, exist_ok=True)
os.makedirs(DAT, exist_ok=True)

LAKE = "#bcd6e8"; LAND = "#f4f2ec"
C_EXIST = "#444444"; C_NEW = "#c1121f"
C_POP = "#1d6f42"; C_ACT = "#1f5fa8"; C_UND = "#e08a00"
plt.rcParams.update({"font.size": 9, "figure.dpi": 140})

GAMMA_HEAD = 1.3
SPACING = 1.6


def lonlat(cg, n):
    x, y = cg.xy[n]
    return (dc.LON0 + x / dc._KM_PER_DEG_LON, dc.LAT0 + y / dc._KM_PER_DEG_LAT)


def draw_basemap(ax, cg):
    poly = np.array(dc.LAND_POLYGON)
    ax.set_facecolor(LAKE)
    ax.fill(poly[:, 0], poly[:, 1], color=LAND, zorder=0, ec="#9bb7c9", lw=0.8)
    for line in dc.EXISTING_LINES:
        xs = [dc.EXISTING_NODES[k][0] for k in line]
        ys = [dc.EXISTING_NODES[k][1] for k in line]
        ax.plot(xs, ys, color=C_EXIST, lw=1.6, alpha=0.55, zorder=2,
                solid_capstyle="round")
    ax.set_xlim(-87.97, -87.5); ax.set_ylim(41.62, 42.10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect(1.0 / np.cos(np.radians(41.85)))


def draw_network(ax, cg, G, color=C_NEW, widthkey="D", base=0.6, alpha=0.9):
    segs, ws = [], []
    for u, v, d in G.edges(data=True):
        segs.append([lonlat(cg, u), lonlat(cg, v)])
        ws.append(d[widthkey] if widthkey in d else 1.0)
    ws = np.array(ws)
    lw = base + 2.6 * (ws / ws.max()) if len(ws) else base
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw,
                                     alpha=alpha, zorder=3))


def draw_demand(ax, cg, labels=False):
    for k, attr in dc.DEMAND_NODES.items():
        lo, la = attr["lonlat"]
        roles = {"pop": attr["pop"], "act": attr["activity"], "und": attr["underserved"]}
        dom = max(roles, key=roles.get)
        col = {"pop": C_POP, "act": C_ACT, "und": C_UND}[dom]
        sz = 18 + 60 * dc.demand_weight(attr)
        ax.scatter([lo], [la], s=sz, c=col, ec="white", lw=0.6, zorder=5)
        if labels:
            ax.annotate(k.replace("_", " "), (lo, la), fontsize=5.0,
                        xytext=(2, 2), textcoords="offset points", zorder=6)
    for k, (lo, la) in dc.EXISTING_NODES.items():
        ax.scatter([lo], [la], s=14, marker="s", c=C_EXIST, ec="white", lw=0.5, zorder=5)


def fig1(cg):
    fig, ax = plt.subplots(figsize=(6.4, 7.0))
    draw_basemap(ax, cg)
    segs = [[lonlat(cg, int(i)), lonlat(cg, int(j))] for i, j in zip(cg.ei, cg.ej)]
    ax.add_collection(LineCollection(segs, colors="#cfcabb", lw=0.3, alpha=0.6, zorder=1))
    draw_demand(ax, cg, labels=True)
    leg = [Line2D([0], [0], color=C_EXIST, lw=2, label="existing 'L' corridor"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_POP, markersize=8, label="population centroid"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_ACT, markersize=8, label="activity centre"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_UND, markersize=8, label="transit desert"),
           Line2D([0], [0], color="#cfcabb", lw=1, label="candidate edge")]
    ax.legend(handles=leg, loc="lower left", fontsize=6.5, framealpha=0.9)
    ax.set_title(f"Constraint geometry: {cg.N} candidate nodes, {cg.E} edges\n"
                 "Lake Michigan = hard no-build (blue)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1_candidate_graph.png")); plt.close(fig)


def fig2(cg, res):
    snaps = res["snapshots"]
    pick = [snaps[0], snaps[len(snaps)//3], snaps[2*len(snaps)//3], snaps[-1]]
    fig, axs = plt.subplots(1, 4, figsize=(13.5, 4.2))
    for ax, (step, D) in zip(axs, pick):
        draw_basemap(ax, cg)
        Dn = D / D.max()
        segs = [[lonlat(cg, int(i)), lonlat(cg, int(j))] for i, j in zip(cg.ei, cg.ej)]
        ax.add_collection(LineCollection(segs, colors=C_NEW, linewidths=0.2 + 3.0 * Dn,
                                         alpha=0.75, zorder=3))
        for k, attr in dc.DEMAND_NODES.items():
            lo, la = attr["lonlat"]; ax.scatter([lo], [la], s=8, c="k", zorder=5)
        ax.set_title(f"step {step}", fontsize=9)
    fig.suptitle("Physarum conductivity field adapting toward the demand-weighted "
                 f"network (gamma={GAMMA_HEAD})", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_evolution.png")); plt.close(fig)


def fig3(cg, G):
    fig, ax = plt.subplots(figsize=(6.6, 7.2))
    draw_basemap(ax, cg)
    draw_network(ax, cg, G, color=C_NEW)
    draw_demand(ax, cg, labels=True)
    m = P.network_metrics(G, cg)
    ax.set_title("Slime-mold grown addition to the Chicago 'L' "
                 f"(gamma={GAMMA_HEAD})\n"
                 f"new track ~ {m['TL']:.0f} km, fault tolerance {m['FT']:.2f}, "
                 f"mean trip {m['MD']:.1f} km", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_grown_network.png")); plt.close(fig)


def fig4(cg, nets):
    fig, axs = plt.subplots(1, 4, figsize=(14, 4.6))
    order = ["MST", "Steiner", "Hub-and-spoke", f"Physarum g={GAMMA_HEAD}"]
    for ax, name in zip(axs, order):
        G = nets[name]
        draw_basemap(ax, cg)
        col = C_NEW if "Physarum" in name else "#7a4eab"
        draw_network(ax, cg, G, color=col,
                     widthkey="D" if "Physarum" in name else "none", base=1.0)
        for k, attr in dc.DEMAND_NODES.items():
            lo, la = attr["lonlat"]; ax.scatter([lo], [la], s=8, c="k", zorder=5)
        m = P.network_metrics(G, cg)
        ax.set_title(f"{name}\nTL={m['TL']:.0f}km  MD={m['MD']:.1f}  FT={m['FT']:.2f}", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_baselines.png")); plt.close(fig)


def fig5(rows, mst_ref):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    tl0, md0 = mst_ref["TL"], mst_ref["MD"]
    for r in rows:
        x = r["TL"] / tl0
        axs[0].scatter(x, r["FT"], s=70, c=r["color"], zorder=3)
        axs[0].annotate(r["label"], (x, r["FT"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
        axs[1].scatter(x, r["MD"] / md0, s=70, c=r["color"], zorder=3)
        axs[1].annotate(r["label"], (x, r["MD"]/md0), fontsize=7, xytext=(4, 3), textcoords="offset points")
    axs[0].set_xlabel("total length / MST length"); axs[0].set_ylabel("fault tolerance FT")
    axs[0].set_title("Cost vs. robustness"); axs[0].grid(alpha=0.3)
    axs[1].set_xlabel("total length / MST length"); axs[1].set_ylabel("mean trip distance / MST")
    axs[1].set_title("Cost vs. transport efficiency"); axs[1].grid(alpha=0.3); axs[1].invert_yaxis()
    fig.suptitle("Tero cost-benefit framework (MST = reference at x=1)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_pareto.png")); plt.close(fig)


def fig6(sweep):
    g = [s["gamma"] for s in sweep]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(g, [s["TL"] for s in sweep], "o-", color=C_NEW, label="total length TL (km)")
    ax.set_xlabel("adaptation exponent gamma"); ax.set_ylabel("TL (km)", color=C_NEW)
    ax.tick_params(axis="y", labelcolor=C_NEW)
    ax2 = ax.twinx()
    ax2.plot(g, [s["FT"] for s in sweep], "s--", color=C_ACT, label="fault tolerance")
    ax2.plot(g, [s["MD"]/sweep[0]["MD"] for s in sweep], "^:", color=C_POP, label="MD / MD(gmin)")
    ax2.set_ylabel("FT  and  normalised MD")
    ax.axvline(GAMMA_HEAD, color="gray", ls=":", alpha=0.6)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=7, loc="center right")
    ax.set_title("Regime control: mesh (g<1) -> balanced -> tree (g>1)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig6_gamma_sweep.png")); plt.close(fig)


def repair(cg, G):
    """Attach any demand terminal not in the main component via the cheapest path."""
    if G.number_of_edges() == 0:
        return G
    main = max(nx.connected_components(G), key=len)
    Gc = cg.edge_graph(weight="cost")
    for d in cg.demand_terms:
        if d in main:
            continue
        dist, path = nx.single_source_dijkstra(Gc, d, weight="cost")
        targets = [(dist[n], n) for n in main if n in dist]
        if not targets:
            continue
        _, tgt = min(targets)
        nodes = path[tgt]
        for a, b in zip(nodes[:-1], nodes[1:]):
            ek = P._edge_lookup(cg, a, b)
            G.add_edge(a, b, length=float(cg.length[ek]), cost=float(cg.cost[ek]), D=1e-3)
        main = max(nx.connected_components(G), key=len)
    return G


def stage1():
    print("Building candidate graph ...")
    cg = P.CandidateGraph(spacing_km=SPACING, weight_coeffs=(0.4, 0.3, 0.3))
    print(f"  N={cg.N} nodes, E={cg.E} edges, {len(cg.terminals)} terminals")
    print("Running headline Physarum (gamma=%.1f) ..." % GAMMA_HEAD)
    res = P.physarum_solve(cg, gamma=GAMMA_HEAD, steps=500, batch=12,
                           I0_scale=6.0, seed=3, record_every=25)
    Ghead = repair(cg, P.physarum_network(cg, res["D"], frac=0.10)[0])
    print("Baselines ...")
    nets = {"MST": P.baseline_mst(cg), "Steiner": P.baseline_steiner(cg),
            "Hub-and-spoke": P.baseline_hub_spoke(cg), f"Physarum g={GAMMA_HEAD}": Ghead}
    metrics = {name: P.network_metrics(G, cg) for name, G in nets.items()}
    print("Figures 1-4 ...")
    fig1(cg); fig2(cg, res); fig3(cg, Ghead); fig4(cg, nets)
    np.save(os.path.join(DAT, "_headline_D.npy"), res["D"])
    with open(os.path.join(DAT, "_metrics.json"), "w") as f:
        json.dump(metrics, f, indent=2)
    with open(os.path.join(DAT, "results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["network", "TL_km", "TC_cost", "MD_km", "FT", "n_edges", "connected"])
        for k, v in metrics.items():
            w.writerow([k, f"{v['TL']:.1f}", f"{v['TC']:.1f}", f"{v['MD']:.2f}",
                        f"{v['FT']:.3f}", v["n_edges"], v["connected"]])
    print("stage1 done.")
    for k, v in metrics.items():
        print(f"{k:18s} TL={v['TL']:6.1f}km MD={v['MD']:5.1f}km FT={v['FT']:.2f} "
              f"edges={v['n_edges']:3d} conn={v['connected']}")


def stage2():
    cg = P.CandidateGraph(spacing_km=SPACING, weight_coeffs=(0.4, 0.3, 0.3))
    with open(os.path.join(DAT, "_metrics.json")) as f:
        metrics = json.load(f)
    print("Multi-seed statistics ...")
    seed_m = []
    for sd in range(10, 18):
        r = P.physarum_solve(cg, gamma=GAMMA_HEAD, steps=400, batch=12, I0_scale=6.0, seed=sd)
        G = repair(cg, P.physarum_network(cg, r["D"], frac=0.10)[0])
        seed_m.append(P.network_metrics(G, cg))
    agg = {k: [float(np.mean([m[k] for m in seed_m])), float(np.std([m[k] for m in seed_m]))]
           for k in ["TL", "TC", "MD", "FT"]}
    print("Gamma sweep ...")
    sweep = []
    for g in [0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.2]:
        r = P.physarum_solve(cg, gamma=g, steps=400, batch=12, I0_scale=6.0, seed=7)
        G = repair(cg, P.physarum_network(cg, r["D"], frac=0.10)[0])
        m = P.network_metrics(G, cg); m["gamma"] = g
        sweep.append(m)
    print("Figures 5-6 ...")
    pareto_rows = [dict(label="MST", color="#7a4eab", **metrics["MST"]),
                   dict(label="Steiner", color="#2a9d8f", **metrics["Steiner"]),
                   dict(label="Hub-spoke", color="#888888", **metrics["Hub-and-spoke"])]
    pareto_rows += [dict(label=f"g={s['gamma']}", color=C_NEW, **s) for s in sweep]
    fig5(pareto_rows, metrics["MST"]); fig6(sweep)
    out = {"candidate_graph": dict(N=cg.N, E=cg.E, terminals=len(cg.terminals), spacing_km=SPACING),
           "headline_gamma": GAMMA_HEAD, "networks": metrics,
           "headline_multiseed": agg, "gamma_sweep": sweep}
    with open(os.path.join(DAT, "results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("headline multi-seed (mean,std):", agg)
    print("stage2 done.")


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "1"
    stage1() if stage == "1" else stage2()

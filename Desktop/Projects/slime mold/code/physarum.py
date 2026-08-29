import numpy as np
import networkx as nx
from scipy.spatial import Delaunay
from matplotlib.path import Path
import data_chicago as dc


def _point_to_segment_km(px, py, ax, ay, bx, by):
    abx, aby = bx - ax, by - ay
    denom = abx * abx + aby * aby
    if denom == 0:
        return np.hypot(px - ax, py - ay)
    t = np.clip(((px - ax) * abx + (py - ay) * aby) / denom, 0.0, 1.0)
    cx, cy = ax + t * abx, ay + t * aby
    return np.hypot(px - cx, py - cy)


class CandidateGraph:
    def __init__(self, spacing_km=1.6, downtown_premium=1.6, row_discount=0.6,
                 weight_coeffs=(0.4, 0.3, 0.3), seed=0):
        self.spacing = spacing_km
        self.downtown_premium = downtown_premium
        self.row_discount = row_discount
        self.a, self.b, self.c = weight_coeffs
        self.rng = np.random.default_rng(seed)
        self._build()

    def _land_path(self):
        poly_xy = [dc.project(lo, la) for lo, la in dc.LAND_POLYGON]
        return Path([(float(x), float(y)) for x, y in poly_xy])

    def _existing_segments_km(self):
        segs = []
        for line in dc.EXISTING_LINES:
            for u, v in zip(line[:-1], line[1:]):
                ux, uy = dc.project(*dc.EXISTING_NODES[u])
                vx, vy = dc.project(*dc.EXISTING_NODES[v])
                segs.append((float(ux), float(uy), float(vx), float(vy)))
        return segs

    def _cost_multiplier(self, mx, my):
        dmin = min(_point_to_segment_km(mx, my, *s) for s in self._segs)
        if dmin < 0.7:
            return self.row_discount
        lx, ly = dc.project(*dc.EXISTING_NODES["Loop"])
        if np.hypot(mx - lx, my - ly) < 2.5:
            return self.downtown_premium
        return 1.0

    def _build(self):
        path = self._land_path()
        self._segs = self._existing_segments_km()
        vx = np.array([dc.project(lo, la)[0] for lo, la in dc.LAND_POLYGON])
        vy = np.array([dc.project(lo, la)[1] for lo, la in dc.LAND_POLYGON])
        x0, x1, y0, y1 = vx.min(), vx.max(), vy.min(), vy.max()
        s = self.spacing
        pts = []
        row = 0
        y = y0
        while y <= y1:
            xoff = 0.0 if row % 2 == 0 else s / 2.0
            x = x0 + xoff
            while x <= x1:
                pts.append((x, y))
                x += s
            y += s * np.sqrt(3) / 2.0
            row += 1
        pts = np.array(pts)
        inside = path.contains_points(pts, radius=-1e-9)
        lattice = pts[inside]
        node_xy = list(map(tuple, lattice))
        self.terminal_index = {}
        self.node_key = {}

        def add_exact(key, lon, lat):
            x, y = dc.project(lon, lat)
            nid = len(node_xy)
            node_xy.append((float(x), float(y)))
            self.terminal_index[key] = nid
            self.node_key[nid] = key

        for k, (lo, la) in dc.EXISTING_NODES.items():
            add_exact(k, lo, la)
        for k, attr in dc.DEMAND_NODES.items():
            add_exact(k, *attr["lonlat"])
        self.xy = np.array(node_xy)
        self.N = len(self.xy)
        tri = Delaunay(self.xy)
        emax = 2.3 * s
        edges = set()
        for simplex in tri.simplices:
            for a, b in ((simplex[0], simplex[1]), (simplex[1], simplex[2]), (simplex[0], simplex[2])):
                i, j = (a, b) if a < b else (b, a)
                edges.add((i, j))
        ei, ej, length, cost = [], [], [], []
        for i, j in edges:
            xi, yi = self.xy[i]
            xj, yj = self.xy[j]
            L = np.hypot(xi - xj, yi - yj)
            if L > emax:
                continue
            mx, my = 0.5 * (xi + xj), 0.5 * (yi + yj)
            if not path.contains_point((mx, my)):
                continue
            ei.append(i); ej.append(j); length.append(L)
            cost.append(L * self._cost_multiplier(mx, my))
        self.ei = np.array(ei)
        self.ej = np.array(ej)
        self.length = np.array(length)
        self.cost = np.array(cost)
        self.E = len(self.ei)
        self.terminals = list(self.terminal_index.values())
        w = {}
        for k in dc.EXISTING_NODES:
            w[self.terminal_index[k]] = 1.0 if k == "Loop" else 0.5
        for k, attr in dc.DEMAND_NODES.items():
            w[self.terminal_index[k]] = dc.demand_weight(attr, self.a, self.b, self.c)
        self.term_weight = np.array([w[t] for t in self.terminals])
        self.demand_terms = [self.terminal_index[k] for k in dc.DEMAND_NODES]
        self.existing_terms = [self.terminal_index[k] for k in dc.EXISTING_NODES]
        G = nx.Graph(); G.add_nodes_from(range(self.N))
        G.add_edges_from(zip(self.ei.tolist(), self.ej.tolist()))
        assert nx.is_connected(G), "candidate graph is disconnected"
        self._G_full = G

    def edge_graph(self, weight="length"):
        G = nx.Graph(); G.add_nodes_from(range(self.N))
        for i, j, l, c in zip(self.ei, self.ej, self.length, self.cost):
            G.add_edge(int(i), int(j), length=float(l), cost=float(c))
        return G


def physarum_solve(cg, gamma=1.0, mu=1.0, dt=0.10, steps=350,
                   batch=30, I0_scale=1.0, sat=True, seed=1, record_every=0):
    rng = np.random.default_rng(seed)
    N, E = cg.N, cg.E
    ei, ej = cg.ei, cg.ej
    Leff = cg.cost.copy()
    D = 1.0 + 0.05 * rng.standard_normal(E)
    D = np.clip(D, 0.2, None)
    terms = np.array(cg.terminals)
    wterm = cg.term_weight / cg.term_weight.sum()
    ground = cg.terminal_index["Loop"]
    free = np.array([n for n in range(N) if n != ground])
    snapshots = []
    flux_hist = []
    for step in range(steps):
        cond = D / Leff
        Lap = np.zeros((N, N))
        np.add.at(Lap, (ei, ei), cond)
        np.add.at(Lap, (ej, ej), cond)
        np.add.at(Lap, (ei, ej), -cond)
        np.add.at(Lap, (ej, ei), -cond)
        A = Lap[np.ix_(free, free)]
        srcs = rng.choice(len(terms), size=batch, p=wterm)
        snks = rng.choice(len(terms), size=batch, p=wterm)
        B = np.zeros((N, batch))
        valid = []
        for k in range(batch):
            s, t = terms[srcs[k]], terms[snks[k]]
            if s == t:
                continue
            I0 = I0_scale * np.sqrt(cg.term_weight[srcs[k]] * cg.term_weight[snks[k]])
            B[s, k] += I0; B[t, k] -= I0; valid.append(k)
        if not valid:
            continue
        Bf = B[free][:, valid]
        Pf = np.linalg.solve(A, Bf)
        Pfull = np.zeros((N, len(valid)))
        Pfull[free] = Pf
        dP = Pfull[ei] - Pfull[ej]
        Q = cond[:, None] * dP
        Qmag = np.abs(Q).mean(axis=1)
        if sat:
            fQ = Qmag ** gamma / (1.0 + Qmag ** gamma)
        else:
            fQ = Qmag ** gamma
        D = D + dt * (fQ - mu * D)
        D = np.clip(D, 1e-6, None)
        flux_hist.append(float(Qmag.mean()))
        if record_every and (step % record_every == 0 or step == steps - 1):
            snapshots.append((step, D.copy()))
    return dict(D=D, Leff=Leff, flux_hist=flux_hist, snapshots=snapshots, ground=ground)


def physarum_network(cg, D, frac=0.08, relative=True):
    thr = frac * D.max() if relative else frac
    keep = D >= thr
    G = nx.Graph(); G.add_nodes_from(range(cg.N))
    for k in np.where(keep)[0]:
        i, j = int(cg.ei[k]), int(cg.ej[k])
        G.add_edge(i, j, length=float(cg.length[k]), cost=float(cg.cost[k]), D=float(D[k]))
    termset = set(cg.terminals)
    changed = True
    while changed:
        changed = False
        for n in [n for n in G.nodes if G.degree(n) == 1 and n not in termset]:
            G.remove_node(n); changed = True
    G.remove_nodes_from([n for n in list(G.nodes) if G.degree(n) == 0 and n not in termset])
    return G, keep


_EDGE_MAP = None


def _edge_lookup(cg, u, v):
    global _EDGE_MAP
    if _EDGE_MAP is None or _EDGE_MAP[0] is not cg:
        d = {}
        for k in range(cg.E):
            i, j = int(cg.ei[k]), int(cg.ej[k])
            d[(i, j)] = k; d[(j, i)] = k
        _EDGE_MAP = (cg, d)
    return _EDGE_MAP[1][(u, v)]


def _metric_closure_tree(cg, weight="cost"):
    G = cg.edge_graph(weight="length")
    for u, v, d in G.edges(data=True):
        d["w"] = d["length"] if weight == "length" else float(cg.cost[_edge_lookup(cg, u, v)])
    terms = cg.terminals
    H = nx.Graph(); paths = {}
    for s in terms:
        dist, path = nx.single_source_dijkstra(G, s, weight="w")
        for t in terms:
            if t > s and t in dist:
                H.add_edge(s, t, w=dist[t]); paths[(s, t)] = path[t]
    T = nx.minimum_spanning_tree(H, weight="w")
    out = nx.Graph(); out.add_nodes_from(range(cg.N))
    for s, t in T.edges():
        p = paths[(s, t)] if (s, t) in paths else paths[(t, s)]
        for a, b in zip(p[:-1], p[1:]):
            ek = _edge_lookup(cg, a, b)
            out.add_edge(a, b, length=float(cg.length[ek]), cost=float(cg.cost[ek]))
    return out


def baseline_mst(cg):
    return _metric_closure_tree(cg, weight="cost")


def baseline_steiner(cg):
    G = cg.edge_graph(weight="cost")
    T = nx.algorithms.approximation.steiner_tree(G, cg.terminals, weight="cost", method="mehlhorn")
    out = nx.Graph(); out.add_nodes_from(range(cg.N))
    for u, v in T.edges():
        ek = _edge_lookup(cg, u, v)
        out.add_edge(u, v, length=float(cg.length[ek]), cost=float(cg.cost[ek]))
    return out


def baseline_hub_spoke(cg, hub_key="Loop"):
    G = cg.edge_graph(weight="cost")
    hub = cg.terminal_index[hub_key]
    out = nx.Graph(); out.add_nodes_from(range(cg.N))
    dist, path = nx.single_source_dijkstra(G, hub, weight="cost")
    for d in cg.demand_terms:
        if d not in path:
            continue
        for a, b in zip(path[d][:-1], path[d][1:]):
            ek = _edge_lookup(cg, a, b)
            out.add_edge(a, b, length=float(cg.length[ek]), cost=float(cg.cost[ek]))
    return out


def network_metrics(G, cg):
    edges = list(G.edges(data=True))
    TL = sum(d["length"] for _, _, d in edges)
    TC = sum(d["cost"] for _, _, d in edges)
    sub_nodes = max(nx.connected_components(G), key=len) if G.number_of_edges() else set()
    connected = all(t in sub_nodes for t in cg.demand_terms)
    Gc = G.subgraph(sub_nodes).copy()
    tin = [t for t in cg.terminals if t in Gc]
    dsum, dcnt = 0.0, 0
    lengths = dict(nx.all_pairs_dijkstra_path_length(Gc, weight="length"))
    for i, s in enumerate(tin):
        for t in tin[i + 1:]:
            if t in lengths.get(s, {}):
                dsum += lengths[s][t]; dcnt += 1
    MD = dsum / dcnt if dcnt else np.inf
    E = Gc.number_of_edges()
    if E == 0:
        FT = 0.0
    else:
        termset = set(tin); crit = 0
        for u, v in list(nx.bridges(Gc)):
            Gc.remove_edge(u, v)
            side = nx.node_connected_component(Gc, u)
            Gc.add_edge(u, v, length=0.0)
            if any(t in side for t in termset) and any(t not in side for t in termset):
                crit += 1
        FT = 1.0 - crit / E
    return dict(TL=TL, TC=TC, MD=MD, FT=FT, connected=connected,
                n_edges=G.number_of_edges(), n_term=len(tin))

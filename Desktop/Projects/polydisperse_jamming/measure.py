"""measure.py — structural observables of a contact network.

Split from `packing.py` so that everything here is a PURE FUNCTION of a contact graph
and can be tested on Windows without Yade. `packing.contact_graph` builds the inputs.

Observables, and why each is here
---------------------------------
`rattler_fraction`     grains carrying no load. The M3 candidate: broader PSDs should
                       bury small grains in the interstices where they are rattlers.
`backbone`             the load-bearing subset after iterated rattler removal.
`z_by_size`            coordination resolved by grain size. If polydispersity works
                       through "small grains stop carrying load", this is where it
                       shows, and a mean Z cannot see it.
`force_moments`        the contact-force distribution. Its exponential tail above the
                       mean is a published signature (anchor A6).
`fabric_tensor`        contact-orientation anisotropy. Distinguishes an isotropically
                       compressed pack from a sheared one.
`spanning_cluster`     does the load-bearing network percolate the periodic cell? This
                       is the discriminating negative for "jammed": a pack below phi_J
                       must NOT span.
"""

from __future__ import annotations

import numpy as np

# --- the rattler threshold is a DEFINITION, and it changes the answer -------------
#
# Iteratively removing grains with fewer than k contacts is k-core decomposition of the
# contact graph. For a frictional pack at mean Z ~ 4.86 the k = 4 core is EMPTY -- the
# whole pack is classified as rattlers -- while a single non-iterated pass removes only
# 15% of it. That is not a bug in the packing: the k-core percolation threshold for a
# random graph with k = 4 sits near mean degree 5.1, so a marginal frictional pack is on
# the wrong side of it. Contact graphs are spatially correlated rather than random, so
# the collapse is not guaranteed, but it is observed here and it must not be hidden.
#
# So the threshold is NOT fixed at one value and quoted. Three are reported:
#
#   k = 2  TRUE rattlers: a grain with 0 or 1 contacts cannot transmit any force.
#          This is the conservative definition and it is the PRIMARY one, because it
#          makes no stability assumption that depends on friction.
#   k = 3  the frictional stability threshold (a frictional sphere in 3D can be
#          constrained by 3 contacts).
#   k = 4  the frictionless stability threshold (d+1). Correct for the mu = 0 anchor,
#          and the one that collapses for marginal frictional packs.
#
# Every result that depends on the choice must be graded at all three, the same way the
# lab grades a regularisation parameter corner by corner.
RATTLER_THRESHOLDS = (2, 3, 4)
RATTLER_MIN_CONTACTS = 2          # PRIMARY definition; see above


def iterated_rattlers(counts: dict[int, int], neigh: dict[int, set[int]],
                      min_contacts: int = RATTLER_MIN_CONTACTS) -> set[int]:
    """Grains with fewer than `min_contacts` load-bearing contacts, found iteratively.

    Removing a rattler lowers its neighbours' counts, which can make THEM rattlers, so
    a single pass systematically under-counts. The under-count grows with
    polydispersity -- exactly the axis this project varies -- so a one-pass count would
    manufacture a CV trend out of nothing.
    """
    counts = dict(counts)
    alive = set(counts)
    changed = True
    while changed:
        changed = False
        for b in list(alive):
            if counts[b] < min_contacts:
                alive.discard(b)
                for nb in neigh.get(b, ()):
                    if nb in alive:
                        counts[nb] -= 1
                changed = True
    return set(counts) - alive


def backbone(counts: dict[int, int], neigh: dict[int, set[int]],
             min_contacts: int = RATTLER_MIN_CONTACTS) -> set[int]:
    """The load-bearing subset: everything that is not a rattler."""
    return set(counts) - iterated_rattlers(counts, neigh, min_contacts)


def rattler_fraction(counts: dict[int, int], neigh: dict[int, set[int]],
                     min_contacts: int = RATTLER_MIN_CONTACTS) -> float:
    if not counts:
        return float("nan")
    return len(iterated_rattlers(counts, neigh, min_contacts)) / len(counts)


def z_by_size(counts: dict[int, int], neigh: dict[int, set[int]],
              diameters: dict[int, float], n_bins: int = 5,
              exclude_rattlers: bool = False) -> dict[str, list[float]]:
    """Mean coordination number resolved by grain diameter.

    Returns bin centres (in units of the mean diameter), mean Z per bin, and the count
    per bin. `n` per bin is ALWAYS returned so a bin with three grains cannot be read
    as a trend.
    """
    ids = sorted(counts)
    if not ids:
        return {"d_centre": [], "z_mean": [], "n": []}
    if exclude_rattlers:
        bb = backbone(counts, neigh)
        ids = [i for i in ids if i in bb]
        if not ids:
            return {"d_centre": [], "z_mean": [], "n": []}
    d = np.array([diameters[i] for i in ids], dtype=float)
    z = np.array([counts[i] for i in ids], dtype=float)
    dbar = float(d.mean())

    # Quantile edges COLLAPSE when diameters tie, and the monodisperse CV = 0 case ties
    # every one of them. Deduplicate, and fall back to linear spacing when the quantile
    # edges are degenerate but the data is not -- otherwise a discrete size set (4 big,
    # 16 small) silently returns a single bin and any size trend in it is invisible.
    edges = np.unique(np.quantile(d, np.linspace(0, 1, n_bins + 1)))
    if edges.size < 3 and np.ptp(d) > 0:
        edges = np.linspace(d.min(), d.max(), min(n_bins, np.unique(d).size) + 1)
    if edges.size < 2:                       # genuinely monodisperse: one bin
        edges = np.array([d.min(), d.max()])
    edges = edges.astype(float)
    edges[-1] = np.nextafter(edges[-1], np.inf) * (1.0 + 1e-9)

    out: dict[str, list[float]] = {"d_centre": [], "z_mean": [], "n": []}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (d >= lo) & (d < hi)
        if not m.any():
            continue
        out["d_centre"].append(float(d[m].mean() / dbar))
        out["z_mean"].append(float(z[m].mean()))
        out["n"].append(int(m.sum()))
    return out


def force_moments(forces: np.ndarray) -> dict[str, float]:
    """Moments and tail statistics of the contact normal-force distribution.

    `tail_slope` is the least-squares slope of ln P(f/<f>) over 1 <= f/<f> <= 4, the
    window in which the published exponential tail is quoted. Returned with the number
    of contacts in that window, because a slope from nine contacts is not a slope.
    """
    f = np.asarray(forces, dtype=float)
    f = f[np.isfinite(f) & (f > 0)]
    if f.size < 20:
        return {"n": int(f.size), "mean": float("nan"), "cv": float("nan"),
                "tail_slope": float("nan"), "n_tail": 0}
    fbar = float(f.mean())
    x = f / fbar
    sel = (x >= 1.0) & (x <= 4.0)
    slope = float("nan")
    if sel.sum() >= 20:
        hist, edges = np.histogram(x[sel], bins=12, density=True)
        c = 0.5 * (edges[:-1] + edges[1:])
        ok = hist > 0
        if ok.sum() >= 4:
            slope = float(np.polyfit(c[ok], np.log(hist[ok]), 1)[0])
    return {"n": int(f.size), "mean": fbar, "cv": float(f.std() / fbar),
            "tail_slope": slope, "n_tail": int(sel.sum())}


def fabric_tensor(normals: np.ndarray) -> dict[str, float]:
    """Second-order fabric tensor F_ij = <n_i n_j> and its anisotropy.

    For an isotropic contact network F -> I/3, so the deviator vanishes. `aniso` is
    the von Mises norm of the deviator, zero for isotropic.
    """
    n = np.asarray(normals, dtype=float)
    if n.ndim != 2 or n.shape[1] != 3 or n.shape[0] < 10:
        return {"aniso": float("nan"), "n": int(n.shape[0]) if n.ndim == 2 else 0}
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    F = (n[:, :, None] * n[:, None, :]).mean(axis=0)
    dev = F - np.eye(3) / 3.0
    return {"aniso": float(np.sqrt(1.5 * np.sum(dev * dev))), "n": int(n.shape[0])}


def rattler_report(counts: dict[int, int], neigh: dict[int, set[int]],
                   thresholds=RATTLER_THRESHOLDS,
                   volumes: dict[int, float] | None = None) -> dict[str, float]:
    """Rattler fraction and non-rattler Z at EVERY threshold, in one dict.

    Reporting all of them is the point: a conclusion that only holds at one k is a
    statement about the definition, not about the material.

    `volumes` (grain id -> volume) is OPTIONAL and added 2026-09-19. When given, the
    VOLUME fraction of the rattler set is reported beside its COUNT fraction.

    Why it matters: the exact decomposition of the packing fraction into a backbone and
    the grains that fell out of it is `phi = phi_B / (1 - psi)` with `psi` the rattler
    VOLUME fraction, and until now **every rattler field on disk was a COUNT**. Rattlers
    are preferentially small, so the two differ by a factor that depends on the size
    distribution -- which is precisely the dependence the project is trying to measure.

    ⚠ Passing nothing changes nothing: with `volumes = None` the output dict is
    byte-identical to what this function returned before, so every existing record and
    every existing caller is unaffected by construction. That is the bit-identity
    argument, and `test_measure.py` asserts it rather than trusting it.
    """
    out: dict[str, float] = {}
    v_tot = None
    if volumes:
        v_tot = sum(volumes.get(g, 0.0) for g in counts)
    for k in thresholds:
        rat = iterated_rattlers(counts, neigh, min_contacts=k)
        bb = set(counts) - rat
        out[f"rattler_frac_k{k}"] = len(rat) / len(counts) if counts else float("nan")
        if v_tot:
            out[f"rattler_vol_frac_k{k}"] = sum(volumes.get(g, 0.0) for g in rat) / v_tot
        if bb:
            inner = sum(1 for b in bb for nb in neigh.get(b, ()) if nb in bb)
            out[f"z_nr_k{k}"] = inner / len(bb)      # each contact counted twice
            out[f"largest_frac_k{k}"] = spanning_cluster(neigh, bb)["largest_frac"]
        else:
            out[f"z_nr_k{k}"] = float("nan")
            out[f"largest_frac_k{k}"] = 0.0
    return out


def rattler_by_class(counts: dict[int, int], neigh: dict[int, set[int]],
                     diameters: dict[int, float], d_split: float,
                     thresholds=RATTLER_THRESHOLDS) -> dict[str, float]:
    """Rattler fraction and coordination resolved into TWO size classes for a binary.

    Grains with diameter < `d_split` are "small", the rest "large". For a binary pack
    `d_split` is the geometric mean of the two design diameters, so the split is exact
    and does not depend on binning. Everything returned is per class:

      n_S, n_L                       grain counts (n_S + n_L = N)
      vol_S                          volume fraction of the solid carried by the small class
      rat_S_k{k}, rat_L_k{k}         rattler fraction BY NUMBER within each class
      ratvol_k{k}                    rattler fraction BY VOLUME over the whole pack
      z_S, z_L                       mean contact count per class (all grains)
      z_SS, z_SL, z_LS, z_LL         mean number of small/large neighbours per small/large grain
      nc_SS, nc_SL, nc_LL            contact counts by pair type (each contact once)

    This exists because the pairwise-filling mechanism makes a prediction the pooled
    rattler fraction cannot test: on the coarse-dominant branch the small grains sit in
    the interstices of a load-bearing large skeleton, so `rat_S` is high and `rat_L`
    low; on the fine-dominant branch the large grains are inclusions in a small-grain
    matrix and BOTH classes carry load. The pooled fraction is number-weighted and so
    reads the small class almost exclusively, which is why it peaks away from the
    density maximum.
    """
    ids = sorted(counts)
    if not ids or d_split <= 0.0:
        return {"n_S": 0, "n_L": 0}
    small = {i for i in ids if diameters[i] < d_split}
    large = set(ids) - small
    vol = {i: diameters[i] ** 3 for i in ids}          # pi/6 cancels in every ratio
    vtot = sum(vol.values())
    out: dict[str, float] = {"n_S": len(small), "n_L": len(large),
                             "vol_S": sum(vol[i] for i in small) / vtot}

    for k in thresholds:
        rat = iterated_rattlers(counts, neigh, min_contacts=k)
        rs, rl = rat & small, rat & large
        out[f"rat_S_k{k}"] = len(rs) / len(small) if small else float("nan")
        out[f"rat_L_k{k}"] = len(rl) / len(large) if large else float("nan")
        out[f"ratvol_k{k}"] = sum(vol[i] for i in rat) / vtot

    def _mean(cls: set[int], sel) -> float:
        if not cls:
            return float("nan")
        return sum(sum(1 for nb in neigh.get(i, ()) if sel(nb)) for i in cls) / len(cls)

    out["z_S"] = _mean(small, lambda nb: True)
    out["z_L"] = _mean(large, lambda nb: True)
    out["z_SS"] = _mean(small, lambda nb: nb in small)
    out["z_SL"] = _mean(small, lambda nb: nb in large)
    out["z_LS"] = _mean(large, lambda nb: nb in small)
    out["z_LL"] = _mean(large, lambda nb: nb in large)
    # each undirected contact appears in two adjacency sets; halve the same-class sums
    out["nc_SS"] = int(round(sum(1 for i in small for nb in neigh.get(i, ()) if nb in small) / 2))
    out["nc_LL"] = int(round(sum(1 for i in large for nb in neigh.get(i, ()) if nb in large) / 2))
    out["nc_SL"] = int(sum(1 for i in small for nb in neigh.get(i, ()) if nb in large))
    return out


def spanning_cluster(neigh: dict[int, set[int]], members: set[int]) -> dict[str, float]:
    """Largest connected component of the load-bearing network, as a fraction.

    In a PERIODIC cell every configuration is topologically "spanning", so percolation
    across a boundary is not the useful test. What is useful, and what this returns, is
    whether the load-bearing network is one connected object (fraction -> 1) or a dust
    of disconnected clumps (fraction << 1). A pack below phi_J gives the latter.
    """
    if not members:
        return {"largest_frac": 0.0, "n_components": 0}
    unseen = set(members)
    sizes = []
    while unseen:
        seed = unseen.pop()
        comp, stack = {seed}, [seed]
        while stack:
            b = stack.pop()
            for nb in neigh.get(b, ()):
                if nb in unseen:
                    unseen.discard(nb)
                    comp.add(nb)
                    stack.append(nb)
        sizes.append(len(comp))
    return {"largest_frac": max(sizes) / len(members), "n_components": len(sizes)}

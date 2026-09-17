"""make_paper_figures.py -- the three manuscript figures, from the same loaders the graders use.

    python paper/make_paper_figures.py           # build fig1..fig3 (.pdf, .png, .csv, .md)
    python paper/make_paper_figures.py --check   # build twice, require byte-identical PDFs

Every plotted point is written to a sibling .csv, and the .md beside each figure carries
the explanation that the house style keeps off the axes. Linear axes start at zero.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(HERE, "figures")
sys.path.insert(0, ROOT)

import sweep_descriptive as sd  # noqa: E402
import grade_mixture as gm  # noqa: E402
import grade_fabric as gf  # noqa: E402
from phase12_mixture import MIN_PHASE_COUNT, N_PROD, X_ANCHOR  # noqa: E402
from psd import BinaryPSD  # noqa: E402

matplotlib.rcParams["pdf.compression"] = 0
matplotlib.rcParams["svg.hashsalt"] = "polydisperse-jamming-paper"
matplotlib.rcParams["font.size"] = 8
matplotlib.rcParams["axes.linewidth"] = 0.7
matplotlib.rcParams["figure.dpi"] = 150
matplotlib.rcParams["mathtext.fontset"] = "cm"

FAM = ("lognormal", "uniform", "bidisperse")
STYLE = {"lognormal": dict(marker="o", color="#1f4e79"),
         "uniform": dict(marker="s", color="#c0392b"),
         "bidisperse": dict(marker="^", color="#2e8b57")}
COL_W = 3.375  # one PRL column, inches


def save(fig, name: str, rows: list[dict], note: str) -> None:
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight",
                    metadata={"CreationDate": None} if ext == "pdf" else None)
    plt.close(fig)
    with open(os.path.join(FIG, f"{name}.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(FIG, f"{name}.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(note.strip() + "\n")


def mean_se(v):
    v = np.asarray(v, float)
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), len(v)


# ------------------------------------------------------------------ figure 1
def fig1(fits):
    rows = []
    fig, (a, b) = plt.subplots(1, 2, figsize=(COL_W, 1.7), gridspec_kw={"width_ratios": [1.6, 1]})
    for fam in FAM:
        g = fits[(fits.family == fam) & (fits.mu_p == sd.MU_REF)]
        xs, ys, es = [], [], []
        for cv in sd.CVS:
            v = g[g.cv_nom == cv].phi_J_lin
            m, se, n = mean_se(v)
            xs.append(cv); ys.append(m); es.append(se)
            rows.append({"panel": "a", "family": fam, "cv": cv, "phi_J": m, "se": se, "n": n})
        a.errorbar(xs, ys, yerr=es, ms=3.5, lw=0.8, capsize=1.5, ls="-", **STYLE[fam])
        shift = sd.paired_shift(fits, fam)
        assert shift is not None, f"paired shift undefined for {fam}"
        m, se, n = shift
        b.bar(fam[:4], m, yerr=se, color=STYLE[fam]["color"], width=0.6, capsize=2, lw=0)
        rows.append({"panel": "b", "family": fam, "cv": "0->0.5", "phi_J": m, "se": se, "n": n})
    a.set_xlim(0, 0.55); a.set_ylim(0, 0.7)
    a.set_xlabel("CV"); a.set_ylabel(r"$\phi_J$")
    a.set_xticks([0, 0.25, 0.5])
    b.set_ylim(0, 0.08); b.set_ylabel(r"$\Delta\phi_J$")
    b.set_yticks([0, 0.04, 0.08])
    for ax, lab in ((a, "(a)"), (b, "(b)")):
        ax.text(0.03, 0.95, lab, transform=ax.transAxes, va="top")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig1_phiJ", rows, r"""
# fig1_phiJ — the jamming fraction against the width of the size distribution

(a) `φ_J` extrapolated to `p̂ → 0` from a four-rung pressure ladder, mean ± se over
three seeds, at `μ_p = 0.3`, N = 1500, for three size-distribution families at matched
number mean and matched CV: lognormal (circles), uniform (squares, zero skewness), and
number-symmetric bidisperse (triangles, zero skewness, all weight at the two extremes).
(b) The paired CV 0 → 0.5 change per family; the two ends share a seed and the
difference is taken per seed before averaging.

Y axes start at zero by house rule so the eye reads the size of the effect: the family
spread in (b) is about 5 % of `φ_J` itself. The bidisperse point at CV 0.5 is a size-ratio-3
mixture with 3.6 % of its volume in fines, i.e. that family is a path through the
(ratio, composition) plane, not a representative bidisperse packing.
""")


# ------------------------------------------------------------------ figure 2
def fig2(sel, fab_rows):
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(2 * COL_W, 1.7))
    getters = (("rattler", lambda r: r["rattler"]["rattler_frac_k2"], r"$f_{\rm r}$", (0, 0.6)),
               ("Z_nr", lambda r: r["rattler"]["z_nr_k2"], r"$Z_{\rm nr}$", (0, 6)),
               ("T_nr", None, r"$T_{\rm nr}$", (0, 1.2)))
    for k, (name, getter, ylab, ylim) in enumerate(getters):
        for fam in FAM:
            xs, ys, es = [], [], []
            for cv in sd.CVS:
                if getter is not None:
                    v = [getter(r) for r in sel
                         if (r.get("family") or "lognormal") == fam and r["cv"] == cv
                         and r["p_hat_target"] == 1e-3]
                else:
                    v = [t for c, t in gf.tstat([r for r in fab_rows
                                                  if (r.get("family") or "lognormal") == fam
                                                  and r["p_hat_target"] == gf.P_PRIMARY], "nr")
                         if c == cv]
                m, se, n = mean_se(v)
                xs.append(cv); ys.append(m); es.append(se)
                rows.append({"panel": "abc"[k], "observable": name, "family": fam, "cv": cv,
                             "mean": m, "se": se, "n": n})
            ax[k].errorbar(xs, ys, yerr=es, ms=3.5, lw=0.8, capsize=1.5, ls="-", **STYLE[fam])
        ax[k].set_xlim(0, 0.55); ax[k].set_ylim(*ylim)
        ax[k].set_xlabel("CV"); ax[k].set_ylabel(ylab)
        ax[k].set_xticks([0, 0.25, 0.5])
        ax[k].text(0.03, 0.95, f"({'abc'[k]})", transform=ax[k].transAxes, va="top")
    ax[2].axhline(1.0, color="0.5", lw=0.6, ls=":")
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig2_backbone", rows, r"""
# fig2_backbone — polydispersity removes grains from the network; it does not restructure it

All panels at `p̂ = 10⁻³`, `μ_p = 0.3`, N = 1500, mean ± se over seeds; same family symbols
as Fig. 1. (a) The rattler fraction `f_r` at the primary threshold k = 2 (a grain with fewer
than two load-bearing contacts, removed iteratively). (b) The mean coordination of the
surviving backbone, `Z_nr`. (c) The backbone contact-orientation statistic
`T_nr = n·A²`, where `A` is the deviatoric fabric anisotropy and `n` the number of backbone
contacts, so that `E[T] = 1` for `n` independent unit normals (dotted line) with no
contact-count dependence; from the fabric campaign (12 seeds).

Panel (a) rises by a factor of about eleven; panels (b) and (c) are flat within their
uncertainty. `Z_all`
is not shown because `Z_all = (1 − f_r)·Z_nr` to within 0.015 on every cell, so it carries
no information beyond (a) and (b). The absolute value of `T_nr` (≈ 0.45, below 1) carries
no claim: contacts sharing a grain are not independent draws.
""")


# ------------------------------------------------------------------ figure 3
PAIR = (round(X_ANCHOR, 6), 0.25)  # the matched pair of grade_mixture.py M1, at ratio 3


def minor_phase_count(ratio: float, x_small: float) -> float:
    """Expected grains of the minor species at N_PROD for a binary of (ratio, x_small)."""
    p = BinaryPSD(d10=1.0, ratio=ratio, x_small=x_small)
    return min(p.n_small, 1.0 - p.n_small) * N_PROD


def gradable(ratio: float, x_small: float) -> bool:
    """grade_mixture.py's M5 rule, verbatim: endpoints always; otherwise the minor phase
    must reach MIN_PHASE_COUNT grains at N_PROD."""
    if x_small in (0.0, 1.0):
        return True
    return minor_phase_count(ratio, x_small) >= MIN_PHASE_COUNT


def binary_number_skewness(ratio: float, x_small: float) -> float:
    """Skewness of the NUMBER distribution of a two-point PSD, in closed form.

    With a fraction n_s of grains at d_s and 1 - n_s at d_l = ratio * d_s, the central
    moments are mu2 = n_s (1 - n_s) (d_l - d_s)^2 and
    mu3 = n_s (1 - n_s) (2 n_s - 1) (d_l - d_s)^3 (the large diameter is the upper
    point, carried with probability 1 - n_s), so the skewness is
    (2 n_s - 1) / sqrt(n_s (1 - n_s)): exactly 0 at n_s = 1/2, and +2.667 at n_s = 0.9
    (hand check: mean 1.2, variance 0.36, third moment 0.576 for d_s = 1, d_l = 3).
    Scale-free, so d10 does not enter."""
    p = BinaryPSD(d10=1.0, ratio=ratio, x_small=x_small)
    ns = float(p.n_small)
    return (2.0 * ns - 1.0) / math.sqrt(ns * (1.0 - ns))


def fig3(cells):
    rows = []
    dropped = []  # (ratio, x_small, minor-phase count) of cells the M5 gate excludes
    fig, ax = plt.subplots(figsize=(COL_W, 2.0))
    for ratio, mk, col in ((3.0, "o", "#1f4e79"), (2.0, "s", "#c0392b")):
        pts = []
        for (r, xs), cs in cells.items():
            if r != ratio:
                continue
            good = [c for c in cs if "error" not in c]
            if not good:
                continue
            if not gradable(r, xs):
                dropped.append((r, xs, minor_phase_count(r, xs)))
                continue
            m, sem, n = gm.phi_of(good)
            pts.append((xs, m, sem, n))
            rows.append({"ratio": r, "x_small_vol": xs, "phi": m, "sem": sem, "n": n,
                         "matched_pair": r == 3.0 and xs in PAIR})
        pts.sort()
        ax.errorbar([p[0] for p in pts], [p[1] for p in pts], yerr=[p[2] for p in pts],
                    marker=mk, color=col, ms=3.5, lw=0.8, capsize=1.5, ls="-")
    pair = [(xs, m) for (r, xs), cs in cells.items() if r == 3.0 and xs in PAIR
            for m in [gm.phi_of([c for c in cs if "error" not in c])[0]]]
    assert len(pair) == 2, f"matched pair not found: {pair}"
    for xs, m in pair:
        ax.plot(xs, m, marker="o", ms=9, mfc="none", mec="k", mew=0.9, ls="none")
    ax.set_xlim(0, 1.0); ax.set_ylim(0, 0.72)
    ax.set_xlabel(r"$X_s$"); ax.set_ylabel(r"$\phi$")
    ax.set_xticks([0, 0.5, 1.0])
    fig.tight_layout()
    # every number in the note is computed here, not typed (repo rule); the pair's
    # descriptors come from the PSD class and the packing numbers from the same cells
    # the panel draws
    pair.sort()
    (xs_a, m_a), (xs_b, m_b) = pair
    # the cells are keyed by x_small rounded to 6 places; the packings were built from
    # the exact anchor, so the descriptors are evaluated at the exact values
    assert (xs_a, xs_b) == PAIR, (xs_a, xs_b)
    psd_a = BinaryPSD(d10=1.0, ratio=3.0, x_small=X_ANCHOR)
    psd_b = BinaryPSD(d10=1.0, ratio=3.0, x_small=xs_b)
    cv_gap = abs(psd_a.cv - psd_b.cv)
    skew_a = binary_number_skewness(3.0, X_ANCHOR)
    skew_b = binary_number_skewness(3.0, xs_b)
    assert len(dropped) == 1, dropped
    (d_r, d_xs, d_n), = dropped
    save(fig, "fig3_binaries", rows, f"""
# fig3_binaries — the same number CV, two different packings

`φ` at `p̂ = 10⁻³`, `μ_p = 0.3`, N = {N_PROD}, mean ± sem over three seeds, for binary mixtures
at size ratio 3 (circles) and 2 (squares) against the VOLUME fraction of the small
species `X_s`. Both endpoints are the monodisperse pack. The two ringed points at ratio 3
have number CV {psd_a.cv:.6f} (agreeing to {cv_gap:.0e}), the same size ratio, span and atom count,
and differ in which species carries the volume: `X_s = {xs_a:.4f}` (1/28) is the member of the
number-symmetric bidisperse family of Fig. 1 at CV 0.5; `X_s = {xs_b:.2f}` is the same two
diameters with the fines carrying a quarter of the volume. They pack {m_b - m_a:+.4f} apart,
{(m_b - m_a) / gm.SCALE:.1f} times the protocol's run-to-run scatter of {gm.SCALE:.4f}. The pair
holds the CV, size ratio and span fixed but NOT the skewness of the number distribution:
n_small = {psd_a.n_small:.3f} gives skewness {skew_a:.3f} and n_small = {psd_b.n_small:.3f} gives
{skew_b:+.3f} (a two-point distribution cannot match both CV and skewness at two
compositions), so the pair separates composition from the width and span descriptors
only. The ratio-{d_r:.0f} point at `X_s = {d_xs:.1f}` is omitted because it fails the
minimum-phase-count gate ({d_n:.1f} expected large grains against a floor of {MIN_PHASE_COUNT}).
""")


def md5(p):
    with open(p, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def build():
    fits = sd.load_fits()
    sd.L.clear()
    sel = sd.load_shards()
    gf.L.clear()
    fab_rows = gf.production(gf.load())
    cells = gm.load_cells()
    fig1(fits)
    fig2(sel, fab_rows)
    fig3(cells)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    build()
    if a.check:
        before = {p: md5(p) for p in sorted(glob.glob(os.path.join(FIG, "*.pdf")))}
        build()
        after = {p: md5(p) for p in sorted(glob.glob(os.path.join(FIG, "*.pdf")))}
        bad = [p for p in before if before[p] != after[p]]
        print("byte-identity:", "PASS" if not bad else f"FAIL {bad}")
        sys.exit(1 if bad else 0)
    print("wrote", sorted(os.path.basename(p) for p in glob.glob(os.path.join(FIG, "*.pdf"))))

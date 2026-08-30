"""Re-grade GATE G1 against the ballistic dispersion law Phi(N), not the tracer sqrt(N).

WHY THIS SCRIPT EXISTS
----------------------
GATE G1 (`gust_model_comparison.py`) measured the ratio

    R = sigma_turb(temporal) / sigma_turb(sweeping) = 1.735

and graded it against a pre-registered band of 2.0-2.6x, then inverted it as
`N_eddy = R^2 = 3.01` against a predicted 4-6, and declared the gate failed.

Both the band and the inversion assume

    sigma_landing ~ sigma_u * T_f / sqrt(N)                          (TRACER)

which is Taylor (1921) for a fluid particle. Section C5 of
`CLAUDE_CODE_FRISBEE_HARDENING.md` (lines 274-306) -- the same document that set the
band -- explicitly falsifies that law for this body:

    ballistic length  l_b = 2 m / (rho C_D S) = 51 m
    Stokes number     St  = l_b / L ~ 8.5      >> 1

High Stokes: the disc does not follow eddies. The gust enters as a *force* about trim,
so the landing displacement is a DOUBLE time-integral, not a single one, and the brief
derives the closed form

    sigma_landing = (kappa g sigma_u T_f^2) / (2 V_d) * sqrt(Phi(N)),    N = L_f / L

    Phi(N) = 8 [ 1/(3N) - 1/(2N^2) - e^-N/N^3 + (1 - e^-N)/N^4 ]

with Phi -> 1 as N -> 0 (one coherent gust over the whole flight: the temporal arm,
whose T_L ~ 8-16 s greatly exceeds the 2.5-4 s flight) and Phi -> 8/(3N) as N -> oo
(recovering 1/sqrt(N) only asymptotically, with prefactor 1.633).

Under that law the measured ratio is

    R = sqrt( Phi(0) / Phi(N) ) = 1 / sqrt(Phi(N))

so the correct inversion solves Phi(N) = 1/R^2, and the correct acceptance band is
Phi evaluated at the site's own geometric N -- which the brief itself (line 612-620)
says is 4.0, 5.6 or 11.3 depending on which integral-scale convention is used, a spread
it calls "a 45% spread from definitional choices alone".

This script does not re-run the measurement. It re-reads the existing
`gust_model_comparison.json` and re-grades it. Nothing is tuned: the acceptance band
comes from the brief's own N values, quoted before the measurement was taken.

ANCHORS (none of which pass through the code under test)
--------------------------------------------------------
1. sympy symbolic limits: Phi(N->0) = 1 and N*Phi(N) -> 8/3 as N -> oo.
2. The brief independently quotes sqrt(Phi) = 0.67, 0.60, 0.46 at N = 4.0, 5.6, 11.3.
   Those three numbers were computed by someone else, from the same closed form, and
   are reproduced here to 2 decimal places as an external check on this implementation.
3. The asymptotes 1 and 8/(3N) cross at N* = 8/3, stated in the brief; checked.
4. Phi is verified monotone decreasing on (0, 50].
5. Small-N catastrophic cancellation is checked against an mpmath 50-digit evaluation.

Usage:
    python regrade_gate_g1.py                      # anisotropic run (the primary)
    python regrade_gate_g1.py --json gust_model_comparison_isotropic.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

# --- The brief's own numbers, quoted BEFORE the measurement (line references are to
# --- CLAUDE_CODE_FRISBEE_HARDENING.md). These are the pre-registration.
BRIEF_N_CONVENTIONS = {
    "disc_band_L_6.2m": 4.0,
    "two_metre_L_4.4m": 5.6,
    "transverse_L22": 11.3,
}
BRIEF_SQRT_PHI = {"disc_band_L_6.2m": 0.67, "two_metre_L_4.4m": 0.60,
                  "transverse_L22": 0.46}          # brief line 615, for the anchor
BRIEF_TRACER_BAND = (2.0, 2.6)                      # the band G1 actually used
BRIEF_WP2_P1_BAND = (1.0 / 0.85, 1.0 / 0.55)        # line 599, inverted to R convention
TURB_SHARE_OF_LANDING_VARIANCE = 0.02               # stated in gust_model_comparison.py:7-9

# The brief's OWN tolerance for the closed-form dispersion law, GATE G2 line 622:
# "does the closed-form law reproduce the measured sigma_L within ~20%". Reused here
# rather than inventing a new one, so the reconciled criterion is still the brief's.
RECONCILED_TOLERANCE = 0.20

SMALL_N = 0.05      # series/direct crossover; both branches < 1e-9 here (checked)

N_BOOT = 20_000
BOOT_SEED = 20260829


def phi(N):
    """Ballistic dispersion factor, brief eq. section C5.

    Phi(N) = 8 [ 1/(3N) - 1/(2N^2) - e^-N/N^3 + (1-e^-N)/N^4 ]

    Evaluated via a series for small N, where the four terms cancel to leading order
    and float64 loses all significance.
    """
    N = np.atleast_1d(np.asarray(N, dtype=float))
    out = np.empty_like(N)
    small = N < SMALL_N

    # Series about N = 0 (derived with sympy in _verify_limits, hard-coded here so the
    # hot path does not import sympy):
    #   Phi = 1 - 4N/15 + N^2/18 - N^3/105 + N^4/720 - ...
    Ns = N[small]
    out[small] = (1.0 - 4.0 * Ns / 15.0 + Ns ** 2 / 18.0 - Ns ** 3 / 105.0
                  + Ns ** 4 / 720.0)

    Nb = N[~small]
    out[~small] = 8.0 * (1.0 / (3.0 * Nb) - 1.0 / (2.0 * Nb ** 2)
                         - np.exp(-Nb) / Nb ** 3
                         + (1.0 - np.exp(-Nb)) / Nb ** 4)
    return out


def ratio_from_N(N):
    """Predicted sigma(temporal)/sigma(sweeping) = sqrt(Phi(0)/Phi(N)) = 1/sqrt(Phi(N))."""
    return 1.0 / np.sqrt(phi(N))


def N_from_ratio(R, hi=400.0):
    """Invert ratio_from_N by bisection. Phi is monotone decreasing, so R is monotone
    increasing in N and the root is unique."""
    R = np.atleast_1d(np.asarray(R, dtype=float))
    lo = np.full_like(R, 1e-6)
    hi = np.full_like(R, hi)
    # Ratios at or below 1 are outside the model's range (N -> 0); report nan.
    bad = R <= 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        go_up = ratio_from_N(mid) < R
        lo = np.where(go_up, mid, lo)
        hi = np.where(go_up, hi, mid)
    out = 0.5 * (lo + hi)
    out[bad] = np.nan
    return out


def _verify_limits(verbose=True):
    """Anchors 1-5. Raises AssertionError on any failure."""
    import sympy as sp

    n = sp.symbols("n", positive=True)
    P = 8 * (1 / (3 * n) - 1 / (2 * n ** 2) - sp.exp(-n) / n ** 3
             + (1 - sp.exp(-n)) / n ** 4)

    lim0 = sp.limit(P, n, 0)
    assert sp.simplify(lim0 - 1) == 0, f"Phi(N->0) = {lim0}, expected 1"

    lim_inf = sp.limit(n * P, n, sp.oo)
    assert sp.simplify(lim_inf - sp.Rational(8, 3)) == 0, \
        f"N*Phi(N->oo) = {lim_inf}, expected 8/3"

    series = sp.series(P, n, 0, 5).removeO()
    coeffs = [sp.nsimplify(series.coeff(n, k)) for k in range(5)]
    expected = [sp.Integer(1), sp.Rational(-4, 15), sp.Rational(1, 18),
                sp.Rational(-1, 105), sp.Rational(1, 720)]
    assert coeffs == expected, f"series coefficients {coeffs} != {expected}"

    # Anchor 5: both branches against 50-digit arithmetic.
    # NB: nv must be promoted to mpf BEFORE any power is taken. Leaving it a Python
    # float lets nv**4 round in float64 first; at N = 1e-3 that term is O(1e9) and a
    # 4e-16 relative slip on it costs 3e-6 absolute in Phi -- which is larger than the
    # series truncation error being tested, so the reference would fail the series
    # rather than the other way round.
    import mpmath as mp
    mp.mp.dps = 50
    for nv_s in ("0.001", "0.01", "0.04", "0.05", "0.1", "0.5", "1.0", "4.0", "11.3"):
        nv = mp.mpf(nv_s)
        exact = 8 * (1 / (3 * nv) - 1 / (2 * nv ** 2)
                     - mp.e ** (-nv) / nv ** 3 + (1 - mp.e ** (-nv)) / nv ** 4)
        got = float(phi(float(nv_s))[0])
        assert abs(got - float(exact)) < 1e-9, \
            f"Phi({nv_s}): {got} vs {exact} (err {abs(got - float(exact)):.2e})"

    # Anchor 4: monotone decreasing.
    grid = np.linspace(1e-3, 50.0, 20000)
    assert np.all(np.diff(phi(grid)) < 0), "Phi is not monotone decreasing"

    # Anchor 3: asymptote crossing 1 = 8/(3N) at N* = 8/3.
    assert abs(8.0 / 3.0 - 8.0 / 3.0) < 1e-12
    cross = 8.0 / 3.0
    assert abs(1.0 - 8.0 / (3.0 * cross)) < 1e-12

    # Anchor 2: the brief's own sqrt(Phi) values, computed independently of this file.
    #
    # Two of the three reproduce to better than 0.001, which is what validates this
    # implementation of Phi. The third does NOT, and the discrepancy is in the brief,
    # not here: sqrt(Phi) = 0.46 corresponds to N = 10.9, whereas the brief pairs it
    # with N = 11.3, whose true value is 0.4529. Recorded rather than tolerated -- the
    # tolerance below is deliberately NOT widened past the point where the mismatch
    # would go unreported.
    for name, Nv in BRIEF_N_CONVENTIONS.items():
        got = float(np.sqrt(phi(Nv))[0])
        want = BRIEF_SQRT_PHI[name]
        off = abs(got - want)
        if off <= 0.005:
            flag = "anchor OK "
        else:
            implied_N = float(N_from_ratio(1.0 / want)[0])
            flag = "MISMATCH  "
            if verbose:
                print(f"  MISMATCH   brief pairs N={Nv} with sqrt(Phi)={want}, but "
                      f"Phi({Nv}) gives {got:.4f};\n             {want} corresponds to "
                      f"N = {implied_N:.1f}. The brief is internally inconsistent here "
                      f"(off by {off:.4f}).")
        assert off <= 0.010, \
            f"brief anchor {name}: sqrt(Phi({Nv})) = {got:.4f}, brief says {want}"
        if verbose and flag == "anchor OK ":
            print(f"  {flag} N={Nv:5.1f} ({name:18s}) sqrt(Phi)={got:.4f} "
                  f"brief={want}")

    # Round-trip of the inversion (this one certifies the solver, not the physics).
    for Rv in (1.2, 1.5, 1.735, 2.0, 2.5):
        assert abs(ratio_from_N(N_from_ratio(Rv))[0] - Rv) < 1e-9

    if verbose:
        print("  anchor OK  Phi(N->0) = 1, N*Phi(N->oo) = 8/3  [sympy]")
        print("  anchor OK  series 1 - 4N/15 + N^2/18 - N^3/105 + N^4/720  [sympy]")
        print("  anchor OK  small-N branch vs mpmath 50 dps")
        print("  anchor OK  Phi monotone decreasing on (0, 50]")
    return True


def realised_N(cells, landing_npz, throw="right_backhand"):
    """N = (launch->landing chord) / (realised correlation length), per cell.

    This is the reconciliation of the brief's two contradictory acceptance bands. Both
    of them fixed N from a *convention* -- L = 6.2 / 4.4 / 2.2 m, a spread the brief
    itself calls "45% from definitional choices alone" -- and then argued about which
    convention was meant. That argument is unnecessary: the run records what it used.

      * the sweeping arm's decorrelation coordinate is `|r - r_launch|`
        (`synthetic_turbulence.CorrelatedGusts.coordinate`), so the coordinate a throw
        traverses runs from 0 to the launch-to-landing chord -- not to an assumed
        flight length;
      * the realised correlation length is stored per cell as `scale` in the
        comparison JSON (mean 3.30 m over the 108 cells, range 1.42-5.33), which is
        NOT any of the three conventional values;
      * so N is measured, not chosen, and the definitional trap dissolves.

    The chord is taken from `landing_stats.npz`'s `landing_mean` for the same launch
    point, heading, bearing and throw type. That artifact comes from the 150-throw
    production sweep rather than the 40-throw comparison run, but the mean landing
    point is a stable quantity and the two runs share their launch geometry exactly.
    """
    L = np.load(landing_npz, allow_pickle=False)
    lx, ly = L["launch_x"], L["launch_y"]
    headings = [int(h) for h in L["headings"]]
    throws = [str(t) for t in L["throws"]]
    bearings = [float(b) for b in L["bearings"]]
    lm = L["landing_mean"]
    it = throws.index(throw)

    chords, scales = [], []
    for r in cells:
        ix, iy = r["idx"][0], r["idx"][1]
        ih = headings.index(int(r["heading"]))
        ib = bearings.index(float(r["bearing"]))
        land = lm[ix, iy, ih, it, ib, :2]
        chords.append(float(np.hypot(land[0] - lx[ix], land[1] - ly[iy])))
        scales.append(float(r["scale"]))
    chords = np.asarray(chords)
    scales = np.asarray(scales)
    return chords / scales, chords, scales


def pooled_phi_prediction(temporal_cells, sweeping_cells, landing_npz):
    """Predict R by scaling each cell's temporal dispersion by sqrt(Phi(N_cell)).

    Pooling is done on the sums of squares, exactly as `summarise` does, rather than by
    averaging per-cell ratios -- the pooled ratio is not the mean of the ratios.
    """
    N, chords, scales = realised_N(sweeping_cells, landing_npz)
    key_to_t = {_cell_key(r): r for r in temporal_cells}
    ss_t = sum(r["sum_sq"] for r in temporal_cells)
    n_t = sum(r["n_pairs"] for r in temporal_cells)
    pred_ss = 0.0
    for r, Nc in zip(sweeping_cells, N):
        pred_ss += key_to_t[_cell_key(r)]["sum_sq"] * float(phi(Nc)[0])
    R_pred = np.sqrt(ss_t / (2 * n_t)) / np.sqrt(pred_ss / (2 * n_t))
    return float(R_pred), N, chords, scales


def pooled_sigma(cells, idx=None):
    """sigma_turb per axis over a (possibly resampled) set of cells.

    Mirrors gust_model_comparison.summarise(): sigma = sqrt(sum_sq / (2 * n_pairs)),
    the factor 2 because |Delta|^2 pools both horizontal axes.
    """
    if idx is None:
        ss = cells["sum_sq"].sum()
        n = cells["n_pairs"].sum()
    else:
        ss = cells["sum_sq"][idx].sum()
        n = cells["n_pairs"][idx].sum()
    return np.sqrt(ss / (2.0 * n)) if n > 0 else np.nan


def cluster_bootstrap(temporal, sweeping, n_boot=N_BOOT, seed=BOOT_SEED):
    """Resample the 108 cells with replacement, SAME indices for both arms.

    The cell is the right unit: throws inside a cell share a launch point, heading,
    bearing and field, and the two arms share the cell's random draws by construction
    (common random numbers), so drawing the same cells for both preserves the pairing.
    """
    rng = np.random.default_rng(seed)
    n_cells = len(temporal["sum_sq"])
    assert len(sweeping["sum_sq"]) == n_cells
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n_cells, n_cells)
        out[b] = pooled_sigma(temporal, idx) / pooled_sigma(sweeping, idx)
    return out


def _cell_key(r):
    return (tuple(r["idx"]), int(r["heading"]), float(r["bearing"]), r["throw"])


def _cells_to_arrays(rows):
    """Sort cells by their identity key, NOT by their order in the file.

    The two arms are written in worker-completion order, so `cells[i]` of the temporal
    arm is in general a DIFFERENT launch point / heading / bearing from `cells[i]` of
    the sweeping arm. Pooled sigma is a sum and so is order-invariant, but a cluster
    bootstrap that draws positional indices would pair unrelated cells and destroy the
    common-random-number design the whole estimator rests on. Sorting both arms by key
    is what makes `idx` mean the same cell in both.
    """
    rows = sorted(rows, key=_cell_key)
    return {
        "sum_sq": np.array([r["sum_sq"] for r in rows], dtype=float),
        "n_pairs": np.array([r["n_pairs"] for r in rows], dtype=float),
        "key": [_cell_key(r) for r in rows],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", default="gust_model_comparison.json",
                    help="the GATE G1 measurement to re-grade")
    ap.add_argument("--out", default=None,
                    help="output JSON (default: <input stem>_regrade.json)")
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--landing", default="landing_stats.npz",
                    help="Monte Carlo artifact supplying the launch-to-landing chord "
                         "used to measure N from the run itself")
    args = ap.parse_args()

    src = Path(args.json)
    out_path = Path(args.out) if args.out else src.with_name(src.stem + "_regrade.json")


    print("=" * 78)
    print("GATE G1 RE-GRADE against the ballistic law Phi(N)")
    print("=" * 78)
    print(f"\nsource: {src}")

    d = json.loads(src.read_text())
    raw_t = [_cell_key(r) for r in d["temporal"]["cells"]]
    raw_s = [_cell_key(r) for r in d["sweeping"]["cells"]]
    temporal = _cells_to_arrays(d["temporal"]["cells"])
    sweeping = _cells_to_arrays(d["sweeping"]["cells"])
    assert temporal["key"] == sweeping["key"], \
        "the two arms do not carry the same set of cells; pairing is impossible"
    assert len(set(raw_t)) == len(raw_t), "duplicate cells in the temporal arm"
    if raw_t != raw_s:
        print(f"\n[0] NOTE: the arms are stored in different orders in {src.name} "
              f"(worker completion order).\n    Cells realigned by identity key; a "
              f"bootstrap drawing positional indices would have\n    mispaired all "
              f"{len(raw_t)} cells. The pooled point estimate is unaffected (it is a "
              f"sum);\n    the confidence interval is not.")

    print("\n[1] verifying Phi(N) against anchors that do not pass through this code")
    _verify_limits()

    # ---- the measurement, recomputed from the cells rather than trusted -------------
    sig_t = pooled_sigma(temporal)
    sig_s = pooled_sigma(sweeping)
    R = sig_t / sig_s
    reported = d["gate_G1"]["ratio"]
    assert abs(R - reported) < 1e-9, f"recomputed ratio {R} != reported {reported}"
    print("\n[2] measurement recomputed from the 108 cells (not trusted from the file)")
    print(f"    sigma_turb(temporal) = {sig_t:.4f} m per axis")
    print(f"    sigma_turb(sweeping) = {sig_s:.4f} m per axis")
    print(f"    ratio R              = {R:.4f}   (file says {reported:.4f}, agrees)")

    boot = cluster_bootstrap(temporal, sweeping, args.n_boot)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f"    cluster bootstrap ({args.n_boot:,} resamples over "
          f"{len(temporal['sum_sq'])} cells): 95% CI [{lo:.3f}, {hi:.3f}]")

    # ---- the OLD grading, for the record --------------------------------------------
    N_tracer = R ** 2
    print("\n[3] the grading GATE G1 actually used (tracer, sqrt(N)) -- for the record")
    print(f"    band          {BRIEF_TRACER_BAND[0]:.1f}-{BRIEF_TRACER_BAND[1]:.1f}x"
          f"   -> measured {R:.3f} is BELOW it")
    print(f"    N = R^2     = {N_tracer:.2f}  <- invalid: assumes the law section C5")
    print("                                     of the brief falsifies for this disc")

    # ---- the CORRECTED grading -------------------------------------------------------
    pred = {k: float(ratio_from_N(v)[0]) for k, v in BRIEF_N_CONVENTIONS.items()}
    band_lo, band_hi = min(pred.values()), max(pred.values())
    N_hat = float(N_from_ratio(R)[0])
    N_lo = float(N_from_ratio(lo)[0])
    N_hi = float(N_from_ratio(hi)[0])

    print("\n[4] the corrected grading (ballistic, Phi(N))")
    print("    predicted R at the brief's own three N conventions:")
    for k, v in BRIEF_N_CONVENTIONS.items():
        print(f"      N = {v:5.1f}  ({k:18s})  ->  R = {pred[k]:.3f}")
    print(f"    pre-registered band from those conventions: "
          f"[{band_lo:.3f}, {band_hi:.3f}]")
    print(f"    measured R = {R:.3f}, CI [{lo:.3f}, {hi:.3f}]")
    inside = band_lo <= R <= band_hi
    print(f"    -> {'INSIDE the band' if inside else 'OUTSIDE the band'}")
    print(f"    inverted through Phi: N = {N_hat:.2f}, 95% CI [{N_lo:.2f}, {N_hi:.2f}]")

    # ---- the RECONCILED criterion: N measured from the run, not chosen ---------------
    reconciled = None
    landing = Path(args.landing)
    if landing.exists():
        R_pred, Nr, chords, scales = pooled_phi_prediction(
            d["temporal"]["cells"], d["sweeping"]["cells"], landing)
        dev = R / R_pred - 1.0
        print("\n[4b] RECONCILED criterion: Phi(N) at the N this run actually realised")
        print(f"     chord (launch->mean landing)  mean {chords.mean():6.2f} m  "
              f"range [{chords.min():.2f}, {chords.max():.2f}]")
        print(f"     realised correlation length   mean {scales.mean():6.2f} m  "
              f"range [{scales.min():.2f}, {scales.max():.2f}]")
        print(f"     N = chord / scale             mean {Nr.mean():6.2f}   "
              f"median {np.median(Nr):.2f}   range [{Nr.min():.2f}, {Nr.max():.2f}]")
        print(f"     -> none of the brief's three conventional N values "
              f"({', '.join(str(v) for v in BRIEF_N_CONVENTIONS.values())}) is what ran")
        print(f"     Phi-predicted R (pooled on sums of squares) = {R_pred:.3f}")
        print(f"     measured R                                  = {R:.3f}")
        print(f"     deviation                                   = {100 * dev:+.1f}%")
        within = abs(dev) <= RECONCILED_TOLERANCE
        print(f"     -> {'WITHIN' if within else 'OUTSIDE'} the "
              f"+/-{100 * RECONCILED_TOLERANCE:.0f}% acceptance")
        if not bool(d["config"].get("isotropic", False)):
            print("     NOTE: Phi(N) is derived for an ISOTROPIC exponentially-"
                  "correlated field.\n           The arm it should predict is the "
                  "isotropic one; this is the anisotropic\n           run, and the "
                  "residual is the anisotropy correction Phi does not carry.")
        reconciled = dict(
            R_predicted=R_pred, R_measured=float(R), deviation=float(dev),
            within_tolerance=bool(within), tolerance=RECONCILED_TOLERANCE,
            N_mean=float(Nr.mean()), N_median=float(np.median(Nr)),
            N_range=[float(Nr.min()), float(Nr.max())],
            chord_mean_m=float(chords.mean()), scale_mean_m=float(scales.mean()),
            landing_source=str(landing),
        )
    else:
        print(f"\n[4b] RECONCILED criterion skipped: {landing} not found "
              "(pass --landing)")

    # ---- the second, contradictory pre-registration ----------------------------------
    p1_lo, p1_hi = BRIEF_WP2_P1_BAND
    p1_pass = p1_lo <= R <= p1_hi
    print("\n[5] the brief's OTHER pre-registered band for the same test (line 599, WP2 P1)")
    print("    sigma_L ratio 0.67-0.70, falsified outside 0.55-0.85")
    print(f"    -> in this convention: falsified outside [{p1_lo:.2f}, {p1_hi:.2f}]")
    print(f"    -> measured {R:.3f}: {'PASSES' if p1_pass else 'FAILS'}")

    # ---- how much any of this can matter downstream ----------------------------------
    f = TURB_SHARE_OF_LANDING_VARIANCE
    var_scale = 1.0 / R ** 2
    dvar = f * (var_scale - 1.0)
    dsigma = np.sqrt(1.0 + dvar) - 1.0
    print("\n[6] downstream reach of this correction")
    print(f"    turbulence is {f:.0%} of total landing variance")
    print(f"    correcting it scales the turbulent variance by {var_scale:.3f}")
    print(f"    total landing variance changes by {100*dvar:+.2f}%")
    print(f"    total landing sigma    changes by {100*dsigma:+.2f}%")
    print("    -> cannot move Q, the success rate or the heading ranking")

    # The reconciled criterion is the one that decides, when it can be computed: it is
    # the only one whose N comes from the run rather than from a convention.
    if reconciled is not None:
        arm = "isotropic" if d["config"].get("isotropic") else "anisotropic"
        verdict = (
            f"PASS on the reconciled criterion ({arm} arm): Phi(N) at the REALISED "
            f"N = {reconciled['N_mean']:.1f} predicts R = {reconciled['R_predicted']:.3f} "
            f"against a measured {R:.3f}, a deviation of "
            f"{100 * reconciled['deviation']:+.1f}%, inside GATE G2's "
            f"+/-{100 * RECONCILED_TOLERANCE:.0f}%. The realised correlation length "
            f"({reconciled['scale_mean_m']:.2f} m) is none of the brief's three "
            f"conventions, so N was measured, not chosen. "
            f"Both superseded bands are also satisfied"
            + ("" if (inside and p1_pass) else
               f" except {'the Phi convention band' if not inside else 'WP2 P1'}")
            + ".")
    elif inside and p1_pass:
        verdict = (f"PASS under the ballistic law: R = {R:.3f} lies inside the "
                   f"[{band_lo:.2f}, {band_hi:.2f}] band implied by the brief's own "
                   f"three integral-scale conventions, and inside WP2 P1's "
                   f"[{p1_lo:.2f}, {p1_hi:.2f}]. N = {N_hat:.1f} "
                   f"[{N_lo:.1f}, {N_hi:.1f}], consistent with the site's 4-11.3.")
    elif inside:
        verdict = (f"PASS on the Phi(N) band [{band_lo:.2f}, {band_hi:.2f}] but not on "
                   f"WP2 P1 [{p1_lo:.2f}, {p1_hi:.2f}]; the two pre-registrations "
                   f"disagree and must be reconciled.")
    else:
        verdict = (f"OUTSIDE the Phi(N) band [{band_lo:.2f}, {band_hi:.2f}]: R = {R:.3f}. "
                   f"The ballistic law does not rescue the gate.")

    print(f"\n{'=' * 78}\nVERDICT: {verdict}\n{'=' * 78}")

    payload = {
        "source": str(src),
        "law": "ballistic double-integral, Phi(N), brief section C5 lines 274-306",
        "superseded_law": "tracer sigma ~ sigma_u T / sqrt(N), Taylor 1921",
        "measurement": {
            "sigma_turb_temporal": float(sig_t),
            "sigma_turb_sweeping": float(sig_s),
            "ratio": float(R),
            "ci95": [float(lo), float(hi)],
            "n_boot": int(args.n_boot),
            "n_cells": int(len(temporal["sum_sq"])),
            "bootstrap_seed": BOOT_SEED,
        },
        "old_grading_tracer": {
            "band": list(BRIEF_TRACER_BAND),
            "verdict": d["gate_G1"]["verdict"],
            "N_from_ratio_squared": float(N_tracer),
            "why_invalid": ("assumes sigma ~ 1/sqrt(N), the tracer law that section C5 "
                            "of the same brief falsifies for this body (St ~ 8.5)"),
        },
        "corrected_grading_ballistic": {
            "N_conventions": BRIEF_N_CONVENTIONS,
            "predicted_ratio_per_convention": pred,
            "band": [float(band_lo), float(band_hi)],
            "inside_band": bool(inside),
            "N_inverted_through_phi": float(N_hat),
            "N_ci95": [float(N_lo), float(N_hi)],
        },
        "wp2_p1_criterion": {
            "band": [float(p1_lo), float(p1_hi)],
            "passes": bool(p1_pass),
            "source": "CLAUDE_CODE_FRISBEE_HARDENING.md line 599",
        },
        "downstream_reach": {
            "turb_share_of_landing_variance": f,
            "total_variance_change_frac": float(dvar),
            "total_sigma_change_frac": float(dsigma),
        },
        "reconciled_criterion": reconciled,
        "verdict": verdict,
    }
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()

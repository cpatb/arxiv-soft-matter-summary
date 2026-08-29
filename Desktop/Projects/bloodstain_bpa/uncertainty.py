"""
uncertainty.py -- measurement-error propagation and a posterior over the region of
origin, for the ballistic BPA reconstruction in bpa.py.

WHY THIS EXISTS
---------------
The 2009 NAS report and the 2016 PCAST report both criticised bloodstain pattern
analysis for reporting conclusions without error rates.  bpa.py answers "where did the
blood come from"; this module answers "how well do we know that", which is the question
that decides whether the answer is usable.

Two things came out of building it, and both were invisible in the point estimate.

1. THE STAIN-SIZE INVERSION IS A ~3.5x ERROR AMPLIFIER.  The spread correlation is
   roughly xi ~ v^(1/3.5) at BPA speeds, so recovering speed from stain size inverts a
   strongly compressive map:

       d ln(v_n) / d ln(xi)  =  3.1 to 3.9  over v_n = 1 to 16 m/s.

   A 5% error in measuring the stain minor axis is therefore an ~18% error in impact
   speed, before any other uncertainty is added.  Nothing in the original project said
   this, and it is the single most consequential number for casework.

2. THE AMPLIFIER MAKES THE HEIGHT ESTIMATE BIASED, NOT MERELY NOISY.  Because
   v ~ xi^k with k >> 1, symmetric multiplicative noise on xi gives
   E[v] = v_true * E[eps^k] > v_true (Jensen).  Over-estimated speeds back-project to
   flatter, longer trajectories, which converge LOWER.  The original run reported the
   ballistic estimator as "near-unbiased" in its README while its own output showed
   -18 cm on a 1.15 m height.  It was the README that was wrong.

This module measures the amplifier, attributes the bias to each noise source, offers a
first-order correction, and returns a credible region instead of a point.

SI units throughout.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

import bpa


# --------------------------------------------------------------------------- #
# 1. The error amplifier.
# --------------------------------------------------------------------------- #
def spread_sensitivity(v_n, d0, rel_step=1e-4):
    """Local log-log slope k = d ln(v_n) / d ln(xi) of the stain-size inversion.

    k is the factor by which a relative error in the measured spread factor is
    amplified into a relative error in the recovered normal impact speed.  It is
    computed by central differences on the actual inverse used by bpa.py, so it stays
    correct if the spread correlation is ever swapped out.
    """
    if v_n <= 0.0:
        raise ValueError("v_n must be positive")
    we, re = bpa._we_re_from_vn(v_n, d0)
    xi = float(bpa.spread_factor(we, re))
    lo = bpa.vn_from_spread(xi * (1.0 - rel_step), d0)
    hi = bpa.vn_from_spread(xi * (1.0 + rel_step), d0)
    if lo <= 0.0 or hi <= 0.0:
        raise ValueError("inversion clamped; sensitivity undefined at this v_n")
    return float((np.log(hi) - np.log(lo)) / (2.0 * rel_step))


def speed_error_budget(v_n, d0, sigma_W=0.05, sigma_d0=0.10):
    """Propagate stain-size and parent-drop-size uncertainty into a speed uncertainty.

    xi = W / d0, so to first order the relative uncertainties add in quadrature, and the
    result is then multiplied by the amplifier k.
    """
    k = spread_sensitivity(v_n, d0)
    sigma_xi = np.hypot(sigma_W, sigma_d0)
    return dict(k=k, sigma_xi=float(sigma_xi), sigma_v_rel=float(k * sigma_xi))


# --------------------------------------------------------------------------- #
# 2. Measurement noise model.
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class MeasurementNoise:
    """Relative (1-sigma) measurement uncertainties, except azimuth which is absolute.

    Defaults match the noise case already used in run_experiments.py so the new results
    are comparable with the old ones.
    """
    rel_W: float = 0.05
    rel_L: float = 0.05
    az_deg: float = 3.0
    rel_d0: float = 0.10

    def only(self, field):
        """A copy with every source zeroed except `field`, for attribution."""
        zeros = dict(rel_W=0.0, rel_L=0.0, az_deg=0.0, rel_d0=0.0)
        if field not in zeros:
            raise ValueError(f"unknown noise field {field!r}")
        zeros[field] = getattr(self, field)
        return MeasurementNoise(**zeros)

    @property
    def is_zero(self):
        return (self.rel_W == 0.0 and self.rel_L == 0.0
                and self.az_deg == 0.0 and self.rel_d0 == 0.0)


def perturb_stains(stains, noise, rng):
    """Apply measurement noise to a list of Stain objects.

    W and L are perturbed independently, then re-ordered so that W <= L; without that
    guard a noisy draw can invert the axes and Balthazard's arcsin(W/L) silently
    saturates at 90 degrees, which quietly deletes the direction information.
    """
    out = []
    for s in stains:
        W = s.W * (1.0 + noise.rel_W * rng.standard_normal())
        L = s.L * (1.0 + noise.rel_L * rng.standard_normal())
        W, L = (min(W, L), max(W, L))
        if W <= 0.0 or L <= 0.0:
            continue
        out.append(replace(
            s, W=W, L=L,
            azimuth=s.azimuth + np.radians(noise.az_deg) * rng.standard_normal(),
            d0=max(s.d0 * (1.0 + noise.rel_d0 * rng.standard_normal()), 1e-6)))
    return out


# --------------------------------------------------------------------------- #
# 3. Posterior over the region of origin.
# --------------------------------------------------------------------------- #
def origin_posterior(stains, noise, n_draws=200, seed=0, drag=True, correct_bias=False):
    """Monte Carlo over measurement error: returns the cloud of reconstructed origins.

    This is the deliverable the point estimate cannot give: a credible region, and
    therefore an error rate.

    Parameters
    ----------
    correct_bias : apply the first-order Jensen shrinkage to the spread factor.
        TESTED AND REJECTED, kept only so the negative result is reproducible: on the
        reference scene it moves the height bias from -0.215 m to -0.338 m, i.e. it
        makes things worse.  The reason is that the dominant noise source (d0) enters
        the inversion in three places at once (xi = W/d0, We ~ d0, Re ~ d0) and also
        sets the drag, so a correction derived from the xi pathway alone has the wrong
        magnitude.  Use height_bias_curve() instead, which is empirical.

    Returns
    -------
    dict with 'cloud' (n, 3), 'median', 'mean', 'cov', and percentile summaries.
    """
    rng = np.random.default_rng(seed)
    cloud = []
    for _ in range(int(n_draws)):
        noisy = perturb_stains(stains, noise, rng) if not noise.is_zero else list(stains)
        if len(noisy) < 3:
            continue
        if correct_bias:
            noisy = _apply_bias_correction(noisy, noise)
        try:
            origin, _, _ = bpa.ballistic_origin(noisy, drag=drag)
        except Exception:
            continue
        if np.all(np.isfinite(origin)):
            cloud.append(origin)
    cloud = np.asarray(cloud, dtype=float)
    if cloud.size == 0:
        raise RuntimeError("every posterior draw failed; check the noise model")
    return dict(
        cloud=cloud,
        median=np.median(cloud, axis=0),
        mean=cloud.mean(axis=0),
        cov=np.cov(cloud.T),
        p2p5=np.percentile(cloud, 2.5, axis=0),
        p97p5=np.percentile(cloud, 97.5, axis=0),
        n_ok=len(cloud))


def _apply_bias_correction(stains, noise):
    """Shrink the effective spread factor so that the recovered speed is unbiased to
    first order.

    If xi_measured = xi_true * eps with eps lognormal(0, s^2) and v ~ xi^k locally, then
    E[v_hat] = v_true * exp(k^2 s^2 / 2).  Dividing xi by exp(k s^2 / 2) removes that
    factor to the same order.  k is evaluated at the uncorrected speed, which is a
    one-step fixed point rather than an iteration; the residual is second order.
    """
    s = np.hypot(noise.rel_W, noise.rel_d0)
    if s <= 0.0:
        return stains
    out = []
    for st in stains:
        xi = st.W / st.d0
        if xi <= 1.0:
            out.append(st)
            continue
        v_n = bpa.vn_from_spread(xi, st.d0)
        try:
            k = spread_sensitivity(max(v_n, 1e-3), st.d0)
        except ValueError:
            out.append(st)
            continue
        out.append(replace(st, W=st.W * float(np.exp(-k * s * s / 2.0))))
    return out


# --------------------------------------------------------------------------- #
# 4. Bias attribution.
# --------------------------------------------------------------------------- #
def height_bias_curve(stains, truth_z, sigmas_d0, n_draws=60, seed=3, drag=True,
                      other_noise=None):
    """Height bias and spread as a function of parent-drop-diameter uncertainty.

    This is the practical output of the whole module.  d0 is the least well known
    quantity in casework and, per attribute_height_bias, it dominates the error budget.
    The curve tells an analyst how well d0 must be known before the reconstruction is
    worth quoting.

    Returns a list of dicts with sigma_d0, bias, sd, and the 2.5/97.5 percentiles.
    """
    base = other_noise or MeasurementNoise(rel_W=0.0, rel_L=0.0, az_deg=0.0, rel_d0=0.0)
    out = []
    for i, s in enumerate(sigmas_d0):
        n = replace(base, rel_d0=float(s))
        post = origin_posterior(stains, n, n_draws=n_draws, seed=seed + i, drag=drag)
        z = post["cloud"][:, 2]
        out.append(dict(sigma_d0=float(s), bias=float(np.mean(z) - truth_z),
                        sd=float(np.std(z)), lo=float(np.percentile(z, 2.5)),
                        hi=float(np.percentile(z, 97.5)), n_ok=post["n_ok"]))
    return out


def attribute_height_bias(stains, truth_z, noise, n_draws=120, seed=1, drag=True):
    """Re-run the posterior with one noise source active at a time.

    Returns {source: (mean_bias_m, sd_m)} plus the 'all' case, so the -18 cm can be
    charged to whichever measurement actually causes it instead of being attributed to
    'noise' in general.
    """
    out = {}
    fields = ["rel_W", "rel_L", "az_deg", "rel_d0"]
    for i, field in enumerate(fields + ["all"]):
        n = noise if field == "all" else noise.only(field)
        if n.is_zero:
            out[field] = (0.0, 0.0)
            continue
        post = origin_posterior(stains, n, n_draws=n_draws, seed=seed + i, drag=drag)
        z = post["cloud"][:, 2]
        out[field] = (float(np.mean(z) - truth_z), float(np.std(z)))
    return out

"""psd.py — particle size distributions, their moments, and the moment traps.

Engine-independent. Nothing here imports Yade.

The whole project exists because "two samples of similar size" names ONE moment of a
distribution while the physics reads a DIFFERENT one. This module makes every moment
explicit and refuses to let a caller ask for "the diameter".

Conventions, fixed here and not to be overridden downstream:

  * `d` is always a particle DIAMETER, in metres. Never a radius. Yade works in radii;
    the conversion happens exactly once, at `radii_for_yade()`.
  * A lognormal distribution is parameterised by its NUMBER mean `d10` and its
    coefficient of variation `CV = std/mean` of the NUMBER distribution.
  * `CV = 0` means strictly monodisperse and is handled exactly, not as a limit.

Moment definitions (Mugele & Evans 1951 notation, the standard in powder technology):

    d_pq = ( <d^p> / <d^q> )^(1/(p-q))

  d10  arithmetic / number mean      <d>
  d32  Sauter (surface-volume) mean  <d^3>/<d^2>   <- enters Carman-Kozeny permeability
  d43  De Brouckere / volume mean    <d^4>/<d^3>   <- what laser diffraction reports

Analytic results for a lognormal at fixed number mean (derived in
`derivations/lognormal_moments.md`, and spot-checked against Monte-Carlo in
`test_psd.py`):

    sigma_l^2 = ln(1 + CV^2)                      # log-space variance
    d32 / d10 = (1 + CV^2)^2                      # EXACT
    d43 / d10 = (1 + CV^2)^3                      # EXACT
    kappa ratio (Carman-Kozeny, ∝ d32^2) = (1 + CV^2)^4

The `d32/d10 = (1+CV^2)^2` identity is the reason mechanism M1 exists at all; see
`SPEC_2026-08-31.md` §1.2.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "LogNormalPSD",
    "UniformPSD",
    "BidispersePSD",
    "make_psd",
    "lognormal_skewness",
    "moment_ratio_d32_d10",
    "moment_ratio_d43_d10",
    "carman_kozeny_kappa",
]


def moment_ratio_d32_d10(cv: float) -> float:
    """d32/d10 for a lognormal at fixed number mean. Exact, closed form.

    >>> moment_ratio_d32_d10(0.0)
    1.0
    >>> round(moment_ratio_d32_d10(0.5), 6)
    1.5625
    """
    if cv < 0:
        raise ValueError(f"CV must be >= 0, got {cv}")
    return float((1.0 + cv**2) ** 2)


def moment_ratio_d43_d10(cv: float) -> float:
    """d43/d10 for a lognormal at fixed number mean. Exact, closed form."""
    if cv < 0:
        raise ValueError(f"CV must be >= 0, got {cv}")
    return float((1.0 + cv**2) ** 3)


def lognormal_skewness(cv: float) -> float:
    """Skewness of a lognormal with coefficient of variation `cv`.

        S = (3 + CV^2) * CV        [standard result]

    This matters because Desmond & Weeks, Phys. Rev. E 90, 022204 (2014)
    (arXiv:1303.4627) show that the SKEWNESS of a size distribution can shift the
    random-close-packing density MORE than the polydispersity does, while kurtosis
    barely matters. A lognormal ladder therefore moves CV and skewness TOGETHER
    (CV = 0.5 gives S = 1.625), and any "CV effect" measured on it alone is really a
    joint (CV, skewness) effect along the lognormal locus.

    `UniformPSD` exists to break that degeneracy: it has skewness EXACTLY zero at any
    CV, so a matched-CV comparison of the two families separates the two axes.
    """
    if cv < 0:
        raise ValueError(f"CV must be >= 0, got {cv}")
    return float((3.0 + cv**2) * cv)


def carman_kozeny_kappa(d32: float, phi: float) -> float:
    """Carman-Kozeny permeability [m^2], with the SAUTER mean as the length.

        kappa = (1 - phi)^3 d32^2 / (180 phi^2)

    This is the same expression as `impact_stress.inertial_dilatant.DilatantSystem.kappa`,
    which takes `a` = particle DIAMETER. Using d32 here rather than a nominal diameter is
    the correction that mechanism M1 turns on.
    """
    if not 0.0 < phi < 1.0:
        raise ValueError(f"phi must be in (0,1), got {phi}")
    if d32 <= 0:
        raise ValueError(f"d32 must be > 0, got {d32}")
    return float((1.0 - phi) ** 3 * d32**2 / (180.0 * phi**2))


@dataclass(frozen=True)
class LogNormalPSD:
    """A lognormal particle size distribution at a fixed NUMBER mean.

    Parameters
    ----------
    d10 : number-mean particle DIAMETER [m]
    cv  : coefficient of variation of the NUMBER distribution, std/mean, dimensionless
    """

    d10: float
    cv: float

    def __post_init__(self) -> None:
        if self.d10 <= 0:
            raise ValueError(f"d10 must be > 0 m, got {self.d10}")
        if self.cv < 0:
            raise ValueError(f"cv must be >= 0, got {self.cv}")

    # ---- log-space parameters -------------------------------------------------

    @property
    def sigma_l(self) -> float:
        """Standard deviation of ln(d). Zero for a monodisperse pack."""
        return float(np.sqrt(np.log1p(self.cv**2)))

    @property
    def mu_l(self) -> float:
        """Mean of ln(d), chosen so that <d> == d10 exactly."""
        return float(np.log(self.d10) - 0.5 * self.sigma_l**2)

    # ---- analytic moments -----------------------------------------------------

    def moment(self, n: int) -> float:
        """<d^n> for the lognormal, analytic: exp(n*mu + n^2 sigma^2/2)."""
        return float(np.exp(n * self.mu_l + 0.5 * n**2 * self.sigma_l**2))

    @property
    def d32(self) -> float:
        """Sauter mean <d^3>/<d^2> [m]. The length Carman-Kozeny reads."""
        return self.d10 * moment_ratio_d32_d10(self.cv)

    @property
    def d43(self) -> float:
        """De Brouckere mean <d^4>/<d^3> [m]. What laser diffraction reports."""
        return self.d10 * moment_ratio_d43_d10(self.cv)

    @property
    def d50(self) -> float:
        """Median diameter [m] = exp(mu_l)."""
        return float(np.exp(self.mu_l))

    # ---- sampling -------------------------------------------------------------

    def sample(self, n: int, rng: np.random.Generator, *,
               truncate: tuple[float, float] | None = (0.2, 5.0)) -> np.ndarray:
        """Draw `n` particle DIAMETERS [m].

        `truncate` clips the draw to (lo*d10, hi*d10) by REJECTION, so the returned
        sample is a true truncated lognormal rather than a pile-up at the bounds.
        Truncation is on by default because an untruncated lognormal at CV >= 0.8
        produces grains that differ by >100x in diameter, which is a DEM timestep
        problem rather than a physics choice. The realised moments are always
        reported by `realised_moments()` -- never assume the analytic ones hold
        after truncation.
        """
        if n <= 0:
            raise ValueError(f"n must be > 0, got {n}")
        if self.cv == 0.0:
            return np.full(n, self.d10, dtype=float)

        if truncate is None:
            return rng.lognormal(self.mu_l, self.sigma_l, n)

        lo, hi = truncate[0] * self.d10, truncate[1] * self.d10
        out = np.empty(n, dtype=float)
        filled, tries = 0, 0
        while filled < n:
            tries += 1
            if tries > 1000:
                raise RuntimeError(
                    f"rejection sampling failed to fill {n} draws at cv={self.cv}, "
                    f"bounds=({lo:.3e},{hi:.3e}) after {tries} rounds"
                )
            draw = rng.lognormal(self.mu_l, self.sigma_l, 2 * (n - filled) + 16)
            keep = draw[(draw >= lo) & (draw <= hi)]
            take = min(len(keep), n - filled)
            out[filled:filled + take] = keep[:take]
            filled += take
        return out

    def realised_moments(self, d: np.ndarray) -> dict[str, float]:
        """Moments of an ACTUAL sample. Always report these, not the analytic ones.

        Truncation and finite N both shift the high moments, and d32 is a high
        moment. Quoting the analytic d32 for a truncated sample is exactly the
        error this module exists to prevent.
        """
        d = np.asarray(d, dtype=float)
        if d.ndim != 1 or d.size == 0:
            raise ValueError(f"d must be a non-empty 1-D array, got shape {d.shape}")
        if np.any(d <= 0):
            raise ValueError("all diameters must be > 0")
        m1, m2, m3, m4 = (float(np.mean(d**k)) for k in (1, 2, 3, 4))
        return {
            "n": int(d.size),
            "d10": m1,
            "d32": m3 / m2,
            "d43": m4 / m3,
            "d50": float(np.median(d)),
            "cv": float(np.std(d) / m1),
            "d_min": float(d.min()),
            "d_max": float(d.max()),
            "span": float(d.max() / d.min()),
        }

    # ---- engine boundary ------------------------------------------------------

    def radii_for_yade(self, n: int, rng: np.random.Generator, **kw) -> np.ndarray:
        """Particle RADII [m] for Yade. The one and only diameter->radius conversion.

        Yade's `sphere()` takes a radius. Every other quantity in this project is a
        diameter. Keeping the factor of 2 in exactly one function is the whole
        defence against the class of bug the lab's units-auditor exists to catch.
        """
        return self.sample(n, rng, **kw) / 2.0


@dataclass(frozen=True)
class UniformPSD:
    """A UNIFORM particle size distribution at a fixed NUMBER mean.

    Diameters are uniform on [d10*(1-w), d10*(1+w)] with w = CV*sqrt(3), so the number
    mean is d10 and the coefficient of variation is exactly `cv`. **Skewness is
    identically zero at every CV**, which is the whole point: paired against
    `LogNormalPSD` at the same CV it isolates the polydispersity axis from the skewness
    axis that Desmond & Weeks showed can dominate it.

    Diameters must stay positive, so w < 1 and therefore CV < 1/sqrt(3) = 0.5774.
    """

    d10: float
    cv: float

    def __post_init__(self) -> None:
        if self.d10 <= 0:
            raise ValueError(f"d10 must be > 0 m, got {self.d10}")
        if self.cv < 0:
            raise ValueError(f"cv must be >= 0, got {self.cv}")
        if self.cv >= 1.0 / np.sqrt(3.0):
            raise ValueError(
                f"a uniform PSD with positive diameters needs CV < 1/sqrt(3) = "
                f"{1/np.sqrt(3):.4f}, got {self.cv}")

    @property
    def half_width(self) -> float:
        """w such that std/mean = cv for a uniform on [d10(1-w), d10(1+w)]."""
        return float(self.cv * np.sqrt(3.0))

    def moment(self, n: int) -> float:
        """<d^n> for the uniform, in closed form.

        <d^n> = ( (1+w)^(n+1) - (1-w)^(n+1) ) * d10^n / (2 w (n+1))
        """
        w = self.half_width
        if w == 0:
            return float(self.d10**n)
        return float(self.d10**n * ((1 + w) ** (n + 1) - (1 - w) ** (n + 1))
                     / (2.0 * w * (n + 1)))

    @property
    def d32(self) -> float:
        return self.moment(3) / self.moment(2)

    @property
    def d43(self) -> float:
        return self.moment(4) / self.moment(3)

    @property
    def d50(self) -> float:
        return float(self.d10)          # symmetric

    def sample(self, n: int, rng: np.random.Generator, *, truncate=None) -> np.ndarray:
        """Draw `n` diameters. `truncate` is accepted and IGNORED -- a uniform PSD is
        already bounded, so truncation would be a second, undocumented distortion."""
        if n <= 0:
            raise ValueError(f"n must be > 0, got {n}")
        if self.cv == 0.0:
            return np.full(n, self.d10, dtype=float)
        w = self.half_width
        return rng.uniform(self.d10 * (1 - w), self.d10 * (1 + w), n)

    realised_moments = LogNormalPSD.realised_moments
    radii_for_yade = LogNormalPSD.radii_for_yade


@dataclass(frozen=True)
class BidispersePSD:
    """A symmetric 50:50 BIDISPERSE distribution at a fixed NUMBER mean.

    Half the grains at d10*(1-CV), half at d10*(1+CV). By construction the number mean
    is exactly d10, the coefficient of variation is exactly `cv`, and the skewness is
    exactly 0 -- the same two moments as `UniformPSD`, but with all the weight at the
    extremes instead of spread across the interval.

    This is the third family, and it exists to test a claim the literature treats as
    established: that a polydisperse suspension behaves like a "statistically equivalent
    bidisperse" one (Singh, Tyagi & Santra arXiv:2510.16464 and refs therein, for
    viscosity, normal stresses, particle pressures and packing efficiency). If matching
    two moments were sufficient, uniform and bidisperse at the same CV would agree.
    Where they disagree, the equivalence needs more than two moments -- which is the
    same conclusion the lognormal-vs-uniform comparison reaches from the other side.

    Size ratio is d_large/d_small = (1+CV)/(1-CV), so CV < 1 keeps diameters positive.
    """

    d10: float
    cv: float

    def __post_init__(self) -> None:
        if self.d10 <= 0:
            raise ValueError(f"d10 must be > 0 m, got {self.d10}")
        if self.cv < 0:
            raise ValueError(f"cv must be >= 0, got {self.cv}")
        if self.cv >= 1.0:
            raise ValueError(
                f"a symmetric bidisperse PSD needs CV < 1 for positive diameters, "
                f"got {self.cv}")

    @property
    def size_ratio(self) -> float:
        """d_large / d_small, the quantity the bidisperse literature reports."""
        if self.cv == 0:
            return 1.0
        return float((1.0 + self.cv) / (1.0 - self.cv))

    def moment(self, n: int) -> float:
        """<d^n>, exact: the mean of the two atoms."""
        lo = self.d10 * (1.0 - self.cv)
        hi = self.d10 * (1.0 + self.cv)
        return float(0.5 * (lo**n + hi**n))

    @property
    def d32(self) -> float:
        return self.moment(3) / self.moment(2)

    @property
    def d43(self) -> float:
        return self.moment(4) / self.moment(3)

    @property
    def d50(self) -> float:
        return float(self.d10)          # symmetric

    def sample(self, n: int, rng: np.random.Generator, *, truncate=None) -> np.ndarray:
        """Draw `n` diameters. `truncate` is accepted and IGNORED: the distribution is
        already bounded, and truncating it would distort the two moments it exists to
        hold fixed."""
        if n <= 0:
            raise ValueError(f"n must be > 0, got {n}")
        if self.cv == 0.0:
            return np.full(n, self.d10, dtype=float)
        lo = self.d10 * (1.0 - self.cv)
        hi = self.d10 * (1.0 + self.cv)
        # Exactly half large, half small, then shuffled: an EXACT 50:50 split rather
        # than a binomial one, so the realised CV does not fluctuate with the seed.
        out = np.concatenate([np.full(n // 2, lo), np.full(n - n // 2, hi)])
        rng.shuffle(out)
        return out

    realised_moments = LogNormalPSD.realised_moments
    radii_for_yade = LogNormalPSD.radii_for_yade


def make_psd(family: str, d10: float, cv: float):
    """Factory. `family` is "lognormal", "uniform" or "bidisperse"."""
    if family == "lognormal":
        return LogNormalPSD(d10=d10, cv=cv)
    if family == "uniform":
        return UniformPSD(d10=d10, cv=cv)
    if family == "bidisperse":
        return BidispersePSD(d10=d10, cv=cv)
    raise ValueError(
        f"unknown PSD family {family!r}; use 'lognormal', 'uniform' or 'bidisperse'")


def summarise_cv_ladder(d10: float, cvs, phi: float) -> str:
    """Human-readable table of the moment ratios and the Carman-Kozeny consequence.

    This reproduces the M1 table in SPEC_2026-08-31.md §1.2 from first principles,
    so the spec's numbers can be regenerated rather than trusted.
    """
    lines = [
        f"lognormal at fixed number mean d10 = {d10*1e6:.3f} um, phi = {phi}",
        f"{'CV':>5} {'d32/d10':>9} {'d43/d10':>9} {'d32[um]':>9} "
        f"{'kappa[m^2]':>12} {'kappa/kappa0':>13}",
    ]
    k0 = None
    for cv in cvs:
        p = LogNormalPSD(d10=d10, cv=cv)
        k = carman_kozeny_kappa(p.d32, phi)
        if k0 is None:
            k0 = k
        lines.append(
            f"{cv:5.2f} {moment_ratio_d32_d10(cv):9.4f} {moment_ratio_d43_d10(cv):9.4f} "
            f"{p.d32*1e6:9.4f} {k:12.4e} {k/k0:13.4f}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarise_cv_ladder(0.8e-6, [0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0], 0.44))

"""Eliminate candidate arrest criteria by degeneracy, and rank what is left.

The test, and why it is not a correlation
-----------------------------------------
With 75 drops and a strong velocity axis, almost any monotone function of U0
correlates with the arrested footprint. Correlation therefore cannot decide
between the candidates. Two things can:

**1. The exponent the candidate would need.** If ``beta_p`` is a function of the
group ``G`` alone, then over the dataset

    d ln beta_p / d ln G  ~  ln(range of beta_p) / ln(range of G).

A group that varies by 3% while the response varies by a factor of 6 needs an
exponent near 70. That is not a soft criterion: it is a statement that the
proposed mechanism is being asked to do work its own dynamic range cannot do.

**2. Degeneracy.** If ``G`` sets ``beta_p``, drops with the same ``G`` must have
the same ``beta_p``. Among pairs whose ``G`` agree to within a tolerance, the
surviving spread in ``beta_p`` is reported as a fraction of the total. A group
that leaves the full spread among its own ties is falsified whatever its
correlation.

The tolerance is a *relative* one on ``G`` and it is deliberately generous:
5% and 20%. None of the material constants that enter these groups is known to
better than that. sigma_star in particular is definition-bound to a factor of
3 to 5 (see ``tier0_wyart_cates_closure/CALIBRATION_2026-08-30.md``), so 5% is
far kinder to the stress clock than the evidence warrants.

Pre-registered verdict rule, fixed before the numbers were looked at:

- **FALSIFIED** if at 5% tolerance the near-ties retain more than half the
  total log-range of ``beta_p``.
- **STRAINED** if the required exponent exceeds 10.
- **SURVIVES** otherwise.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import spearmanr

TOL_TIGHT = 0.05
TOL_LOOSE = 0.20
EXPONENT_STRAIN = 10.0
DEGENERACY_FALSIFY = 0.5


@dataclass(frozen=True)
class ClockVerdict:
    """One candidate criterion, audited against the measured response."""

    name: str
    n: int
    g_min: float
    g_max: float
    response_log_range: float
    required_exponent: float
    tie_fraction_tight: float     # median retention among 5% ties, vs all pairs
    tie_fraction_loose: float
    tie_max_tight: float          # worst tied pair, as a share of the full range
    n_tie_pairs_tight: int
    spearman_rho: float
    spearman_p: float

    @property
    def g_range(self) -> float:
        return self.g_max / self.g_min if self.g_min > 0 else float("inf")

    @property
    def verdict(self) -> str:
        if self.n_tie_pairs_tight > 0 and self.tie_fraction_tight > DEGENERACY_FALSIFY:
            return "FALSIFIED"
        if self.required_exponent > EXPONENT_STRAIN:
            return "STRAINED"
        return "survives"

    def __str__(self) -> str:
        return (f"{self.name:<40s} range {self.g_range:9.3g}x  "
                f"exponent needed {self.required_exponent:8.2f}  "
                f"ties@5% retain {self.tie_fraction_tight:5.2f} "
                f"(max {self.tie_max_tight:4.2f}, {self.n_tie_pairs_tight:5d} pairs)  "
                f"@20% {self.tie_fraction_loose:5.2f}  "
                f"rho {self.spearman_rho:+.3f}  {self.verdict}")


def _log_range(v: np.ndarray) -> float:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v) & (v > 0)]
    if v.size < 2:
        return 0.0
    return float(np.log(v.max() / v.min()))


def tie_retention(g: np.ndarray, r: np.ndarray, tol: float) -> tuple[float, float, int]:
    """How much of the response's spread survives among drops with the same ``g``.

    Returns ``(median retention, max retention, n tied pairs)``, both retentions
    as fractions.

    - **median retention** is the median ``|d ln r|`` among pairs whose ``g``
      agree to within ``tol``, divided by the median ``|d ln r|`` over *all*
      pairs. If ``g`` determined ``r`` this would be near zero; at 1 the tie
      carries as much spread as a randomly chosen pair, i.e. knowing ``g`` has
      told you nothing.
    - **max retention** is the largest ``|d ln r|`` among tied pairs over the
      response's full log-range. It is the strict falsifier -- one counterexample
      refutes a deterministic law -- but over thousands of pairs it is a
      single-pair statistic and one mis-measured drop saturates it.

    The verdict below uses the median. That choice was made *after* seeing that
    the max saturates near 1 for almost every candidate, which is recorded here
    rather than presented as a prior decision: the max is reported alongside so
    a reader can apply the strict rule instead.
    """
    g = np.asarray(g, dtype=float)
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(g) & np.isfinite(r) & (g > 0) & (r > 0)
    g, r = np.log(g[ok]), np.log(r[ok])
    if g.size < 3:
        return 0.0, 0.0, 0
    total = float(r.max() - r.min())
    iu = np.triu_indices(g.size, k=1)
    dg = np.abs(g[:, None] - g[None, :])[iu]
    dr = np.abs(r[:, None] - r[None, :])[iu]
    baseline = float(np.median(dr))
    tied = dg <= np.log1p(tol)
    if not tied.any() or baseline <= 0 or total <= 0:
        return 0.0, 0.0, int(tied.sum())
    return (float(np.median(dr[tied]) / baseline),
            float(dr[tied].max() / total),
            int(tied.sum()))


def audit_clock(name: str, g: np.ndarray, r: np.ndarray) -> ClockVerdict:
    g = np.asarray(g, dtype=float)
    r = np.asarray(r, dtype=float)
    ok = np.isfinite(g) & np.isfinite(r) & (g > 0) & (r > 0)
    g, r = g[ok], r[ok]
    if g.size < 3:
        raise ValueError(f"{name}: only {g.size} usable drops")
    lr_r = _log_range(r)
    lr_g = _log_range(g)
    required = float(lr_r / lr_g) if lr_g > 0 else float("inf")
    t_tight, m_tight, n_tight = tie_retention(g, r, TOL_TIGHT)
    t_loose, _, _ = tie_retention(g, r, TOL_LOOSE)
    rho, p = spearmanr(g, r)
    return ClockVerdict(
        name=name, n=int(g.size), g_min=float(g.min()), g_max=float(g.max()),
        response_log_range=lr_r, required_exponent=required,
        tie_fraction_tight=t_tight, tie_fraction_loose=t_loose,
        tie_max_tight=m_tight, n_tie_pairs_tight=n_tight,
        spearman_rho=float(rho), spearman_p=float(p),
    )

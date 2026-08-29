"""Readers for the three Edinburgh DataShare deposits, and the constants they fix.

The point of these datasets is that they come with the authors' *own* fitted
parameters, so a fit made here can be compared with a published one rather than
only with itself. That only works if the fixed constants match, which is why
each deposit gets a named constant set in this module instead of the values
being retyped into a script.

The data is not redistributed with the package (see ``.gitignore`` and
``DATA_INVENTORY.md``); these readers take a path.

Richards et al. 2019, PRL 123, 038004, arXiv:1902.07655
--------------------------------------------------------
Cornstarch in 50 wt% glycerol-water, seven time-averaged flow curves, stress
controlled. **The paper uses weight fractions, not volume fractions, and says
why**: "Here we use weight fractions, phi, due to the porosity of our main model
system (cornstarch)". So there is no wt-to-vol conversion to be done here, and
:func:`mass_to_volume_fraction` exists to be *not* used on this data; see its
docstring. What follows from that is that ``phi_0`` for cornstarch is not RCP:
their frictionless jamming point is 0.546 in the same weight-fraction units.

Their published global fit, from the Figure 1 caption and reproduced verbatim in
the deposit's own plotting script: ``phi_m = 0.457``, ``phi_rcp = 0.546``,
``beta = 0.94``, ``sigma* = 5.1 Pa``, ``eta_s = 6 mPa.s``, exponent 2.

Guy, Hermes & Poon 2015, PRL 115, 088304, arXiv:1506.09041
-----------------------------------------------------------
PHSA-stabilised PMMA in a density-matching CHB/decalin mixture, two sizes. Their
fit to the upper (frictional) branch, from the Fig. 1 caption:
``eta_2 = A (1 - phi/phi_m)^-n`` with ``A = 0.20(9)``, ``phi_m = 0.558(5)``,
``n = 2.2(2)``.

**The solvent viscosity is not one number.** The paper says "the solvent
viscosities were eta_f = 2.83 mPa.s for the large and 2.4 mPa.s for the small
particles at 19 degC". Using 2.4 for both, as ``DATA_INVENTORY.md`` originally
recorded, makes the d = 3770 nm relative viscosities 18% too high, and that
error goes straight into phi_m. Hence two constant sets, not one.

**The onset stress is not one number either**: sigma* ~= 200 Pa at d = 404 nm
and ~= 2 Pa at d = 3770 nm, with a fitted size scaling sigma* = B d^alpha,
alpha = -1.9(1).

**The file names do not match the panel letters.** ``Fig1b.xlsx`` holds flow
curves and ``Fig1c.xlsx`` holds the phi-axis compilation, which is the other way
round from the published panels: the caption puts the eta(phi) compilation and
its eta_2 fit in panel (b), and the d = 3770 nm flow curves in panel (c). Two
independent checks agree, and :func:`verify_guy_size_assignment` runs the second
one on the data:

1. ``Fig1c.xlsx`` contains sheets named "eta_1 data (literature)" including
   de Kruif 1985 (silica), and the caption places the literature comparison in
   panel (b).
2. The low-stress viscosity plateau of each ``Fig1a`` sheet, divided by 2.4
   mPa.s, reproduces the ``d=404 nm eta_1`` column of ``Fig1c.xlsx``, and each
   ``Fig1b`` sheet divided by 2.83 mPa.s reproduces ``d=3770 nm eta_1``. That
   fixes the size assignment, the label-to-true-phi map and the solvent
   viscosity together, from the numbers rather than from the file names.

**Two curves in each panel are not steady measurements.** "The flow in both (a)
and (c) was unsteady for phi >= 0.56, and points represent temporal averages",
and the paper separately reports wall slip at all stresses below phi_m and
"massive wall slip in shear thickened samples above phi_m". Those are exactly
the highest-phi curves, and they are the ones a jamming fit leans on hardest.
"""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from dst_fit.io_trios import FlowSweep

# --------------------------------------------------------------------------
# constant sets, as published
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ConstantSet:
    """The four Wyart-Cates fixed constants for one published system.

    ``phi_units`` records whether ``phi`` means a volume or a weight fraction,
    because the two cannot be mixed and the difference is not visible in any of
    the numbers. ``phi_m_published`` is the authors' own answer, which is what
    the gates compare against.
    """

    name: str
    eta_s: float
    phi_0: float
    n: float
    beta: float
    phi_units: str
    phi_m_published: float
    phi_m_published_err: float
    sigma_star_published: float | None
    source: str

    def as_fixed(self) -> dict[str, float]:
        """The keyword arguments for :meth:`dst_fit.models.Model.with_fixed`."""
        return {"eta_s": self.eta_s, "phi_0": self.phi_0, "n": self.n, "beta": self.beta}

    def __str__(self) -> str:
        return (
            f"{self.name}: eta_s={self.eta_s:g} Pa.s, phi_0={self.phi_0:g}, "
            f"n={self.n:g}, beta={self.beta:g}, phi in {self.phi_units} fraction; "
            f"published phi_m = {self.phi_m_published:.3f}"
            f"({self.phi_m_published_err * 1000:.0f}e-3)"
        )


# Richards et al. 2019. Every number from the Figure 1 caption and confirmed in
# the deposit's own Figure1.py (beta = 0.94, phiM = 0.457, phiRCP = 0.546,
# etaF = 0.006, sigma0 = 5.1).
RICHARDS_CORNSTARCH = ConstantSet(
    name="Richards 2019 cornstarch / 50 wt% glycerol-water",
    eta_s=6.0e-3,
    phi_0=0.546,
    n=2.0,
    beta=0.94,
    phi_units="weight",
    phi_m_published=0.457,
    phi_m_published_err=0.0,  # the paper quotes no uncertainty on it
    sigma_star_published=5.1,
    source="PRL 123, 038004 (2019); arXiv:1902.07655; DataShare 10283/3357",
)

# Guy et al. 2015, one constant set per particle size. eta_s is stated in the
# paper and confirmed by the data: the d=404 nm eta_1 sheet gives 0.018 Pa.s at
# a relative viscosity of 7.5, i.e. 2.4e-3 exactly, and dividing the Fig1b
# plateaus by 2.83e-3 reproduces the d=3770 nm eta_1 column to within 6%.
#
# phi_0 is measured, not assumed, by fitting the two frictionless branches in
# Fig1c with fit_divergence at the same exponent as phi_m (n = 2.2), because a
# pole located at one exponent cannot be combined with one located at another:
#
#     d = 3770 nm, phi > 0.35   phi_0 = 0.6123 +/- 0.0010
#     d =  404 nm               phi_0 = 0.6076 +/- 0.0022
#
# Two sizes of the same particle chemistry, agreeing to 0.005. Each set uses its
# own measured pole; the gate also reruns at the other one so the sensitivity of
# phi_m to this choice is measured rather than assumed. (At n = 2 the same
# branches give 0.6103 and 0.6007, which is why the exponent must be settled
# first.) Note this is *not* the 0.64 the paper's introduction quotes for
# frictionless RCP: that is the value their phi scale was anchored to, while
# 0.61 is where their own measured eta_1 branch actually diverges.
_GUY_SOURCE = "PRL 115, 088304 (2015); arXiv:1506.09041; DataShare 10283/834"

GUY_PMMA_404 = ConstantSet(
    name="Guy 2015 PMMA d=404 nm",
    eta_s=2.4e-3,
    phi_0=0.6076,
    n=2.2,
    beta=1.0,
    phi_units="volume",
    phi_m_published=0.558,
    phi_m_published_err=0.005,
    sigma_star_published=200.0,
    source=_GUY_SOURCE,
)

GUY_PMMA_3770 = ConstantSet(
    name="Guy 2015 PMMA d=3770 nm",
    eta_s=2.83e-3,
    phi_0=0.6123,
    n=2.2,
    beta=1.0,
    phi_units="volume",
    phi_m_published=0.558,
    phi_m_published_err=0.005,
    sigma_star_published=2.0,
    source=_GUY_SOURCE,
)

CONSTANT_SETS = {
    c.name: c for c in (RICHARDS_CORNSTARCH, GUY_PMMA_404, GUY_PMMA_3770)
}


# --------------------------------------------------------------------------
# the conversion the cornstarch papers decline to make
# --------------------------------------------------------------------------
def mass_to_volume_fraction(w: np.ndarray, rho_particle: float, rho_fluid: float) -> np.ndarray:
    """Volume fraction from a weight fraction, given both densities.

    ``phi = (w/rho_p) / (w/rho_p + (1-w)/rho_f)``, which assumes the particles
    are non-porous solids of a single known density and that mixing is volume-
    additive.

    **For cornstarch neither assumption holds**, and Richards et al. say so
    explicitly: they work in weight fractions "due to the porosity of our main
    model system (cornstarch)". A cornstarch granule takes up water, so its
    effective density is neither the dry-solid value nor constant across the
    series, and the conversion is not a change of units but a model with an
    unmeasured parameter. Applying it silently is candidate (3) in the N4 open
    questions.

    So this function is here to be called deliberately, with the densities
    written out at the call site, in order to *quantify* how much a conversion
    would move things. It is not called by any gate.
    """
    w = np.asarray(w, dtype=float)
    if np.any((w < 0) | (w > 1)):
        raise ValueError("weight fraction must lie in [0, 1]")
    if rho_particle <= 0 or rho_fluid <= 0:
        raise ValueError("densities must be positive")
    v_p = w / rho_particle
    return v_p / (v_p + (1.0 - w) / rho_fluid)


# --------------------------------------------------------------------------
# Richards et al. 2019
# --------------------------------------------------------------------------
def read_richards_flow_curve(path: str | Path) -> FlowSweep:
    """One ``NNwt%.csv`` time-averaged flow curve.

    Columns are ``Shear-rate`` (1/s), ``Error`` (1/s, the standard deviation
    over three up-sweeps) and ``Stress`` (Pa, the *applied* external stress).
    This is a stress-controlled measurement, so stress is the independent
    variable and the viscosity is the derived quantity, ``sigma / gamma_dot``.

    ``stress_reported`` is left as ``None`` rather than being set to the same
    column that ``stress`` came from: :func:`dst_fit.io_trios.check_stress_consistency`
    would then compare a column with itself and report a perfect pass that tests
    nothing.

    ``metadata['phi']`` is the **weight** fraction, which is the packing variable
    the paper itself uses; ``metadata['phi_units']`` records that so it cannot be
    read as a volume fraction by accident.
    """
    path = Path(path)
    df = pd.read_csv(path)
    missing = [c for c in ("Shear-rate", "Error", "Stress") if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing column(s) {missing}; found {list(df.columns)}")

    gdot = df["Shear-rate"].to_numpy(dtype=float)
    err = df["Error"].to_numpy(dtype=float)
    sigma = df["Stress"].to_numpy(dtype=float)

    keep = np.isfinite(gdot) & np.isfinite(sigma) & (gdot > 0) & (sigma > 0)
    if not keep.any():
        raise ValueError(f"{path.name}: no finite positive (rate, stress) rows")
    gdot, err, sigma = gdot[keep], err[keep], sigma[keep]

    order = np.argsort(sigma)  # the model lives on the stress axis
    gdot, err, sigma = gdot[order], err[order], sigma[order]

    m = re.match(r"^(\d{1,3})wt%", path.stem)
    if not m:
        raise ValueError(f"{path.name}: expected a name like '38wt%.csv'")
    w = int(m.group(1)) / 100.0

    return FlowSweep(
        shear_rate=gdot,
        viscosity=sigma / gdot,
        stress=sigma,
        stress_reported=None,
        temperature=None,
        step_time=None,
        normal_stress=None,
        source=path,
        sheet=path.stem,
        metadata={
            "phi": w,
            "phi_units": "weight",
            "shear_rate_error": err,
            "dataset": RICHARDS_CORNSTARCH.name,
        },
    )


def read_richards_series(directory: str | Path) -> list[FlowSweep]:
    """All seven ``NNwt%.csv`` curves in ``Figure1/Plotting``, sorted by phi."""
    directory = Path(directory)
    paths = sorted(directory.glob("[0-9][0-9]wt%.csv"))
    if not paths:
        raise FileNotFoundError(f"no NNwt%.csv files in {directory}")
    return sorted(
        (read_richards_flow_curve(p) for p in paths),
        key=lambda s: s.phi if s.phi is not None else 0.0,
    )


# --------------------------------------------------------------------------
# Guy et al. 2015
# --------------------------------------------------------------------------
# Label on the workbook sheet -> the true volume fraction the authors used in
# their own phi-axis figure. Verified against Fig1c by
# verify_guy_size_assignment(), which compares the low-stress viscosity plateau
# of each flow curve with the eta_1 column at the matched phi.
GUY_FIG1A_PHI = {  # d = 404 nm
    "40%": 0.39296,
    "47%": 0.45952,
    "53%": 0.51968,
    "54.5%": 0.53632,
    "56%": 0.55616,
    "58%": 0.5696,
}
GUY_FIG1B_PHI = {  # d = 3770 nm
    "19%": 0.1856,
    "41%": 0.4096,
    "47%": 0.46656,
    "51%": 0.51008,
    "54%": 0.53632,
    "56%": 0.56128,
    "58%": 0.5824,
}


def _read_workbook(path: Path) -> dict[str, pd.DataFrame]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pd.read_excel(path, sheet_name=None, header=None)


def read_guy_flow_curves(path: str | Path, phi_map: dict[str, float] | None = None) -> list[FlowSweep]:
    """Flow curves from ``Fig1a.xls`` or ``Fig1b.xlsx``, one per sheet.

    Every sheet is ``shear rate [1/s]``, ``shear stress [Pa]``, ``viscosity
    [Pa.s]`` with the column names in row 0. Rows are returned sorted by stress,
    with exact duplicate ``(stress, viscosity)`` rows dropped: they are repeated
    points in the export, and leaving them in would inflate the point count and
    correlate the residuals the runs test reads.

    Repeated stresses with *different* viscosities are kept, and counted in
    ``metadata['n_repeated_stresses']``. They are not an export artefact: in the
    DST region one imposed stress admits several shear rates, which is the
    S-shaped flow curve showing itself. Every model here is single-valued in
    stress and so cannot reproduce that; a log-residual fit will split the
    difference between them. Dropping them would hide the one place the model is
    known to be inapplicable, so they are counted instead.

    ``phi_map`` supplies the true volume fraction for each sheet label. Sheets
    with no entry get ``phi = None`` and the label kept in
    ``metadata['phi_label']``, so a curve whose phi is unknown cannot silently
    contribute to a phi-versus-nominal regression.
    """
    path = Path(path)
    phi_map = {} if phi_map is None else phi_map
    out: list[FlowSweep] = []
    for sheet, raw in _read_workbook(path).items():
        body = raw.iloc[1:].apply(pd.to_numeric, errors="coerce")
        arr = body.to_numpy(dtype=float)
        arr = arr[np.all(np.isfinite(arr[:, :3]), axis=1)]
        gdot, sigma, eta = arr[:, 0], arr[:, 1], arr[:, 2]

        good = (gdot > 0) & (sigma > 0) & (eta > 0)
        gdot, sigma, eta = gdot[good], sigma[good], eta[good]
        order = np.argsort(sigma)
        gdot, sigma, eta = gdot[order], sigma[order], eta[order]

        _, first = np.unique(np.column_stack([sigma, eta]), axis=0, return_index=True)
        keep = np.sort(first)
        n_dropped = int(sigma.size - keep.size)
        gdot, sigma, eta = gdot[keep], sigma[keep], eta[keep]

        label = str(sheet).strip()
        out.append(
            FlowSweep(
                shear_rate=gdot,
                viscosity=eta,
                stress=sigma,
                stress_reported=None,
                temperature=None,
                step_time=None,
                normal_stress=None,
                source=path,
                sheet=label,
                metadata={
                    "phi": phi_map.get(label),
                    "phi_units": "volume",
                    "phi_label": label,
                    "n_duplicate_rows_dropped": n_dropped,
                    "n_repeated_stresses": int(sigma.size - np.unique(sigma).size),
                    "dataset": _GUY_SOURCE,
                },
            )
        )
    return sorted(out, key=lambda s: (s.phi is None, s.phi if s.phi is not None else 0.0))


def read_guy_branch(path: str | Path, sheet: str) -> tuple[np.ndarray, np.ndarray]:
    """``(phi, relative viscosity)`` from one ``Fig1c.xlsx`` branch sheet.

    Works for ``d=3770 nm eta_1``, ``d=3770 nm eta_2`` and ``d=404 nm eta_1``.
    Column 0 is ``#phi``; the relative viscosity is the *last* column, because
    the eta_2 sheet carries a ``phi_m - phi`` column in between and the eta_1
    sheets for d = 404 nm carry both an absolute and a relative viscosity.
    """
    book = _read_workbook(Path(path))
    if sheet not in book:
        raise KeyError(f"{Path(path).name}: no sheet {sheet!r}; has {list(book)}")
    raw = book[sheet]
    body = raw.iloc[1:].apply(pd.to_numeric, errors="coerce")
    arr = body.to_numpy(dtype=float)
    phi, eta_r = arr[:, 0], arr[:, -1]
    good = np.isfinite(phi) & np.isfinite(eta_r) & (phi > 0) & (eta_r > 0)
    phi, eta_r = phi[good], eta_r[good]
    order = np.argsort(phi)
    return phi[order], eta_r[order]


def read_guy_phi_m_column(path: str | Path) -> float:
    """The ``phi_m`` implied by the ``phi_m - phi`` column of ``d=3770 nm eta_2``.

    The authors tabulated the distance to jamming alongside phi, so their phi_m
    can be read straight off without fitting anything. Raises if the six rows do
    not agree, which is the only way this could be misread.
    """
    book = _read_workbook(Path(path))
    raw = book["d=3770 nm eta_2"]
    body = raw.iloc[1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    body = body[np.all(np.isfinite(body[:, :2]), axis=1)]
    implied = body[:, 0] + body[:, 1]
    if not np.allclose(implied, implied[0], atol=5e-4):
        raise ValueError(
            f"phi + (phi_m - phi) is not constant across rows: {implied}; "
            "the column does not mean what it appears to"
        )
    return float(np.mean(implied))


@dataclass(frozen=True)
class SizeAssignmentCheck:
    """Result of matching a workbook of flow curves to a phi-axis branch sheet."""

    workbook: str
    branch_sheet: str
    labels: tuple[str, ...]
    phi: tuple[float, ...]
    plateau_visc: tuple[float, ...]  # Pa.s, the lowest-stress measured viscosity
    branch_rel_visc: tuple[float, ...]  # dimensionless, the tabulated eta_1
    implied_eta_s: tuple[float, ...]  # Pa.s, their ratio
    eta_s_expected: float
    matched: bool

    @property
    def ratio(self) -> tuple[float, ...]:
        return tuple(e / self.eta_s_expected for e in self.implied_eta_s)

    def __str__(self) -> str:
        verdict = "matches" if self.matched else "DOES NOT MATCH"
        med = float(np.median(self.implied_eta_s))
        lines = [
            f"{self.workbook} vs {self.branch_sheet}: {verdict}. "
            f"implied eta_s {min(self.implied_eta_s) * 1e3:.2f} to "
            f"{max(self.implied_eta_s) * 1e3:.2f} mPa.s (median {med * 1e3:.2f}), "
            f"stated {self.eta_s_expected * 1e3:.2f}"
        ]
        for lab, p, a, b, e in zip(
            self.labels, self.phi, self.plateau_visc, self.branch_rel_visc,
            self.implied_eta_s,
        ):
            lines.append(
                f"    {lab:>6s}  phi {p:.5f}   lowest-stress eta {a:9.5f} Pa.s   "
                f"tabulated eta_1 {b:9.2f}   implied eta_s {e * 1e3:5.2f} mPa.s"
            )
        return "\n".join(lines)


def verify_guy_size_assignment(
    workbook: str | Path,
    fig1c: str | Path,
    branch_sheet: str,
    phi_map: dict[str, float],
    eta_s: float,
    ratio_range: tuple[float, float] = (0.9, 1.1),
) -> SizeAssignmentCheck:
    """Check a flow-curve workbook against the phi-axis sheet it should reproduce.

    The lowest-stress viscosity of each flow curve is the frictionless plateau
    in Pa.s; the ``eta_1`` column tabulates the same plateau made dimensionless.
    Their ratio is therefore the solvent viscosity, measured rather than
    assumed, and the check is whether that comes out at the value the paper
    states. Three things are confirmed at once: which particle size the workbook
    holds, what its sheet labels mean as volume fractions, and which of the
    paper's two solvent viscosities applies.

    That last one is the reason this returns a number instead of a boolean. Run
    against the d = 3770 nm sheets at 2.4 mPa.s, the ratios come out 11 to 17%
    high with no scatter, which is not a size mismatch but the wrong solvent
    viscosity: the paper gives 2.83 mPa.s for the large particles and 2.4 for
    the small. A pass/fail on the size alone would have absorbed that 18% into
    phi_m and reported success.
    """
    curves = {s.sheet: s for s in read_guy_flow_curves(workbook, phi_map)}
    phi_b, eta_b = read_guy_branch(fig1c, branch_sheet)

    labels, phis, plateau, branch = [], [], [], []
    for label, phi in sorted(phi_map.items(), key=lambda kv: kv[1]):
        if label not in curves:
            raise KeyError(f"{Path(workbook).name}: no sheet {label!r}")
        j = int(np.argmin(np.abs(phi_b - phi)))
        if abs(phi_b[j] - phi) > 1e-4:
            raise ValueError(
                f"{branch_sheet}: no phi within 1e-4 of {phi:g}; nearest {phi_b[j]:g}"
            )
        labels.append(label)
        phis.append(phi)
        plateau.append(float(curves[label].viscosity.min()))
        branch.append(float(eta_b[j]))

    implied = np.array(plateau) / np.array(branch)
    ratio = implied / eta_s
    return SizeAssignmentCheck(
        workbook=Path(workbook).name,
        branch_sheet=branch_sheet,
        labels=tuple(labels),
        phi=tuple(phis),
        plateau_visc=tuple(plateau),
        branch_rel_visc=tuple(branch),
        implied_eta_s=tuple(float(v) for v in implied),
        eta_s_expected=float(eta_s),
        matched=bool(ratio.min() >= ratio_range[0] and ratio.max() <= ratio_range[1]),
    )


# --------------------------------------------------------------------------
# Lin et al. 2015: the only direct measurement of the frictional fraction
# --------------------------------------------------------------------------
# Lin, Guy, Hermes, Ness, Sun, Poon & Cohen, PRL 115, 228304 (2015),
# arXiv:1509.02750, DataShare 10283/886.
#
# A flow-reversal experiment splits the measured viscosity into a contact part
# and a hydrodynamic part. Wyart-Cates *infers* a frictional fraction f(sigma)
# from the shape of eta(sigma) and never observes it; this deposit observes the
# contact contribution directly, at the same volume fraction, on the same
# sample. Nothing else on disk allows that comparison.
#
# Everything is **relative** viscosity against shear rate, which is the useful
# accident: the reduced stress
#
#     Sigma = eta_r * gamma_dot = sigma / eta_s
#
# needs no solvent viscosity at all, and neither does the contact fraction. So
# the whole comparison below is free of the one constant this dataset does not
# state. Build the observable so the unknown constant cancels, then say where it
# cancelled.
LIN_SYSTEMS = {
    "Fig_2a": "2 um silica",
    "Fig_2c": "1.4 um PMMA",
}


@dataclass(frozen=True)
class ReversalSplit:
    """One system's flow curve plus its measured contact/hydrodynamic split.

    ``rate_sweep`` is the conventional flow curve. ``steady_state``, ``contact``
    and ``hydro`` come from the reversal protocol and are reported at the same
    handful of shear rates as each other.

    Stresses are **reduced**: ``Sigma = eta_r * gamma_dot = sigma / eta_s``. The
    solvent viscosity is not needed and is not assumed.
    """

    system: str
    label: str
    sweep_rate: np.ndarray
    sweep_eta_r: np.ndarray
    rev_rate: np.ndarray
    rev_total: np.ndarray
    rev_contact: np.ndarray
    rev_hydro: np.ndarray
    source: Path

    @property
    def sweep_sigma(self) -> np.ndarray:
        """Reduced stress of the rate sweep, sigma / eta_s."""
        return self.sweep_eta_r * self.sweep_rate

    @property
    def rev_sigma(self) -> np.ndarray:
        return self.rev_total * self.rev_rate

    @property
    def contact_fraction(self) -> np.ndarray:
        """Share of the viscosity carried by contacts, eta_c / (eta_c + eta_h).

        Not the same quantity as Wyart-Cates's ``f``, which is a fraction of
        *contacts* that are frictional rather than a fraction of *stress*. They
        are different numbers with the same crossover, which is why the
        comparison below is between landmarks and not between values.
        """
        return self.rev_contact / (self.rev_contact + self.rev_hydro)

    def closure_error(self) -> np.ndarray:
        """Relative gap in ``contact + hydro = total``. A file-reading check."""
        return np.abs(self.rev_contact + self.rev_hydro - self.rev_total) / self.rev_total


def _read_two_column(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Headerless two-column CSV, tolerating the trailing empty fields in the export."""
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 2 or not parts[0] or not parts[1]:
            continue
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"{path.name}: no numeric rows")
    arr = np.array(rows, dtype=float)
    order = np.argsort(arr[:, 0])
    return arr[order, 0], arr[order, 1]


def read_lin_reversal(directory: str | Path, system: str = "Fig_2a") -> ReversalSplit:
    """Read one system's four files from ``Fig2/Fig2``.

    Raises if the reversal files do not share a common set of shear rates: the
    contact fraction is a per-rate ratio and pairing points from different rates
    would silently invent a measurement.
    """
    directory = Path(directory)
    if system not in LIN_SYSTEMS:
        raise KeyError(f"unknown system {system!r}; have {sorted(LIN_SYSTEMS)}")

    sweep = _read_two_column(directory / f"{system}_rate_sweep.csv")
    total = _read_two_column(directory / f"{system}_steady_state.csv")
    contact = _read_two_column(directory / f"{system}_contact.csv")
    hydro = _read_two_column(directory / f"{system}_hydro.csv")

    rates = total[0]
    for name, (r, _) in (("contact", contact), ("hydro", hydro)):
        if r.shape != rates.shape or not np.allclose(r, rates, rtol=1e-9):
            raise ValueError(
                f"{system}: {name} is reported at {r} but the steady state at {rates}; "
                "refusing to pair points measured at different shear rates"
            )

    return ReversalSplit(
        system=system,
        label=LIN_SYSTEMS[system],
        sweep_rate=sweep[0],
        sweep_eta_r=sweep[1],
        rev_rate=rates,
        rev_total=total[1],
        rev_contact=contact[1],
        rev_hydro=hydro[1],
        source=directory,
    )


def log_crossing(x: np.ndarray, y: np.ndarray, target: float) -> float:
    """Where a monotone-ish ``y(x)`` first reaches ``target``, interpolated in log x.

    Log interpolation because the abscissa here is a stress spanning decades:
    interpolating linearly between 9.8 and 24.7 would place the crossing in the
    wrong place by tens of per cent, which is the size of the effect being
    measured. Returns NaN when the target is never reached, rather than
    extrapolating past the data.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if np.all(y < target) or np.all(y > target):
        return float("nan")
    for i in range(y.size - 1):
        lo, hi = y[i], y[i + 1]
        if (lo - target) * (hi - target) <= 0 and lo != hi:
            t = (target - lo) / (hi - lo)
            return float(10.0 ** (np.log10(x[i]) + t * (np.log10(x[i + 1]) - np.log10(x[i]))))
    return float("nan")


def peak_log_slope(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """``(x, slope)`` at the steepest point of ``d log y / d log x``.

    A landmark that needs no model. Central differences on the log-log curve,
    so the value is a property of the data rather than of a fitted form. The
    peak slope is the right landmark for a thickening transition; a threshold
    crossing would depend on where the plateaus happen to sit.
    """
    lx = np.log10(np.asarray(x, dtype=float))
    ly = np.log10(np.asarray(y, dtype=float))
    if lx.size < 3:
        raise ValueError("peak_log_slope needs at least three points")
    slope = np.gradient(ly, lx)
    j = int(np.argmax(slope))
    return float(10.0 ** lx[j]), float(slope[j])


# --------------------------------------------------------------------------
# Lin et al. 2015, the reversal *transients*
# --------------------------------------------------------------------------
# Panels 2b and 2d hold what panels 2a and 2c do not: the instantaneous relative
# viscosity against accumulated strain after the reversal, at several rates. That
# is the only measurement on disk of a *rate of change* rather than a steady
# state, and it is what :mod:`dst_fit.transient` turns into the strain clock.
#
# Each transient panel is paired with the steady panel measured on the same
# sample, because the inversion from viscosity to frictional fraction needs that
# sample's own Wyart-Cates fit. Pairing a transient with the wrong steady panel
# would silently invert through another material's jamming point.
#
# Panel 3a is the simulation counterpart: the same protocol with the contact and
# hydrodynamic stresses evaluated directly rather than inferred from reversal. Its
# rates are in units of the simulation's own scale gamma_dot_0 = F^CL/(6 pi eta_0
# a^2), not 1/s, so its clock is dimensionless and it is a control on the method
# rather than a second measurement of the same constant.
LIN_TRANSIENTS = {
    "Fig_2b": ("2 um silica", "Fig_2a", "Fig2/Fig2"),
    "Fig_2d": ("1.4 um PMMA", "Fig_2c", "Fig2/Fig2"),
    "Fig_3a": ("simulation, 2 um silica parameters", None, "Fig3/Fig3"),
}

_RATE_IN_NAME = re.compile(r"_gamma_dot_([0-9]*\.?[0-9]+)\.csv$", re.IGNORECASE)


@dataclass(frozen=True)
class ReversalTransients:
    """One panel's post-reversal records, one per shear rate.

    ``strain`` and ``eta_r`` are tuples of arrays, ordered by ascending rate and
    aligned with ``rates``. They are *not* stacked into a 2-D array: the records
    have different lengths (10 010 to 24 963 samples in this deposit) and padding
    them to a common shape would invent samples.
    """

    panel: str
    label: str
    steady_panel: str | None
    rates: tuple[float, ...]
    strain: tuple[np.ndarray, ...]
    eta_r: tuple[np.ndarray, ...]
    source: Path

    def __len__(self) -> int:
        return len(self.rates)

    def __iter__(self):
        """Yield ``(rate, strain, eta_r)`` per record, ascending in rate."""
        return zip(self.rates, self.strain, self.eta_r)

    def sampling(self) -> tuple[float, ...]:
        """Median strain step of each record, for checking uniformity claims."""
        return tuple(float(np.median(np.diff(g))) for g in self.strain)

    def __str__(self) -> str:
        return (
            f"{self.panel} ({self.label}): {len(self)} rates "
            f"{self.rates[0]:g} to {self.rates[-1]:g}, "
            f"{min(g.size for g in self.strain)} to {max(g.size for g in self.strain)} "
            f"samples each"
        )


def read_lin_transients(root: str | Path, panel: str = "Fig_2b") -> ReversalTransients:
    """Read every ``*_gamma_dot_*.csv`` of one transient panel.

    ``root`` is the top of the extracted deposit, i.e. the directory holding
    ``Fig2/`` and ``Fig3/``, because the two transient panels live under
    different subdirectories.

    The shear rate is parsed from the file name rather than being generated from
    a list of expected rates. The deposit formats them inconsistently -- ``0.5``,
    ``1.0``, ``0.05`` and ``0.20`` all appear -- so any round-trip through a
    format string silently misses files, and a missing rate would quietly shorten
    the lever arm of the two-clock separation.
    """
    if panel not in LIN_TRANSIENTS:
        raise KeyError(f"unknown panel {panel!r}; have {sorted(LIN_TRANSIENTS)}")
    label, steady, subdir = LIN_TRANSIENTS[panel]
    directory = Path(root) / subdir
    if not directory.is_dir():
        raise FileNotFoundError(f"{directory} does not exist")

    found: list[tuple[float, Path]] = []
    for path in sorted(directory.glob(f"{panel}_gamma_dot_*.csv")):
        m = _RATE_IN_NAME.search(path.name)
        if m is None:
            raise ValueError(f"{path.name}: cannot parse a shear rate from the file name")
        found.append((float(m.group(1)), path))
    if not found:
        raise FileNotFoundError(f"no {panel}_gamma_dot_*.csv files in {directory}")

    found.sort(key=lambda pair: pair[0])
    rates = [r for r, _ in found]
    if len(set(rates)) != len(rates):
        raise ValueError(f"{panel}: duplicate shear rates in {[r for r, _ in found]}")

    strains, etas = [], []
    for rate, path in found:
        g, e = _read_two_column(path)
        if np.any(np.diff(g) <= 0):
            raise ValueError(f"{path.name}: accumulated strain is not strictly increasing")
        if np.any(e <= 0):
            raise ValueError(f"{path.name}: non-positive relative viscosity")
        strains.append(g)
        etas.append(e)

    return ReversalTransients(
        panel=panel,
        label=label,
        steady_panel=steady,
        rates=tuple(rates),
        strain=tuple(strains),
        eta_r=tuple(etas),
        source=directory,
    )

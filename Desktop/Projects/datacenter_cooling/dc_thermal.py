"""
dc_thermal.py — reduced-order thermal + energy model for an air-cooled,
raised-floor data hall (hot-aisle / cold-aisle).

Tier-0 of the data-center-cooling optimization project: a fast, differentiable-ish,
library-light model that captures the ONE tension that governs air-cooling energy —
    raise the supply-air temperature to make the chiller efficient,
    but hot-air *recirculation* then pushes the hottest rack past its ASHRAE inlet limit.
CFD (OpenFOAM) comes later to calibrate the recirculation matrix; this model is the
thing you actually optimize over and hand to a company as a first concrete artifact.

Physics, reduced to a linear solve:
  Rack air-temperature rise (server ΔT):      d_i = Q_i / (rho * cp * q_i)
  Recirculation influence matrix B (row sum = recirculated fraction of inlet air):
      T_in  = T_s + B (T_ex - T_s),   T_ex = T_in + d
   => (I - B)(T_in - T_s) = B d   =>   T_in = T_s + (I - B)^{-1} B d
Perfect containment (eta=1) => B = 0 => every rack inlet = supply temp.

References the model is built to match (open literature; verify before formal use):
  VanGilder & Shrivastava, influence/capture-index reduced-order DC models.
  ASHRAE TC 9.9 thermal guidelines (recommended inlet <= 27 C; allowable A1 <= 32 C).
  Chiller COP rises with evaporator/supply temperature; fan power ~ flow^3 (affinity law).
"""
from __future__ import annotations
import numpy as np

# ---- air properties (dry air, ~25 C) ----
RHO = 1.18       # kg/m^3
CP = 1006.0      # J/(kg K)

# ---- recirculation vs containment eta ----
# 'cfd': calibrated 2026-07-23 from the de Vahl Davis-validated buoyant solver in cfd/
#        (see cfd/cfd_calibration.json). phi(eta) = floor + amp*exp(-eta/tau).
#        Key finding: phi does NOT vanish at eta=1 (a ~7% recirculation floor), and the
#        falloff is concave -- the last increments of containment barely help.
# 'linear': the original placeholder base_recirc*(1-eta) -> 0 at eta=1 (over-credits containment).
CFD_RECIRC = dict(floor=0.007, amp=0.228, tau=0.790)

# ---- recirculation vs containment AND airflow provisioning ----
# 'cfd2d_p': calibrated 2026-08-31 from cfd/cfd_provisioning.json -- the same solver,
#        now with a forced through-flow (floor supply grille under the cold aisle,
#        ceiling return over the hot aisle). 48 cells, 0 diverged. This is the fix for
#        DESIGN.md Sec. 8 gap (1): over-provisioning DOES suppress recirculation, by a
#        factor of 5.5 between p = 0 and p = 1 at zero containment, where the model
#        above credits it with nothing at all.
#
# The table is the measurement, not a fit. A three-parameter dilution law
# (phi0 exp(-eta/tau) / (1 + p/k) + floor) was tried and is recorded in
# cfd/cfd_provisioning_fit.json; it is not used, because its `floor` fits to its own
# lower bound and its worst cell is off by four orders of magnitude. The floor it is
# trying to absorb is eta-dependent -- at eta = 0 there is no partition and the hot
# wall conducts straight to the inlet face -- and one additive constant cannot be
# both that and zero.
#
# TWO THINGS THIS TABLE IS NOT
#   1. Its `p` is anchored on the CELL's own buoyant circulation: p = flux / Q_1 with
#      Q_1 = max|psi| of the buoyancy-only solution at eta = 0. Identifying that with
#      the hall's provisioning ratio is an ASSUMPTION and it has not been tested.
#      What transfers is the SHAPE of phi(p) -- that it falls steeply and then
#      saturates -- not the location of p = 1. `P_ANCHOR_SCALE` exists so the
#      identification can be moved and the conclusions re-checked; it is not a knob to
#      tune against an answer.
#   2. It is a Ra = 4e3 laminar cavity. A hall is Ra ~ 1e9-1e12. cfd/run_floor_probe.py
#      shows the eta -> 1 "7% recirculation floor" in CFD_RECIRC is thermal DIFFUSION,
#      not recirculated air: sweeping the tracer's diffusivity at eta = 0.95 moves the
#      signal from 0.0025 (Le = 50) to 0.06945 (Le = 1), where it reproduces the
#      temperature channel to five decimals, while the same sweep at eta = 0 moves it
#      by 0.4%. An advected signal is Lewis-independent; a diffused one is proportional
#      to the diffusivity. So the floor -- and the "containment ROI flattens" reading
#      built on it -- is a low-Rayleigh artefact and does not belong in a hall.
CFD_P_ETA = (0.0, 0.2, 0.4, 0.6, 0.8, 0.95)
CFD_P_P = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
CFD_P_PHI = (
    (0.2411, 0.1216, 0.0862, 0.0435, 0.0272, 0.0226, 0.0211, 0.0213),
    (0.1897, 0.0772, 0.0496, 0.0179, 0.0063, 0.0028, 0.0011, 0.0006),
    (0.1382, 0.0582, 0.0332, 0.0091, 0.0022, 0.0006, 0.0001, 0.0000),
    (0.1163, 0.0530, 0.0272, 0.0060, 0.0012, 0.0003, 0.0000, 0.0000),
    (0.1007, 0.0385, 0.0147, 0.0020, 0.0003, 0.0001, 0.0000, 0.0000),
    (0.0729, 0.0144, 0.0030, 0.0002, 0.0000, 0.0000, 0.0000, 0.0000),
)
P_ANCHOR_SCALE = 1.0   # hall p -> cavity p. See note 1 above. Not a tuning knob.

# ---- server-fan ramp ----
# Server fans hold a floor speed until the inlet reaches a knee near 25 C, then
# ramp to hold junction temperature. Published measurements put the room-wide
# sensitivity at 0.35-0.5 % of server power per K in the 23-30 C ASHRAE band, and
# 7-20 % total over 15 -> 35 C (which is 0.35-1.0 %/K averaged over that wider
# span, consistent with a knee partway through it). 0.4 %/K is taken as the
# mid-band value; T_RAMP is the knee.
#
# This term matters because it lands INSIDE the IT denominator of PUE. Cutting
# CRAH airflow raises inlet temperatures, server fans take up the work, and PUE
# improves while the hall burns more. Any model that optimizes PUE without it will
# recommend under-blowing. Sources are listed in DESIGN.md Sec. 8.
#
# NOT to be confused with ASHRAE's x-factor, which is a hardware FAILURE-RATE
# multiplier versus inlet temperature, not a power term.
SERVER_FAN_RAMP = 0.004      # fraction of rack IT power per K above T_RAMP
T_RAMP = 25.0                # C, inlet temperature at which server fans start ramping


def recirc_fraction(eta, model="cfd", base_recirc=0.28, provisioning=1.0):
    """Geometric recirculated fraction at containment ``eta``.

    ``provisioning`` is used only by the ``cfd2d_p`` model; the older two ignore it,
    which is precisely the defect that model exists to fix.
    """
    eta = float(np.clip(eta, 0.0, 1.0))
    if model == "cfd2d_p":
        return _recirc_grid(eta, float(provisioning))
    if model == "cfd":
        c = CFD_RECIRC
        return c["floor"] + c["amp"] * np.exp(-eta / c["tau"])
    return base_recirc * (1.0 - eta)          # legacy linear placeholder


def _recirc_grid(eta, p):
    """Bilinear interpolation of the measured phi(eta, p) table.

    Clamped, not extrapolated, on both axes. Above p = 4 the table is already flat
    (0.0213 at eta = 0, 0.0000 elsewhere), so holding it is the measurement rather
    than an assumption; below eta = 0 and above eta = 0.95 there is nothing to
    extrapolate into, and mass balance -- not this table -- is what sets the small-p
    behaviour in `influence_matrix`.
    """
    e = np.clip(eta, CFD_P_ETA[0], CFD_P_ETA[-1])
    q = np.clip(p * P_ANCHOR_SCALE, CFD_P_P[0], CFD_P_P[-1])
    ie = int(np.searchsorted(CFD_P_ETA, e, side="right") - 1)
    ie = min(max(ie, 0), len(CFD_P_ETA) - 2)
    iq = int(np.searchsorted(CFD_P_P, q, side="right") - 1)
    iq = min(max(iq, 0), len(CFD_P_P) - 2)
    te = (e - CFD_P_ETA[ie]) / (CFD_P_ETA[ie + 1] - CFD_P_ETA[ie])
    tq = (q - CFD_P_P[iq]) / (CFD_P_P[iq + 1] - CFD_P_P[iq])
    g = CFD_P_PHI
    return float(
        (1 - te) * ((1 - tq) * g[ie][iq] + tq * g[ie][iq + 1])
        + te * ((1 - tq) * g[ie + 1][iq] + tq * g[ie + 1][iq + 1])
    )


class Hall:
    """A generic data hall: rack loads, geometry-driven recirculation, cooling plant."""

    def __init__(self, rack_kw, dT_rack=13.0, base_recirc=0.28, coupling_len=1.5,
                 containment=0.0, cop_ref=5.0, cop_slope=0.15, T_ref=18.0,
                 fan_kw_rated=None, T_allow=27.0, p_min=0.6, recirc_model="cfd",
                 fan_ramp=SERVER_FAN_RAMP, T_ramp=T_RAMP):
        self.recirc_model = recirc_model               # 'cfd' (calibrated) or 'linear' (placeholder)
        self.fan_ramp = float(fan_ramp)                # server-fan power, fraction of IT per K
        self.T_ramp = float(T_ramp)                    # inlet temp at which server fans ramp
        self.Q = np.asarray(rack_kw, float)            # per-rack IT load [kW]
        self.n = self.Q.size
        self.dT_rack = float(dT_rack)                  # server air-side rise [K]
        self.base_recirc = float(base_recirc)          # recirc fraction at zero containment
        self.coupling_len = float(coupling_len)        # rack-spacings over which exhaust couples
        self.containment = float(containment)          # eta in [0,1]; 1 = perfect containment
        self.cop_ref = float(cop_ref)                  # chiller COP at T_ref
        self.cop_slope = float(cop_slope)              # dCOP/dT_supply [1/K]
        self.T_ref = float(T_ref)
        self.T_allow = float(T_allow)                  # max allowable rack-inlet temp [C]
        self.p_min = float(p_min)                      # min supply/demand airflow ratio
        # rated CRAH fan power at provisioning p=1: ~10% of IT as a generic default
        self.fan_kw_rated = float(fan_kw_rated if fan_kw_rated is not None
                                  else 0.10 * self.Q.sum())
        self._W = self._proximity_weights()

    # geometry: how each rack's exhaust is distributed over inlets (row-normalized)
    def _proximity_weights(self):
        idx = np.arange(self.n)
        D = np.abs(idx[:, None] - idx[None, :])
        W = np.exp(-D / self.coupling_len)             # nearer racks recirculate more
        W /= W.sum(axis=1, keepdims=True)              # each row sums to 1
        return W

    def influence_matrix(self, provisioning):
        """B: recirculated fraction reaching each inlet. Worsens if under-provisioned.
        The containment dependence phi(eta) is CFD-calibrated by default (see cfd/).

        Mass balance sets a hard floor. With p = supplied/demanded airflow, the CRAHs
        deliver only p of what the servers ingest, so a fraction (1-p) of every inlet
        stream must be made up from recirculated exhaust -- there is no third air
        source in a sealed hall. Hence rowsum(B) >= 1-p, enforced below. The geometric
        term phi(eta)*(1+under) alone violates this for p < ~0.85 (at p=0.6 it returns
        0.23 where conservation demands 0.40), which let the optimizer buy the p^3 fan
        saving without paying for the recirculation that under-blowing actually causes.
        """
        p = float(provisioning)
        under = max(0.0, 1.0 - p)                       # airflow deficit amplifies recirc
        if self.recirc_model == "cfd2d_p":
            # The p dependence is measured, so it must NOT also be applied by hand:
            # multiplying the table by (1 + under) would count under-provisioning
            # twice. The mass-balance floor still applies -- the cavity has no rack
            # fan, so nothing in the CFD knows that a hall short of supply air must
            # make the difference up from the hot aisle.
            phi = recirc_fraction(self.containment, self.recirc_model,
                                  self.base_recirc, provisioning=p)
        else:
            phi = recirc_fraction(self.containment, self.recirc_model,
                                  self.base_recirc) * (1.0 + under)
        phi = max(phi, under)                           # mass balance, see docstring
        phi = np.clip(phi, 0.0, 0.95)
        return phi * self._W

    def rack_inlet_temps(self, T_supply, provisioning):
        d = np.full(self.n, self.dT_rack)               # server rise (airflow sized to load)
        B = self.influence_matrix(provisioning)
        rise = np.linalg.solve(np.eye(self.n) - B, B @ d)
        return T_supply + rise                          # [C], per rack

    # ---- IT side ----
    def it_kw(self, T_supply, provisioning):
        """IT power including the server-fan ramp. See SERVER_FAN.

        Server fans hold a floor speed until their inlet reaches T_ramp, then speed
        up to hold junction temperature. That power is inside the IT denominator of
        PUE, so a model that omits it can lower PUE by making the racks hotter --
        which is exactly what the airflow variable does.
        """
        Tin = self.rack_inlet_temps(T_supply, provisioning)
        excess = np.maximum(0.0, Tin - self.T_ramp)      # per rack, K above the knee
        return float(np.sum(self.Q * (1.0 + self.fan_ramp * excess)))

    # ---- cooling plant ----
    def cop(self, T_supply):
        return max(1.5, self.cop_ref + self.cop_slope * (T_supply - self.T_ref))

    def fan_kw(self, provisioning):
        return self.fan_kw_rated * provisioning ** 3    # affinity law

    def chiller_kw(self, T_supply, provisioning=1.0):
        # Heat removed is the FULL IT draw, server-fan ramp included: that power is
        # dissipated in the hall like any other watt.
        return self.it_kw(T_supply, provisioning) / self.cop(T_supply)

    def cooling_kw(self, T_supply, provisioning):
        pumps = 0.02 * self.Q.sum()                     # small ~constant parasitic
        return (self.fan_kw(provisioning)
                + self.chiller_kw(T_supply, provisioning) + pumps)

    def pue(self, T_supply, provisioning):
        """Total facility power / IT power.

        Both terms move with the server-fan ramp, and they move the ratio in
        OPPOSITE directions -- which is the whole problem with PUE as an objective.
        Use total_kw() to rank designs by energy.
        """
        it = self.it_kw(T_supply, provisioning)
        return (it + self.cooling_kw(T_supply, provisioning)) / it

    def total_kw(self, T_supply, provisioning):
        """Total facility power [kW]. The honest objective; PUE is a ratio."""
        return self.it_kw(T_supply, provisioning) + self.cooling_kw(T_supply, provisioning)

    def evaluate(self, T_supply, provisioning):
        Tin = self.rack_inlet_temps(T_supply, provisioning)
        it = self.it_kw(T_supply, provisioning)
        return {
            "T_supply": T_supply, "provisioning": provisioning,
            "rack_inlet_C": Tin, "hottest_C": float(Tin.max()),
            "margin_C": float(self.T_allow - Tin.max()),
            "COP": self.cop(T_supply),
            "fan_kW": self.fan_kw(provisioning),
            "chiller_kW": self.chiller_kw(T_supply, provisioning),
            "cooling_kW": self.cooling_kw(T_supply, provisioning),
            "IT_nominal_kW": float(self.Q.sum()),
            "IT_kW": it,
            "server_fan_kW": it - float(self.Q.sum()),
            "total_kW": self.total_kw(T_supply, provisioning),
            "PUE": self.pue(T_supply, provisioning),
        }


def generic_hall(n_racks=20, total_kw=200.0, hot_spots=True, seed=0, **kw):
    """A generic mixed-density hall: uniform base load with a couple of hot rows."""
    rng = np.random.default_rng(seed)
    base = np.full(n_racks, total_kw / n_racks)
    if hot_spots:
        load = base * rng.uniform(0.6, 1.0, n_racks)   # some racks lightly loaded
        hot = rng.choice(n_racks, size=max(2, n_racks // 6), replace=False)
        load[hot] *= 1.8                                # high-density pockets
        load *= total_kw / load.sum()                  # renormalize to target total
    else:
        load = base
    return Hall(load, **kw)

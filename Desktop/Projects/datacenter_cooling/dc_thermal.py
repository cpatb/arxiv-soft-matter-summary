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


def recirc_fraction(eta, model="cfd", base_recirc=0.28):
    eta = float(np.clip(eta, 0.0, 1.0))
    if model == "cfd":
        c = CFD_RECIRC
        return c["floor"] + c["amp"] * np.exp(-eta / c["tau"])
    return base_recirc * (1.0 - eta)          # legacy linear placeholder


class Hall:
    """A generic data hall: rack loads, geometry-driven recirculation, cooling plant."""

    def __init__(self, rack_kw, dT_rack=13.0, base_recirc=0.28, coupling_len=1.5,
                 containment=0.0, cop_ref=5.0, cop_slope=0.15, T_ref=18.0,
                 fan_kw_rated=None, T_allow=27.0, p_min=0.6, recirc_model="cfd"):
        self.recirc_model = recirc_model               # 'cfd' (calibrated) or 'linear' (placeholder)
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
        phi = recirc_fraction(self.containment, self.recirc_model, self.base_recirc) * (1.0 + under)
        phi = max(phi, under)                           # mass balance, see docstring
        phi = np.clip(phi, 0.0, 0.95)
        return phi * self._W

    def rack_inlet_temps(self, T_supply, provisioning):
        d = np.full(self.n, self.dT_rack)               # server rise (airflow sized to load)
        B = self.influence_matrix(provisioning)
        rise = np.linalg.solve(np.eye(self.n) - B, B @ d)
        return T_supply + rise                          # [C], per rack

    # ---- cooling plant ----
    def cop(self, T_supply):
        return max(1.5, self.cop_ref + self.cop_slope * (T_supply - self.T_ref))

    def fan_kw(self, provisioning):
        return self.fan_kw_rated * provisioning ** 3    # affinity law

    def chiller_kw(self, T_supply):
        return self.Q.sum() / self.cop(T_supply)        # heat removed = IT load

    def cooling_kw(self, T_supply, provisioning):
        pumps = 0.02 * self.Q.sum()                     # small ~constant parasitic
        return self.fan_kw(provisioning) + self.chiller_kw(T_supply) + pumps

    def pue(self, T_supply, provisioning):
        return 1.0 + self.cooling_kw(T_supply, provisioning) / self.Q.sum()

    def evaluate(self, T_supply, provisioning):
        Tin = self.rack_inlet_temps(T_supply, provisioning)
        return {
            "T_supply": T_supply, "provisioning": provisioning,
            "rack_inlet_C": Tin, "hottest_C": float(Tin.max()),
            "margin_C": float(self.T_allow - Tin.max()),
            "COP": self.cop(T_supply),
            "fan_kW": self.fan_kw(provisioning),
            "chiller_kW": self.chiller_kw(T_supply),
            "cooling_kW": self.cooling_kw(T_supply, provisioning),
            "IT_kW": float(self.Q.sum()),
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

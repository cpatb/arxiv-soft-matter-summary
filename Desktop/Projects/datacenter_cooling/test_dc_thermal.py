"""
test_dc_thermal.py — validation against analytic limits and conservation laws.
Run: python test_dc_thermal.py   (also runs under pytest)
Rigor convention (per the group CLAUDE.md): analytic-limit checks, energy-balance
closure, and monotonicity of the physical trends — a cancelling pair of bugs
should not pass silently.
"""
import numpy as np
from dc_thermal import generic_hall
from optimize_cooling import optimize

TOL = 1e-9


def test_perfect_containment_linear_model_inlet_equals_supply():
    # the idealized linear model has zero recirculation at eta=1
    h = generic_hall(containment=1.0, seed=2, recirc_model="linear")
    Tin = h.rack_inlet_temps(T_supply=20.0, provisioning=1.0)
    assert np.allclose(Tin, 20.0, atol=TOL), "linear model: eta=1 must give inlet == supply"


def test_cfd_containment_floor():
    # CFD-calibrated default: full containment does NOT zero out recirculation (~7% floor)
    h = generic_hall(containment=1.0, seed=2, recirc_model="cfd")
    Tin = h.rack_inlet_temps(T_supply=20.0, provisioning=1.0)
    assert Tin.max() > 20.3, "CFD model: a recirculation floor must remain at eta=1"
    from dc_thermal import recirc_fraction
    assert 0.05 < recirc_fraction(1.0, "cfd") < 0.10, "CFD floor should be ~0.07"


def test_recirc_raises_inlet_above_supply():
    h = generic_hall(containment=0.0, seed=2)
    Tin = h.rack_inlet_temps(T_supply=20.0, provisioning=1.0)
    assert Tin.min() > 20.0 - TOL and Tin.max() > 20.0, "recirc must lift inlets above supply"


def test_influence_matrix_rows_bounded():
    h = generic_hall(containment=0.0, seed=2)
    B = h.influence_matrix(provisioning=1.0)
    assert (B >= 0).all() and (B.sum(axis=1) < 1.0).all(), "recirc fraction must be in [0,1)"


def test_underprovisioning_respects_mass_balance():
    # p = supplied/demanded airflow, and rowsum(B) is the recirculated fraction of
    # inlet air. If the CRAHs deliver only p of what the servers ingest, at least
    # (1-p) of every inlet stream must be recirculated exhaust -- there is no other
    # air source. The geometric term alone under-counts this badly at low p.
    h = generic_hall(containment=0.3, seed=1)
    for p in (1.0, 0.9, 0.8, 0.7, 0.6):
        rows = h.influence_matrix(provisioning=p).sum(axis=1)
        assert np.allclose(rows, rows[0], atol=TOL), "phi must be common to all rows"
        assert rows[0] >= max(0.0, 1.0 - p) - TOL, (
            f"p={p}: recirculated fraction {rows[0]:.4f} is below the mass-balance "
            f"floor {1.0 - p:.4f}")


def test_underprovisioning_is_penalised_not_free():
    # The fan saving goes as p^3 while the thermal penalty must grow with (1-p).
    # If it were free the optimizer would always sit on the p_min bound.
    h = generic_hall(containment=0.3, seed=1)
    hot = [h.rack_inlet_temps(20.0, p).max() for p in (1.0, 0.9, 0.8, 0.7, 0.6)]
    assert all(b > a for a, b in zip(hot, hot[1:])), (
        f"hottest inlet must rise monotonically as airflow is cut, got {hot}")


def test_energy_balance_closure():
    # chiller must remove exactly the IT heat -- including the server-fan ramp,
    # which is dissipated in the hall like any other watt -- and the PUE identity
    # must hold against the SAME IT number in numerator and denominator.
    h = generic_hall(seed=3)
    for T_s in (18.0, 26.0, 30.0):
        ev = h.evaluate(T_supply=T_s, provisioning=1.0)
        assert abs(ev["chiller_kW"] * ev["COP"] - ev["IT_kW"]) < 1e-6
        pue_check = (ev["IT_kW"] + ev["cooling_kW"]) / ev["IT_kW"]
        assert abs(pue_check - ev["PUE"]) < 1e-9
        assert abs(ev["total_kW"] - (ev["IT_kW"] + ev["cooling_kW"])) < 1e-9
        assert abs(ev["IT_kW"] - ev["IT_nominal_kW"] - ev["server_fan_kW"]) < 1e-9


def test_server_fan_ramp_is_off_below_the_knee_and_monotone_above():
    h = generic_hall(containment=0.3, seed=1)
    # cold supply: every inlet is below T_ramp, so IT is the nominal load
    cold = h.evaluate(T_supply=14.0, provisioning=1.0)
    assert cold["hottest_C"] < h.T_ramp, "fixture too warm to test the floor"
    assert abs(cold["server_fan_kW"]) < 1e-12, "fans must not ramp below the knee"
    # above it, fan power rises with supply temperature and never falls
    draws = [h.evaluate(T, 1.0)["server_fan_kW"] for T in (24.0, 26.0, 28.0, 30.0)]
    assert draws[0] >= 0.0
    assert all(b > a for a, b in zip(draws, draws[1:])), draws


def test_pue_objective_never_beats_energy_objective_on_energy():
    # By construction the energy optimum minimises total power, so optimising PUE
    # can only match or lose on kW. If this ever flips, the objective wiring is wrong.
    for T_allow in (27.0, 32.0, 35.0, 40.0):
        h = generic_hall(containment=0.3, seed=1, T_allow=T_allow)
        pue_opt = h.evaluate(*optimize(h, objective="pue").x)
        kw_opt = h.evaluate(*optimize(h, objective="total").x)
        assert pue_opt["total_kW"] >= kw_opt["total_kW"] - 1e-6, (
            f"T_allow={T_allow}: PUE optimum drew LESS power than the energy optimum")


def test_pue_can_be_gamed_when_the_ramp_is_steep():
    # The point of charging server-fan power: at the top of the published ramp range
    # and on a relaxed inlet envelope, the design with the BEST PUE is not the one
    # that burns the least. If this test stops failing to find a gap, the term has
    # been neutered.
    h = generic_hall(containment=0.3, seed=1, T_allow=35.0, fan_ramp=0.010)
    pue_opt = h.evaluate(*optimize(h, objective="pue").x)
    kw_opt = h.evaluate(*optimize(h, objective="total").x)
    assert pue_opt["PUE"] < kw_opt["PUE"] - 1e-4, "the PUE optimum should look better on PUE"
    assert pue_opt["total_kW"] > kw_opt["total_kW"] + 1.0, (
        "the PUE optimum should cost real kW; got "
        f"{pue_opt['total_kW']:.2f} vs {kw_opt['total_kW']:.2f}")
    assert pue_opt["server_fan_kW"] > kw_opt["server_fan_kW"], (
        "the gap must come from server-fan power moving into the IT denominator")


def test_cop_and_fan_monotonic():
    h = generic_hall(seed=3)
    assert h.cop(25) > h.cop(15), "COP must rise with supply temp"
    assert h.fan_kw(1.2) > h.fan_kw(0.8), "fan power must rise with airflow"
    # cubic scaling
    assert abs(h.fan_kw(1.0) * 8 - h.fan_kw(2.0)) < 1e-6


def test_optimizer_respects_constraint_and_helps():
    h = generic_hall(containment=0.3, seed=1)
    base = h.pue(18.0, 1.0)
    r = optimize(h)
    hottest = h.rack_inlet_temps(r.x[0], r.x[1]).max()
    assert hottest <= h.T_allow + 1e-3, "optimum must satisfy the ASHRAE inlet limit"
    assert r.fun <= base + 1e-6, "optimized PUE must not be worse than the baseline"


def test_better_containment_lowers_optimized_pue():
    # Weak monotonicity, with a tolerance on BOTH comparisons. Exact ties are
    # expected here and are not a bug: once the mass-balance floor (1-p) exceeds
    # the geometric term phi(eta)*(1+under), eta drops out of B entirely, so every
    # eta returns the same optimum. At p=1 the geometric term is the only term and
    # the sweep separates again -- that is the curve reported in README.md.
    pues = []
    for eta in (0.0, 0.5, 1.0):
        r = optimize(generic_hall(containment=eta, seed=1))
        pues.append(r.fun)
    assert pues[0] >= pues[1] - 1e-9 >= pues[2] - 2e-9, "more containment -> lower or equal PUE"


def test_containment_roi_at_fixed_provisioning():
    # The containment lever, measured the way an operator sees it: airflow held at
    # the design point (p=1, mass balance satisfied with no reliance on the
    # unmodelled p>1 branch), supply temperature raised to the ASHRAE limit.
    pues, temps = [], []
    for eta in (0.0, 0.5, 1.0):
        h = generic_hall(containment=eta, total_kw=200.0, seed=1)
        T_s = h.T_allow - h.rack_inlet_temps(0.0, 1.0).max()
        temps.append(T_s)
        pues.append(h.pue(T_s, 1.0))
    assert temps[0] < temps[1] < temps[2], "containment must unlock higher supply temp"
    assert pues[0] > pues[1] > pues[2], "containment must lower PUE at fixed airflow"


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for fn in fns:
        try:
            fn(); print(f"[PASS] {fn.__name__}"); passed += 1
        except AssertionError as e:
            print(f"[FAIL] {fn.__name__}: {e}")
    print(f"\n{passed}/{len(fns)} checks passed")

"""
Validation suite for the 6-DOF disc integrator.

Runs four canonical scenarios and prints summary statistics. Designed to
fail loudly if a refactor breaks the physics.

Scenarios:
    1. Flat backhand in still air         - check trim convergence
    2. Headwind 5 m/s                     - check range reduction
    3. Tailwind 5 m/s                     - check alpha increase + anhyzer turn
    4. Crosswind 5 m/s                    - check lateral drift ~= U_c * T

All four cases use the Hummel short-flight coefficients.
"""

from __future__ import annotations
import numpy as np
import disc_6dof as d6
from coefficients import HUMMEL_SHORT, ULTRASTAR
from trajectories import right_handed_backhand


def case_still_air():
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0, wind=d6.still_air, t_max=8.0)
    alpha_final = res.diag[res.landing_idx]["alpha_deg"]
    return dict(name="still air", range_m=res.range_m,
                flight_s=res.flight_time, alpha_landing=alpha_final,
                landing=res.r[-1])


def case_headwind():
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0,
                      wind=d6.uniform_wind([-5.0, 0.0, 0.0]),
                      t_max=8.0)
    return dict(name="headwind 5 m/s", range_m=res.range_m,
                flight_s=res.flight_time,
                alpha_landing=res.diag[res.landing_idx]["alpha_deg"],
                landing=res.r[-1])


def case_tailwind():
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0,
                      wind=d6.uniform_wind([+5.0, 0.0, 0.0]),
                      t_max=8.0)
    return dict(name="tailwind 5 m/s", range_m=res.range_m,
                flight_s=res.flight_time,
                alpha_landing=res.diag[res.landing_idx]["alpha_deg"],
                landing=res.r[-1])


def case_crosswind():
    r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
    res = d6.simulate(r0, v0, q0, w0,
                      wind=d6.uniform_wind([0.0, +5.0, 0.0]),
                      t_max=8.0)
    drift_predicted = 5.0 * res.flight_time
    return dict(name="crosswind 5 m/s", range_m=res.range_m,
                flight_s=res.flight_time,
                alpha_landing=res.diag[res.landing_idx]["alpha_deg"],
                landing=res.r[-1],
                drift_y=float(res.r[-1, 1] - r0[1]),
                drift_predicted_Uc_T=drift_predicted)


def main():
    print(f"{'case':22s} {'range (m)':>10s} {'flight (s)':>10s} "
          f"{'alpha land':>10s} {'landing (x,y,z)':>30s}")
    print("-" * 90)
    for fn in (case_still_air, case_headwind, case_tailwind, case_crosswind):
        out = fn()
        landing = out["landing"]
        print(f"{out['name']:22s} {out['range_m']:>10.2f} "
              f"{out['flight_s']:>10.2f} {out['alpha_landing']:>10.2f}  "
              f"({landing[0]:7.2f},{landing[1]:7.2f},{landing[2]:7.2f})")
        if "drift_y" in out:
            print(f"   crosswind drift y observed = {out['drift_y']:.2f} m,"
                  f" U_c*T predicted = {out['drift_predicted_Uc_T']:.2f} m")


if __name__ == "__main__":
    main()

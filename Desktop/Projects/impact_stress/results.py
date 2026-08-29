"""
results.py — produce every number and figure in the report, from one run.

    python results.py            # writes out/results.json + out/results.txt + figures/

Nothing in the report is typed by hand; everything traces to a line here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import conventions as cv
import dilatancy as dl
import inertial_dilatant as idn
import materials as mat
import stress as st
from dilatancy import T_STAR_EVAL

OUT = Path(__file__).parent / "out"
FIG = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

# --------------------------------------------------------------------------
# The lab's experiment, as recorded in the status cards and repo configs
# --------------------------------------------------------------------------
LAB = dict(
    D0=2.5e-3,                      # 2.3-2.9 mm from one needle; midpoint
    U_list=[1.15, 1.5, 2.3, 3.0, 4.0, 5.0],
    U_band_B=2.3,                   # velocity band B centre (JUN)
    U_onset_measured={              # fitted u0* from the front collapse, m/s
        "0 mM": 1.473, "10 mM": 1.120, "50 mM": 1.138, "100 mM": 0.455},
    gammadot_c_measured={           # fitted onset, 1/s, on the U0/D0 convention
        "0 mM": 631.6, "10 mM": 478.8, "50 mM": 470.6, "100 mM": 195.5},
    phi=0.44,
)
R_LAB = LAB["D0"] / 2

SUSP = mat.lab_silica(phi=LAB["phi"])
DROP = mat.suspension_as_fluid(SUSP)

# Viscosity is the largest single uncertainty; carry a bracket everywhere.
#   low  : Krieger-Dougherty with the true solvent (0.955 mPa.s) and phi0 = 0.70
#   mid  : the repo's fitted Wyart-Cates f=0 branch  (what the DNS uses)
#   high : Krieger-Dougherty with phi_m = 0.58
# Ordered low -> high, and checked at import so the labels cannot drift
# away from the values (they had, once: the "high" entry was below "mid").
ETA_BRACKET = {
    "low  (Krieger-Dougherty, true solvent 0.955 mPa.s, phi0=0.70)":
        0.955e-3 / (1 - 0.44 / 0.70) ** 2,
    "mid  (Krieger-Dougherty, true solvent 0.955 mPa.s, phi_m=0.58)":
        0.955e-3 / (1 - 0.44 / 0.58) ** 2,
    "high (repo Wyart-Cates fit, eta_s=4.74 mPa.s -- what the DNS uses)":
        SUSP.eta_0,
}
assert list(ETA_BRACKET.values()) == sorted(ETA_BRACKET.values()), \
    f"eta bracket is not ordered low->high: {ETA_BRACKET}"


def drop_with_eta(eta):
    return mat.Fluid(name=DROP.name, rho=DROP.rho, c=DROP.c, eta=eta,
                     gamma=DROP.gamma, provenance=DROP.provenance)


def imp(U, target=None, eta=None, R=R_LAB):
    return st.Impact(drop=drop_with_eta(eta or DROP.eta), U=U, R=R,
                     target=target or st.rigid_target())


# ==========================================================================
def section_1_stress_ladder():
    """Every stress scale of one impact, on one axis."""
    U = LAB["U_band_B"]
    i = imp(U)
    sd = st.shock_detachment(i)
    gf = st.gas_film(i)
    xo = st.wagner_hammer_crossover(i)
    pf = st.peak_force(i)

    rows = []

    def row(name, p, ell, t, note):
        rows.append(dict(stage=name, stress_Pa=float(p), over_rhoU2=float(p) / i.p_dyn,
                         length_m=float(ell), time_s=float(t), note=note))

    row("0 gas cushion (peak film pressure)", gf["p_gas"], gf["l_dimple"],
        gf["t_dimple_star"] * i.tau, "screens the water hammer; scaling only, no prefactor")
    row("1 water hammer (rigid wall)", sd["p"], sd["d_d"], sd["t_d"],
        f"touches {sd['volume_fraction_affected']:.2e} of the drop volume")
    row("2 Wagner pressure, on axis, t*=t*_peak", st.wagner_pressure_axis(i, pf["t_star_max"]),
        i.R, pf["t_max_s"], "(sqrt3/pi) rho U^2 t*^-1/2")
    row("2 Wagner pressure, screened rim", st.wagner_pressure_peak(i, pf["t_star_max"]),
        st.DELTA_OVER_D * float(st.contact_radius(i, pf["t_max_s"])), pf["t_max_s"],
        "5.5x the axis, at all times; ring, not spot")
    row("3 wall shear, contact-area mean", st.wall_shear_mean(i, pf["t_star_max"]),
        float(st.contact_radius(i, pf["t_max_s"])), pf["t_max_s"],
        "1.135 rho U^2 Re_R^-1/2 t*^-1/2 -- from the published drag force")
    row("3 wall shear, screened rim", st.wall_shear_peak(i, pf["t_star_max"]),
        st.DELTA_OVER_D * float(st.contact_radius(i, pf["t_max_s"])), pf["t_max_s"],
        "upper bound: depends on the empirical cut-off")
    row("ref: dynamic pressure rho U^2", i.p_dyn, i.R, i.tau, "the usual (wrong) gauge")
    row("ref: dilatancy onset sigma* (silica)", SUSP.sigma_star, 2.0 * SUSP.a, np.nan,
        "length is the particle DIAMETER 2a = 800 nm (SUSP.a is a radius); "
        "this is what all of the above must be compared against")

    return dict(
        U=U, D0=LAB["D0"],
        dimensionless=dict(Re=i.Re, Re_R=i.Re_R, We=i.We, Oh=i.Oh, Fr=i.Fr,
                           Ma=i.Ma, St_gas=i.St_gas),
        rho=i.drop.rho, eta=i.drop.eta, c_wood=i.drop.c,
        tau_R_over_U_s=i.tau, tau_D_over_U_s=i.tau_D,
        peak_force=pf, shock=sd, crossover=xo, gas=gf,
        ladder=rows,
    )


def section_2_substrates():
    """What changes when the target is not rigid."""
    U = LAB["U_band_B"]
    water = mat.WATER
    targets = {
        "rigid (idealisation)": st.rigid_target(),
        "fused-silica glass slide (real)": st.glass_target(),
        "deep water bath": st.liquid_target(water),
        "deep bath of the SAME suspension": st.liquid_target(DROP),
        "Carbopol bath, tau_y = 1 Pa": st.liquid_target(water, tau_y=1.0, name="Carbopol 0.1%"),
        "Carbopol bath, tau_y = 30 Pa": st.liquid_target(water, tau_y=30.0, name="Carbopol 1%"),
        "PDMS gel E = 100 kPa (Sun et al. 2022)": st.elastic_target(970.0, 100e3),
        "PDMS gel E = 3 MPa": st.elastic_target(1030.0, 3e6),
        "acrylic E = 3 GPa": st.elastic_target(1180.0, 3e9),
    }
    out = {}
    for name, tg in targets.items():
        i = imp(U, target=tg)
        cf = st.compliance_factors(i)
        out[name] = dict(
            applicable_law=cf["applicable_law"],
            acoustic_factor=cf["acoustic_factor"],
            inertial_factor=cf["inertial_factor"],
            Cauchy_number_C=cf["compliance_number_C"],
            yield_number_J=cf["yield_number_J"],
            p_water_hammer_Pa=cf["p_water_hammer"],
            tau_slip_ratio=st.shear_sharing_factor(i)["slip_ratio"],
        )
    # The modulus at which this impact reaches Cauchy number C = 1, i.e. the
    # substrate is soft for the whole impact. This is the module's own
    # criterion (Eq. 15), not a transcription of anyone else's worked example.
    out["_E_star_at_Cauchy_number_unity_Pa"] = float(imp(U).p_dyn)
    out["_note"] = ("Howland et al. PRL 117, 184502 (2016) give E_c = 93 kPa for "
                    "THEIR conditions (R = 1.5 mm, U = 2.36 m/s, ethanol); that is "
                    "not transferable to these drops and is quoted, not recomputed.")
    return out


def section_3_dilatancy_gauges():
    """Inertial gauge vs shear gauge, and the onset velocity each predicts."""
    lo, gm, hi = dl.sigma_star_bracket()
    sigmas = {f"lit low ({lo:.3g} Pa)": lo, f"lit geometric mean ({gm:.3g} Pa)": gm,
              f"lit high ({hi:.3g} Pa)": hi,
              f"repo fit ({SUSP.sigma_star:.3g} Pa)": SUSP.sigma_star}

    gauges = {}
    for U in LAB["U_list"]:
        i = imp(U)
        gauges[f"U={U}"] = {
            "rho_U2_Pa": i.p_dyn,
            "Lambda_P (rho U^2/sigma*)": dl.lambda_inertial(i, SUSP.sigma_star),
            f"tau_mean_Pa (t*={T_STAR_EVAL})": float(st.wall_shear_mean(i, T_STAR_EVAL)),
            "tau_mean_Pa (t*=1, sensitivity)": float(st.wall_shear_mean(i, 1.0)),
            "Lambda_F (<tau>/sigma*)": dl.lambda_shear(i, SUSP.sigma_star, T_STAR_EVAL),
            "Lambda_F (t*=1, sensitivity)": dl.lambda_shear(i, SUSP.sigma_star, 1.0),
            "frictional_area_fraction": float(
                dl.frictional_area_fraction(i, SUSP.sigma_star, T_STAR_EVAL)),
            "xi_contour": float(dl.frictional_contour(i, SUSP.sigma_star, T_STAR_EVAL)),
        }

    onsets = {}
    for sname, s in sigmas.items():
        for ename, eta in ETA_BRACKET.items():
            key = f"{sname} | {ename}"
            onsets[key] = {
                "U_c_mean_criterion": dl.onset_velocity(SUSP, R_LAB, sigma_star=s,
                                                        t_star=T_STAR_EVAL, criterion="mean",
                                                        eta=eta),
                "U_c_area_half": dl.onset_velocity(SUSP, R_LAB, sigma_star=s,
                                                   t_star=T_STAR_EVAL, criterion="area",
                                                   area_threshold=0.5, eta=eta),
                "U_c_mean_at_t_star_1": dl.onset_velocity(SUSP, R_LAB, sigma_star=s,
                                                          t_star=1.0, criterion="mean", eta=eta),
            }

    # the salt prediction: U_c ~ sigma*^(2/3)
    meas = LAB["U_onset_measured"]
    salt = {
        "measured_u0star_ratio_0_over_100mM": meas["0 mM"] / meas["100 mM"],
        "implied_sigma_star_ratio": (meas["0 mM"] / meas["100 mM"]) ** 1.5,
        "note": "U_c ~ sigma*^(2/3) so sigma*(0mM)/sigma*(100mM) = (u0*ratio)^(3/2)",
    }
    return dict(gauges=gauges, onset_velocities=onsets,
                measured_onsets=meas, salt_prediction=salt,
                sigma_star_literature=dl.SIGMA_STAR_LITERATURE)


def section_4_conventions():
    """The shear-rate convention problem, quantified."""
    U = LAB["U_band_B"]

    def eta_of(gd):
        # self-consistent Wyart-Cates viscosity at shear rate gd
        sig = np.logspace(-3, 6, 4000)
        f = SUSP.f_ss(sig)
        eta = SUSP.eta_of_f(f)
        gdg = sig / eta
        return float(np.interp(gd, gdg, eta))

    tab = cv.convention_table(DROP.rho, U, LAB["D0"], SUSP.gamma, eta_of, DROP.eta)

    # what the lab's fitted gammadot_c becomes under each convention
    amp = tab["_amplification_D_over_Lc"]
    corrected = {k: v * amp for k, v in LAB["gammadot_c_measured"].items()}

    # the rheology-side onset, two definitions
    rheo = {
        "gammadot at sigma = sigma* using eta_0 (no thickening yet)": SUSP.gammadot_c_naive,
        "gammadot at f = 0.5 (self-consistent)": SUSP.gammadot_at_f(0.5),
        "gammadot at f = 0.1": SUSP.gammadot_at_f(0.1),
        "sigma at f = 0.5 [Pa]": SUSP.sigma_at_f(0.5),
        "eta_0 [Pa.s]": SUSP.eta_0,
    }
    return dict(conventions_per_s=tab,
                lab_gammadot_c_as_fitted=LAB["gammadot_c_measured"],
                lab_gammadot_c_rescaled_to_Lc=corrected,
                rheology_side=rheo)


def section_5_ID():
    """The inertial-dilatant number, silica vs cornstarch."""
    U = LAB["U_band_B"]
    out = {"definition": "ID = rho kappa U0 / (A mu_w dphi D0), Mobaseri et al. "
                         "PRL 136, 148201 (2026), Eq. 6"}

    # 1. reproduce the paper's own system, at its own conditions
    cs = idn.cornstarch_M26(phi=0.41)
    out["cornstarch_M26_check"] = idn.report(cs, U0=5.1, D0=4.0e-3)

    # 2. the lab's silica, over a range of the one unmeasured input
    sil = {}
    for phi_c in [0.36, 0.38, 0.40, 0.42, 0.43]:
        s = idn.lab_silica_dilatant(phi=0.44, phi_c=phi_c)
        sil[f"phi_c={phi_c}"] = idn.report(s, U0=U, D0=LAB["D0"])
    out["silica_vs_phi_c"] = sil

    # 3. the headline ratio, at matched dphi
    out["ratio_at_matched_dphi"] = idn.compare(
        idn.lab_silica_dilatant(0.44, 0.42), idn.cornstarch_M26(0.39),
        U0=U, D0=LAB["D0"])

    # 4. sensitivity to the friction angle
    out["beta0_sensitivity"] = {
        f"beta0={b} deg": idn.ID(idn.lab_silica_dilatant(0.44, 0.42, beta0_deg=b), U, LAB["D0"])
        for b in [30.0, 37.0, 44.0, 51.0]}

    # 4b. the validity criterion that survived review
    out["poroelastic_validity"] = {
        "cornstarch (M26 conditions)": idn.poroelastic_time_ratio(cs, 5.1, 4.0e-3),
        "lab silica": idn.poroelastic_time_ratio(
            idn.lab_silica_dilatant(0.44, 0.42), U, LAB["D0"]),
        "lab silica with 0.1% entrained gas": idn.poroelastic_time_ratio(
            idn.lab_silica_dilatant(0.44, 0.42), U, LAB["D0"], chi=1e-3 / 1.013e5),
    }

    # 5. phi < phi_c: no dilatant branch at all
    out["if_phi_c_is_0.47_as_Shah2022_suggests"] = idn.report(
        idn.lab_silica_dilatant(0.44, 0.47), U0=U, D0=LAB["D0"])
    return out


def section_6_inversion():
    """Turn the framework around: use the MEASURED onset to infer sigma*.

    This is the direction that is actually useful. Predicting U_c from a
    literature sigma* spans a factor of 9 and is therefore weak. But the
    SCALING U_c ~ sigma*^(2/3) is prefactor-free, so one calibration point
    fixes the constant and the rest of the salt ladder becomes a prediction
    with no further freedom.
    """
    meas = LAB["U_onset_measured"]
    out = {}

    # (a) absolute inversion, under each criterion and viscosity
    for cname in ["mean", "area"]:
        for ename, eta in ETA_BRACKET.items():
            U_ref = dl.onset_velocity(SUSP, R_LAB, sigma_star=SUSP.sigma_star,
                                      t_star=T_STAR_EVAL, criterion=cname,
                                      area_threshold=0.5, eta=eta)
            out[f"sigma_star_implied | {cname} | {ename.split('(')[0].strip()}"] = {
                k: SUSP.sigma_star * (v / U_ref) ** 1.5 for k, v in meas.items()}

    # (b) prefactor-calibrated prediction: fix the constant on 0 mM, predict the rest
    #     U_c/U_c(0mM) = (sigma*/sigma*(0mM))^(2/3), so
    #     sigma*(salt)/sigma*(0mM) = (u0*(salt)/u0*(0mM))^(3/2)
    ratio = {k: (v / meas["0 mM"]) ** 1.5 for k, v in meas.items()}
    out["sigma_star_ratio_to_0mM_prefactor_free"] = ratio
    out["prefactor_free_note"] = (
        "These ratios contain no fitted constant, no viscosity, no drop size and "
        "no choice of criterion - all of those cancel between salt levels. They "
        "are the falsifiable output: measure sigma*(salt) on a rheometer and "
        "check these four numbers."
    )

    # (c) how far above sigma* the shear gauge sits AT the measured onset
    out["Lambda_F_at_the_measured_onset"] = {
        k: dl.lambda_shear(imp(v), SUSP.sigma_star, T_STAR_EVAL) for k, v in meas.items()}
    out["frictional_area_at_the_measured_onset"] = {
        k: float(dl.frictional_area_fraction(imp(v), SUSP.sigma_star, T_STAR_EVAL))
        for k, v in meas.items()}
    return out


def section_7_ID_regime_check():
    """Does the capped ID land the regime boundary where it is observed?"""
    out = {}
    cs = idn.cornstarch_M26(phi=0.41)
    out["cornstarch_ID_vs_U0"] = {
        f"U0={u}": {"ID_eq6": idn.ID(cs, u, 4e-3), "ID_capped": idn.ID_capped(cs, u, 4e-3)}
        for u in [0.5, 1.0, 2.0, 3.0, 5.1, 7.2]}
    out["cornstarch_top_speed_is_7.2_m_s_and_Regime3_IS_observed"] = True

    sil = idn.lab_silica_dilatant(0.44, 0.42)
    out["silica_ID_vs_U0"] = {
        f"U0={u}": {"ID_eq6": idn.ID(sil, u, LAB["D0"]), "ID_capped": idn.ID_capped(sil, u, LAB["D0"])}
        for u in [1.0, 2.3, 5.0, 10.0, 20.0, 30.0, 50.0]}

    # Jerome et al. PRL 117, 098003 (2016) measured -10 to -100 kPa directly.
    out["jerome_measured_suction_Pa"] = [1e4, 1e5]
    out["cavitation_ceiling_Pa"] = idn.P_ATM
    out["comment"] = (
        "The largest pore suction ever directly measured under impact into a "
        "dense suspension is ~1e5 Pa = 1 atm, which is exactly the cavitation "
        "ceiling. That is the empirical basis for capping the Darcy pressure.")
    return out


# ==========================================================================
def figures():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "figure.dpi": 160,
                         "axes.grid": True, "grid.alpha": 0.25,
                         "savefig.bbox": "tight"})
    U = LAB["U_band_B"]
    i = imp(U)

    # ---- Fig 1: the self-similar field, pressure and shear ---------------
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.3))
    xi = np.linspace(0, 1 - st.DELTA_OVER_D, 800)
    for ts, col in zip([0.01, 0.1, 0.34, 1.0], ["#1b4a8f", "#2e8b57", "#d9822b", "#a8324a"]):
        ax[0].plot(xi, st.wagner_pressure_selfsimilar(i, ts, xi) / i.p_dyn, color=col,
                   label=f"$t^*={ts}$")
        ax[1].plot(xi, st.wall_shear_selfsimilar(i, ts, xi), color=col, label=f"$t^*={ts}$")
    ax[0].set(xlabel=r"$\xi=r/d(t)$", ylabel=r"$p/\rho U^2$", yscale="log",
              title="pressure: $(1-\\xi^2)^{-1/2}$, rim is 5.5$\\times$ the axis")
    ax[1].axhline(SUSP.sigma_star, color="k", ls="--", lw=1)
    ax[1].text(0.03, SUSP.sigma_star * 1.25, r"$\sigma^*=44$ Pa", fontsize=8)
    ax[1].set(xlabel=r"$\xi=r/d(t)$", ylabel=r"wall shear stress $\tau$  [Pa]", yscale="log",
              title=r"shear: $\xi/(1-\xi^2)$, zero on axis")
    ax[0].legend(fontsize=7)
    fig.suptitle(r"Self-similar wall stress: 800 nm silica, $\phi$=0.44, "
                 f"$U_0$={U} m/s, $D_0$={LAB['D0']*1e3:.1f} mm", fontsize=9.5, y=1.07)
    fig.savefig(FIG / "fig1_selfsimilar_fields.png")
    plt.close(fig)

    # ---- Fig 2: the stress ladder ----------------------------------------
    lad = section_1_stress_ladder()["ladder"]
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    names = [r["stage"] for r in lad]
    vals = [r["stress_Pa"] for r in lad]
    cols = ["#7a7a7a" if "ref" in n else "#1b4a8f" for n in names]
    cols[-1] = "#a8324a"
    ax.barh(range(len(vals)), vals, color=cols)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("stress [Pa]")
    ax.invert_yaxis()
    ax.set_title(f"Stress ladder of one impact ($U_0$={U} m/s), and $\\sigma^*$", fontsize=9)
    fig.savefig(FIG / "fig2_stress_ladder.png")
    plt.close(fig)

    # ---- Fig 3: the two gauges vs velocity, against sigma* ---------------
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    Us = np.linspace(0.3, 5.5, 300)
    lo, gm, hi = dl.sigma_star_bracket()
    p_dyn = np.array([imp(u).p_dyn for u in Us])
    tau_m = np.array([float(st.wall_shear_mean(imp(u), T_STAR_EVAL)) for u in Us])
    ax.fill_between(Us, lo, hi, color="#a8324a", alpha=0.15,
                    label=rf"$\sigma^*$ literature bracket, {lo:.0f}–{hi:.0f} Pa")
    ax.axhline(gm, color="#a8324a", lw=1.2,
           label=rf"$\sigma^*$ geom. mean {gm:.0f} Pa")
    ax.plot(Us, p_dyn, color="#7a7a7a", lw=2, label=r"inertial gauge $\rho U_0^2$")
    for (lab_, eta), ls in zip(ETA_BRACKET.items(), ["--", "-", ":"]):
        tm = np.array([float(st.wall_shear_mean(imp(u, eta=eta), T_STAR_EVAL)) for u in Us])
        ax.plot(Us, tm, color="#1b4a8f", ls=ls, lw=1.6,
                label=r"shear gauge $\langle\tau\rangle$, " + lab_.split("(")[0].strip())
    for j, (k, v) in enumerate(LAB["U_onset_measured"].items()):
        ax.axvline(v, color="#2e8b57", lw=0.9, alpha=0.7)
        ax.text(v, 6e4 if j % 2 == 0 else 9e3, " " + k, rotation=90, fontsize=6.5,
                color="#2e8b57", va="top")
    ax.text(0.35, 1.4, "green lines: MEASURED jamming onset, per salt level",
            fontsize=6.5, color="#2e8b57")
    ax.set(xlabel=r"impact speed $U_0$ [m/s]", ylabel="stress [Pa]", yscale="log",
           ylim=(1, 1e5),
           title="Which gauge can even cross the onset stress?")
    ax.legend(fontsize=6.5, loc="lower right", framealpha=0.95)
    fig.savefig(FIG / "fig3_gauge_comparison.png")
    plt.close(fig)
    _ = tau_m

    # ---- Fig 4: frictional map in the drop -------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(9.2, 3.0), sharey=True)
    for axx, Uv in zip(axs, [1.0, 2.3, 4.5]):
        ii = imp(Uv)
        ts = np.logspace(-2.5, 0.3, 260)
        xic = dl.frictional_contour(ii, SUSP.sigma_star, ts)
        d = st.contact_radius(ii, ts * ii.tau) * 1e3
        axx.fill_between(ts, xic * d, d, color="#a8324a", alpha=0.35, label=r"$\tau>\sigma^*$")
        axx.plot(ts, d, color="k", lw=1.4, label=r"contact radius $d(t)$")
        axx.plot(ts, xic * d, color="#a8324a", lw=1.4, label=r"$\tau=\sigma^*$ contour")
        axx.set(xscale="log", xlabel=r"$t^*=U_0t/R$",
                title=f"$U_0$={Uv} m/s")
    axs[0].set_ylabel("radius [mm]")
    axs[0].legend(fontsize=6.5, loc="upper left")
    fig.suptitle(r"Where in the contact patch is the suspension frictional? "
                 r"($\sigma^*=44$ Pa)", fontsize=9.5, y=1.10)
    fig.savefig(FIG / "fig4_frictional_annulus.png")
    plt.close(fig)

    # ---- Fig 5: ID vs particle size --------------------------------------
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    a_um = np.logspace(-1, 1.6, 200)
    for dphi, col in zip([0.01, 0.02, 0.05], ["#1b4a8f", "#2e8b57", "#d9822b"]):
        ids, idcap = [], []
        for a in a_um * 1e-6:
            s = idn.DilatantSystem(name="s", phi=0.44, phi_c=0.44 - dphi, a=a, rho=1440.0)
            ids.append(idn.ID(s, U, LAB["D0"]))
            idcap.append(idn.ID_capped(s, U, LAB["D0"]))
        ax.plot(a_um, ids, color=col, lw=1.8, label=rf"Eq. 6, $\Delta\phi$={dphi}")
        ax.plot(a_um, idcap, color=col, lw=1.2, ls="--")
    ax.axhline(1.0, color="k", lw=1, ls="-.")
    ax.text(0.11, 1.3, "ID = 1: liquid-like above, solid-like below", fontsize=7)
    ax.axvline(0.8, color="#a8324a", lw=1.2)
    ax.text(0.86, 1e-7, "lab silica\n800 nm", fontsize=7, color="#a8324a")
    ax.axvline(10.0, color="#7a7a7a", lw=1.2)
    ax.text(9.2, 3e-8, "cornstarch\n10 µm", fontsize=7, color="#7a7a7a", ha="right")
    ax.set(xscale="log", yscale="log", xlabel=r"particle diameter $a$ [µm]",
           ylabel="ID", ylim=(1e-8, 1e2),
           title=r"ID $\propto\kappa\propto a^2$   (dashed: capped at the cavitation limit)")
    ax.legend(fontsize=7, loc="upper left", framealpha=0.95)
    fig.savefig(FIG / "fig5_ID_vs_particle_size.png")
    plt.close(fig)

    return sorted(p.name for p in FIG.glob("*.png"))


# ==========================================================================
def main():
    res = {
        "1_stress_ladder": section_1_stress_ladder(),
        "2_substrates": section_2_substrates(),
        "3_dilatancy_gauges": section_3_dilatancy_gauges(),
        "4_shear_rate_conventions": section_4_conventions(),
        "5_inertial_dilatant_number": section_5_ID(),
        "6_inversion_to_sigma_star": section_6_inversion(),
        "7_ID_regime_check": section_7_ID_regime_check(),
        "_inputs": {"LAB": {k: v for k, v in LAB.items()},
                    "eta_bracket_Pa_s": ETA_BRACKET,
                    "suspension": {"name": SUSP.name, "phi": SUSP.phi,
                                   "sigma_star_Pa": SUSP.sigma_star,
                                   "eta_0_Pa_s": SUSP.eta_0, "rho": SUSP.rho,
                                   "c_wood_m_s": SUSP.c_wood,
                                   "caveats": list(SUSP.caveats)}},
    }
    res["_figures"] = figures()

    def default(o):
        if isinstance(o, (np.floating, np.integer)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    (OUT / "results.json").write_text(json.dumps(res, indent=2, default=default),
                                      encoding="utf-8")
    print(json.dumps(res, indent=2, default=default))


if __name__ == "__main__":
    main()

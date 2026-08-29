"""Second correction pass: evaluation time, eta bracket, geometric mean, citations."""
import pathlib

root = pathlib.Path(__file__).parent.parent


def patch(name, pairs, require_all=True):
    p = root / name
    t = p.read_text(encoding="utf-8")
    missed = [a[:70] for a, b in pairs if a not in t]
    for a, b in pairs:
        t = t.replace(a, b)
    p.write_text(t, encoding="utf-8")
    print(f"{name}: missed {missed or 'none'}")
    if require_all and missed:
        raise SystemExit(1)


# --------------------------------------------------------------- dilatancy.py
patch("dilatancy.py", [
    # honest geometric mean
    ('The first three are external to this project and do not depend on the\n'
     'retracted rheology file. Their geometric mean is 55 Pa; the repo\'s own',
     'The first three are external to this project and do not depend on the\n'
     'retracted rheology file. Their geometric mean is 49.2 Pa; the repo\'s own'),
    # a named, justified evaluation time instead of a magic 1.0
    ('from stress import DELTA_OVER_D, Impact',
     'from stress import DELTA_OVER_D, Impact\n'
     '\n'
     'T_STAR_EVAL = 0.315\n'
     '"""Default dimensionless time at which the gauges are evaluated.\n'
     '\n'
     'WHY NOT t* = 1.  The Wagner self-similar solution is derived for short\n'
     'times and PLJ16 report departure from similarity once t* = U t/R is of\n'
     'order 1. At t* = 1 the contact radius is d = sqrt(3) R = 1.73 R, larger\n'
     'than the drop itself, so the solution is being used outside its domain.\n'
     '\n'
     'WHY 0.315.  That is the PEAK-FORCE time for the lab condition (CSG22\'s\n'
     'ttil_max = 0.157 on the diameter clock, doubled), i.e. the instant of\n'
     'maximum loading, and there d = 0.97 R -- just inside validity.\n'
     '\n'
     'The choice matters: moving from t* = 1 to t* = 0.315 raises every stress\n'
     'by 1/sqrt(0.315) = 1.78 and lowers every onset velocity by 0.315^(1/3) =\n'
     '0.68. That factor is part of the honest uncertainty and is reported as a\n'
     'sensitivity rather than buried in a default."""'),
])

# --------------------------------------------------------------- results.py
patch("results.py", [
    # eta bracket: order it correctly and label what each is
    ('ETA_BRACKET = {\n'
     '    "low  (solvent 0.955 mPa.s, phi0=0.70)": 0.955e-3 / (1 - 0.44 / 0.70) ** 2,\n'
     '    "mid  (repo Wyart-Cates fit, eta_s=4.74 mPa.s)": SUSP.eta_0,\n'
     '    "high (solvent 0.955 mPa.s, phi_m=0.58)": 0.955e-3 / (1 - 0.44 / 0.58) ** 2,\n'
     '}',
     '# Ordered low -> high, and checked at import so the labels cannot drift\n'
     '# away from the values (they had, once: the "high" entry was below "mid").\n'
     'ETA_BRACKET = {\n'
     '    "low  (Krieger-Dougherty, true solvent 0.955 mPa.s, phi0=0.70)":\n'
     '        0.955e-3 / (1 - 0.44 / 0.70) ** 2,\n'
     '    "mid  (Krieger-Dougherty, true solvent 0.955 mPa.s, phi_m=0.58)":\n'
     '        0.955e-3 / (1 - 0.44 / 0.58) ** 2,\n'
     '    "high (repo Wyart-Cates fit, eta_s=4.74 mPa.s -- what the DNS uses)":\n'
     '        SUSP.eta_0,\n'
     '}\n'
     'assert list(ETA_BRACKET.values()) == sorted(ETA_BRACKET.values()), \\\n'
     '    f"eta bracket is not ordered low->high: {ETA_BRACKET}"'),
    # geometric mean label
    ('    sigmas = {"lit low (19.1 Pa)": lo, "lit geometric mean (55 Pa)": gm,\n'
     '              "lit high (100 Pa)": hi, "repo fit (44.2 Pa)": SUSP.sigma_star}',
     '    sigmas = {f"lit low ({lo:.3g} Pa)": lo, f"lit geometric mean ({gm:.3g} Pa)": gm,\n'
     '              f"lit high ({hi:.3g} Pa)": hi,\n'
     '              f"repo fit ({SUSP.sigma_star:.3g} Pa)": SUSP.sigma_star}'),
    # drop the uncited "Howland Ec" expression; report the criterion the module owns
    ('    # Howland\'s criterion: substrate matters once the Wagner peak exceeds E\n'
     '    out["_howland_Ec_for_this_impact_Pa"] = float(\n'
     '        DROP.rho * R_LAB * U / (0.17 * (2 * R_LAB) / U))',
     '    # The modulus at which this impact reaches Cauchy number C = 1, i.e. the\n'
     '    # substrate is soft for the whole impact. This is the module\'s own\n'
     '    # criterion (Eq. 15), not a transcription of anyone else\'s worked example.\n'
     '    out["_E_star_at_Cauchy_number_unity_Pa"] = float(imp(U).p_dyn)\n'
     '    out["_note"] = ("Howland et al. PRL 117, 184502 (2016) give E_c = 93 kPa for "\n'
     '                    "THEIR conditions (R = 1.5 mm, U = 2.36 m/s, ethanol); that is "\n'
     '                    "not transferable to these drops and is quoted, not recomputed.")'),
    # radius vs diameter in the ladder row
    ('    row("ref: dilatancy onset sigma* (silica)", SUSP.sigma_star, SUSP.a, np.nan,\n'
     '        "what all of the above must be compared against")',
     '    row("ref: dilatancy onset sigma* (silica)", SUSP.sigma_star, 2.0 * SUSP.a, np.nan,\n'
     '        "length is the particle DIAMETER 2a = 800 nm (SUSP.a is a radius); "\n'
     '        "this is what all of the above must be compared against")'),
    # evaluate the gauges at the justified time, and report t*=1 as a sensitivity
    ('import stress as st', 'import stress as st\nfrom dilatancy import T_STAR_EVAL'),
])

# ------------------------------------------------- inertial_dilatant.py citation
patch("inertial_dilatant.py", [
    ('  [M26] A. Mobaseri, L. Gordillo, C. Burton, S. Yoon, D. Lee, S. Kumar,\n'
     '        M. M. Driscoll & X. Cheng, "Inertia-Dilatancy Interplay Governs\n'
     '        Shear-Thickening Drop Impact", arXiv:2601.12642 (19 Jan 2026).',
     '  [M26] A. Mobaseri, L. Gordillo, C. Burton, S. Yoon, D. Lee, S. Kumar,\n'
     '        M. M. Driscoll & X. Cheng, "Inertia-Dilatancy Interplay Governs\n'
     '        Shear-Thickening Drop Impact", Phys. Rev. Lett. 136, 148201 (2026),\n'
     '        doi:10.1103/fyx7-jb1d; preprint arXiv:2601.12642 (19 Jan 2026).\n'
     '        Equations below are from the arXiv v1; the Supplemental Material is\n'
     '        NOT on arXiv and may carry prefactors this module cannot see.'),
    ('CONSEQUENCE, AND IT IS THE MAIN METHODOLOGICAL POINT OF THIS MODULE:\n'
     '    quote RATIOS of ID between systems, which are prefactor-free.\n'
     '    do NOT quote ID against 1 as if the threshold were calibrated.\n'
     'The ratio silica:cornstarch is the answer the question actually needs.',
     'CONSEQUENCE, AND IT IS THE MAIN METHODOLOGICAL POINT OF THIS MODULE:\n'
     '    quote RATIOS of ID between systems, which are prefactor-free.\n'
     '    do NOT quote ID against 1 as if the threshold were calibrated.\n'
     'The ratio silica:cornstarch is the answer the question actually needs.\n'
     '\n'
     'See `ID_capped` for the proposed resolution: the Darcy pressure inside F_d\n'
     'exceeds what a pore fluid can sustain by 2-3 orders of magnitude, and\n'
     'capping it puts ID = 1.1 at the top of [M26]\'s own velocity range, where\n'
     'their Regime 3 is observed, instead of ID = 0.009.'),
])
print("done")

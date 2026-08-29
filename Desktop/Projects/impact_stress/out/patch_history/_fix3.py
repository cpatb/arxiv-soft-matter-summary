"""Third correction pass: everything the two adversarial referees established."""
import pathlib

root = pathlib.Path(__file__).parent.parent


def patch(name, pairs):
    p = root / name
    t = p.read_text(encoding="utf-8")
    missed = [a[:70] for a, b in pairs if a not in t]
    for a, b in pairs:
        t = t.replace(a, b)
    p.write_text(t, encoding="utf-8")
    print(f"{name}: missed {missed or 'none'}")
    if missed:
        raise SystemExit(1)


# =========================================================== stress.py
patch("stress.py", [
    # --- the shear/pressure ratio was overstated by 10x -------------------
    ("      * PREFACTOR rho U^2 Re_R^(-1/2).  The shear stress is smaller than the\n"
     "        pressure by ~1/sqrt(Re). For the lab's silica suspension Re_R ~ 100,\n"
     "        so the shear stress is ~10x below the pressure — two orders of\n"
     "        magnitude below rho U^2 once the geometric factor is included.\n"
     "        THAT IS THE WHOLE POINT: rho U^2 is the wrong gauge for a\n"
     "        shear-thickening drop, and it is wrong by ~100x.",
     "      * PREFACTOR rho U^2 Re_R^(-1/2).  The shear stress is smaller than the\n"
     "        pressure by ~1/sqrt(Re). At the lab's Re_R = 121, the area-averaged\n"
     "        shear stress is 5.3x below the on-axis pressure and 9.7x below\n"
     "        rho U^2 at t* = 1 (5.5x at t* = 0.315). CORRECTED AFTER REVIEW: an\n"
     "        earlier version of this docstring said '~100x', which is wrong by an\n"
     "        order of magnitude and contradicted dilatancy.py. The right\n"
     "        statement is that rho U^2 overstates the stress that drives\n"
     "        thickening by roughly one decade, which is the difference between a\n"
     "        gauge whose threshold is reachable and one whose is not."),
    # --- delete the false "convergent integral" claim ---------------------
    ('    This is the number to quote. It is (a) a literature result, prefactor\n'
     '    included, (b) insensitive to the screening cut-off because the integral\n'
     '    converges, and (c) exactly the mean tangential traction the substrate\n'
     '    feels, which is what "the stress of the impact" ought to mean.\n    """',
     '    HOW MUCH TO TRUST THE PREFACTOR. Not as much as "5.35" suggests.\n'
     '      * The integral DOES NOT converge. int xi^2/(1-xi^2) dxi is\n'
     '        logarithmically divergent -- PLJ16 say so in as many words -- so\n'
     '        <tau> depends on the screening length, logarithmically. A factor 1.4\n'
     '        in Delta is 11% in F_d; PLJ16 Table 2\'s own factor-2 spread in\n'
     '        Delta is 24%.\n'
     '      * The published 5.35 is internally inconsistent with the Delta that\n'
     '        PLJ16 report. Their stated Delta/sqrt(ttil) ~ 0.03 (Table 2:\n'
     '        0.03-0.06 over Re 250-5000) integrates to a coefficient of\n'
     '        1.02 here, not 1.135; recovering 1.135 needs Delta = 0.0202\n'
     '        sqrt(ttil), below their whole quoted range.\n'
     '      * Integrating F_d out to d-Delta but dividing by the full area pi d^2\n'
     '        (as the published construction does) is a further +3.3%.\n'
     '    So: use C = 1.0 with a band 0.76-1.14, not 1.1353 to five figures. The\n'
     '    published value is kept as the default so the number is traceable, and\n'
     '    `MEAN_SHEAR_COEFF_BAND` carries the honest spread. U_c goes as\n'
     '    C^(-2/3), so this is +8% to +31% on every onset velocity.\n    """'),
    # --- rim amplification: not time-independent, not 5.5 -----------------
    ('    THE HEADLINE\n'
     '    ------------\n'
     '    p_peak / p_axis = (3/(2 pi)) t*^(-1) / ((sqrt(3)/pi) t*^(-1/2))\n'
     '                    = (sqrt(3)/2) t*^(-1/2)\n'
     '\n'
     '    So the rim is hotter than the centre by sqrt(3)/2 / sqrt(t*), which is a\n'
     '    factor of ~3 at t* = 0.1 and ~9 at t* = 0.01. THE STRESS IS IN A RING, AND\n'
     '    THE RING IS THE WHOLE STORY EARLY ON.\n    """',
     '    THE HEADLINE, AND ITS LIMITS\n'
     '    ----------------------------\n'
     '    The stress is in a RING, and the ring is the whole story early on. But\n'
     '    the amplification factor is not a constant of nature:\n'
     '      * It is time-independent ONLY because Delta is taken proportional to\n'
     '        sqrt(t). That proportionality is an empirical fit, not a mechanism;\n'
     '        PLJ16 say the origin of Delta "needs to be clearly identified".\n'
     '      * Their two candidate mechanisms give qualitatively different answers\n'
     '        (a viscous length gives an amplification that GROWS with time and\n'
     '        depends on Re; an inertial Wagner inner layer gives one that is\n'
     '        UNBOUNDED as t -> 0). Both are excluded by PLJ16\'s own Table 2,\n'
     '        which shows Delta/sqrt(ttil) varying less than 2x across Re 250-5000\n'
     '        and ttil 5e-3 to 0.1. So the empirical sqrt(t) survives as an\n'
     '        unexplained fit, over that window only.\n'
     '      * Inside PLJ16\'s own Delta range (0.03-0.06) the pressure ratio is\n'
     '        3.8-5.4 and the shear shape factor is 14-29. Quote 4-6, not 5.5.\n'
     '      * Below t* ~ 2e-5 the rim is capped by compressibility (the axis\n'
     '        crossover is at 6e-7, so the rim hits the water-hammer ceiling ~31x\n'
     '        later in t*), and the ratio falls toward 1.\n'
     '    Correct statement: "p_peak/p_axis = 4-6 for 5e-3 < t* < 0.1, set by an\n'
     '    empirical screening length of unknown origin with a factor-2 uncertainty."\n    """'),
    # --- the omitted 1/2 rho U^2 on the axis ------------------------------
    ('def wagner_pressure_axis(imp: Impact, t_star):\n'
     '    """On-axis pressure p(0,t) = (sqrt3/pi) rho U^2 t*^(-1/2) [Pa].\n'
     '\n'
     '    = 0.5513 rho U^2 t*^(-1/2). CSG22 Eq. [9] gives the same thing in D-clock\n'
     '    variables (plus a constant rho V^2/2 that is a choice of reference\n'
     '    pressure). It equals rho U^2 at t* = 3/pi^2 = 0.304.\n'
     '    """\n'
     '    return ((np.sqrt(3.0) / np.pi) * imp.drop.rho * imp.U**2\n'
     '            * np.asarray(t_star, dtype=float) ** -0.5)',
     'def wagner_pressure_axis(imp: Impact, t_star, include_stagnation: bool = False):\n'
     '    """On-axis pressure p(0,t) [Pa].\n'
     '\n'
     '        singular part:  (sqrt3/pi) rho U^2 t*^(-1/2) = 0.5513 rho U^2 t*^(-1/2)\n'
     '        CSG22 Eq. [9]:  p(0,t) = (1/2) rho U^2 + [the singular part]\n'
     '\n'
     '    CSG22 call the constant (1/2) rho U^2 "a constant dynamic pressure that\n'
     '    can be absorbed by redefining the reference pressure", and Eq. [10] omits\n'
     '    it. That is fine for the SHAPE and for anything deviatoric, and it is\n'
     '    what the default does. It is NOT fine if you want the gauge pressure on\n'
     '    the wall late in the impact: the omitted term is 9% of the total at\n'
     '    t* = 0.01, 29% at t* = 0.1 and 50% at t* = 0.3. Pass\n'
     '    include_stagnation=True there.\n'
     '\n'
     '    The singular part alone equals rho U^2 at t* = 3/pi^2 = 0.304.\n'
     '    """\n'
     '    p = ((np.sqrt(3.0) / np.pi) * imp.drop.rho * imp.U**2\n'
     '         * np.asarray(t_star, dtype=float) ** -0.5)\n'
     '    return p + 0.5 * imp.p_dyn if include_stagnation else p'),
    # --- glass is not "within 0.1%" ---------------------------------------
    ('def rigid_target(name: str = "rigid solid (glass slide)") -> Target:\n'
     '    return Target(name=name, rho=np.inf, Z=np.inf, E=np.inf,\n'
     '                  provenance="idealisation; glass is within ~0.1% of it")',
     'def rigid_target(name: str = "rigid solid (perfectly rigid)") -> Target:\n'
     '    """The idealisation. NOTE that real glass is NOT this in the acoustic\n'
     '    stage: fused silica has Z = 1.3e7 against the suspension\'s 2.3e6, so\n'
     '    Z2/(Z1+Z2) = 0.85 and a glass slide takes 15% off the water-hammer\n'
     '    pressure. Use `glass_target()` for that. In the incompressible stage\n'
     '    glass really is rigid (Cauchy number 8e-8)."""\n'
     '    return Target(name=name, rho=np.inf, Z=np.inf, E=np.inf,\n'
     '                  provenance="idealisation, exact only in the Wagner stage")\n'
     '\n'
     '\n'
     'def glass_target(name: str = "fused-silica microscope slide") -> Target:\n'
     '    """A real glass slide: finite acoustic impedance, effectively infinite\n'
     '    stiffness. rho = 2200 kg/m^3, c_long = 5900 m/s, E = 70 GPa."""\n'
     '    return Target(name=name, rho=np.inf, Z=2200.0 * 5900.0, E=70e9 / (1 - 0.17**2),\n'
     '                  provenance="fused silica; rho c = 1.30e7 Pa s/m, E = 70 GPa")'),
    # --- the compliance threshold carries a 2/pi --------------------------
    ("            rho U^2 / E*  >~  t*                                          (14)\n"
     "\n"
     "        so the substrate is \"soft\" for t* below the COMPLIANCE NUMBER\n"
     "\n"
     "            C = rho U^2 / E*                                              (15)",
     "            t*  <~  (2/pi) * rho U^2 / E*                                 (14)\n"
     "\n"
     "        so the substrate is \"soft\" for t* below (2/pi) times the CAUCHY\n"
     "        NUMBER (Yokoyama et al. Nat. Commun. 17, 607 (2026) name it that)\n"
     "\n"
     "            C = rho U^2 / E*                                              (15)\n"
     "\n"
     "        The 2/pi comes out of the algebra and must not be dropped: with\n"
     "        p = (3/pi) rho U^2 R / d and d = R sqrt(3 t*), the deflection\n"
     "        w = p d / E* = 3 rho U^2 R/(pi E*) is TIME-INDEPENDENT, while the\n"
     "        gap being closed is h_gap = (3/2) R t*. Setting w = h_gap gives\n"
     "        t* = (2/pi) C.\n"
     "\n"
     "        HOWLAND'S CRITERION IS THE SAME STATEMENT. Howland et al. PRL 117,\n"
     "        184502 (2016) write p_peak >~ E with p_peak = rho R U/(t_ej - t_i),\n"
     "        which is algebraically rho U^2/t*_ej. So 'C >~ t*_ej' IS their\n"
     "        criterion, and their worked example (R = 1.5 mm, U = 2.36 m/s,\n"
     "        t_ej - t_i = 30 us, rho = 789) reproduces their E_c = 93 kPa exactly.\n"
     "        They use E, not E* = E/(1-nu^2); at nu = 0.5 those differ by 1.33.\n"
     "\n"
     "        C IS NECESSARY BUT NOT SUFFICIENT. Yokoyama et al. show the peak\n"
     "        force does not collapse on Ca alone: the inertial-to-Hertzian\n"
     "        crossover needs a second group theta = eta/sqrt(rho E R^2) (the\n"
     "        DROP's viscosity), through Z = theta^0.28 Ca, and they classify the\n"
     "        problem as self-similarity of the second kind, so 0.28 is not\n"
     "        obtainable by dimensional analysis."),
    # --- the "3 becomes some other O(1)" caveat is wrong ------------------
    ("        CAVEAT, stated so it is not forgotten: the contact radius prefactor\n"
     "        (the 3 in d^2 = 3RUt) comes from free-surface pile-up on the DROP\n"
     "        side. With a deformable target there is pile-up on both sides, so that\n"
     "        3 becomes some other O(1) number. Eq. (13) is the leading-order\n"
     "        amplitude change; the geometric prefactor is an O(1) uncertainty on\n"
     "        top of it.",
     "        THE 3 DOES NOT CHANGE, and an earlier version of this note wrongly\n"
     "        said it did. In the linearised Wagner problem each side's free\n"
     "        surface responds with the same shape function at the same d, so the\n"
     "        Wagner condition involves only the RELATIVE approach u1 + u2 = U and\n"
     "        is identical to the one-sided problem. Cimpeanu & Moore, JFM 856,\n"
     "        764 (2018) Eq. (3.17) confirm it: a drop on a rigid wall and a drop\n"
     "        on a deep flat pool give an identical d(t). For a CURVED target it\n"
     "        becomes a reduced radius, 1/R_eff = 1/R1 + 1/R2.\n"
     "\n"
     "        The matching condition is continuity of PRESSURE IMPULSE, rho*phi,\n"
     "        not of the potential. Both sides solve the same flat-disc mixed\n"
     "        boundary-value problem and therefore carry the same universal 2/pi,\n"
     "        which is what makes rho1 u1 = rho2 u2 follow rather than be assumed.\n"
     "        Cimpeanu & Moore Eq. (3.13) give the interface coefficient (1+V)/2,\n"
     "        equal to 1 at a rigid wall and 1/2 at a deep pool: a published\n"
     "        confirmation of the factor 1/2 at rho1 = rho2. The generalisation to\n"
     "        rho1 != rho2 is not published there and is derived here.\n"
     "\n"
     "        DO NOT CONFLATE THIS WITH Ve = V/2. Langley et al. and Bisighini et\n"
     "        al. report the interface receding at V/2, which is the same u1 -- but\n"
     "        substituting V -> V/2 into the Wagner pressure would give 2^(-3/2) =\n"
     "        0.354, or 0.25 on a rho V^2 gauge. The answer is 1/2 because d(t) is\n"
     "        set by the RELATIVE approach U while the amplitude is set by u1, so\n"
     "        only one velocity factor is halved.\n"
     "\n"
     "        AND THE SHEAR STRESS DOES SOMETHING ELSE ENTIRELY -- see\n"
     "        `shear_sharing_factor`. It does not halve; for a bath of the same\n"
     "        liquid it VANISHES."),
])

# add the shear-sharing function and the honest coefficient band
p = root / "stress.py"
t = p.read_text(encoding="utf-8")
t = t.replace(
    'RE_CORRECTION_COEFF = 8.0 * np.sqrt(6.0) / (3.0 * np.pi**1.5)',
    'MEAN_SHEAR_COEFF_BAND = (0.76, 1.14)\n'
    '"""Honest spread on MEAN_SHEAR_COEFF, from PLJ16 Table 2\'s own factor-2\n'
    'range in the screening length Delta (0.03-0.06 sqrt(ttil)) plus the +3.3%\n'
    'area-normalisation inconsistency in the published construction. The default\n'
    '1.135 is the literature value; it sits at the top of this band."""\n\n'
    'RE_CORRECTION_COEFF = 8.0 * np.sqrt(6.0) / (3.0 * np.pi**1.5)')

t += '''

def shear_sharing_factor(imp: Impact) -> dict:
    """How the WALL SHEAR STRESS changes on a non-rigid target.

    THIS IS NOT THE PRESSURE LAW, AND THE DIFFERENCE IS LARGE.
    The wall shear comes from the SLIP between the liquid and the surface. On a
    rigid wall the slip is the full outer radial velocity. On a liquid target
    both sides move, and the tangential velocities are

        u_{r,i} = (2/pi) u_i r / sqrt(d^2 - r^2)

    with u_1 = U rho2/(rho1+rho2) and u_2 = U rho1/(rho1+rho2) as in
    `compliance_factors`. The slip that drives the boundary layer is the
    DIFFERENCE:

        |u_1 - u_2| / U = |rho2 - rho1| / (rho1 + rho2)                  (18)

    CONSEQUENCE, and it is the sharpest prediction in this module:

        a drop impacting a deep bath OF ITSELF feels NO wall shear at all.

    At matched density the interface is a symmetry plane with a slip boundary
    condition (CSG22's Fig. 7a caption says exactly that), so the Stokes layer
    never forms. The pressure halves; the shear goes to zero. Anything that
    depends on wall shear -- including the dilatancy gauge this project uses --
    is switched off, not merely reduced.

    Returned factors are the SLIP ratio; the actual tau ratio carries a further
    two-sided boundary-layer factor <= 1 that is not computed here, so treat
    these as upper bounds on tau/tau_rigid.
    """
    rho1, rho2 = imp.drop.rho, imp.target.rho
    if np.isinf(rho2):
        slip = 1.0
    else:
        slip = abs(rho2 - rho1) / (rho1 + rho2)
    return {
        "slip_ratio": slip,
        "tau_over_tau_rigid_upper_bound": slip,
        "note": ("exactly zero for a bath of the same liquid; the interface is a "
                 "symmetry plane and no Stokes layer forms"),
    }
'''
p.write_text(t, encoding="utf-8")
print("stress.py: shear_sharing_factor + MEAN_SHEAR_COEFF_BAND added")

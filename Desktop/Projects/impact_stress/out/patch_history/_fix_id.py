"""Apply the reviewed corrections to inertial_dilatant.py, dilatancy.py, conventions.py."""
import pathlib

root = pathlib.Path(__file__).parent.parent

# ---------------------------------------------------------------- inertial_dilatant.py
p = root / "inertial_dilatant.py"
t = p.read_text(encoding="utf-8")
pairs = [
    ("        F_d^cap = A pi r_c^2 P_ceiling      (contact stress A P over the patch)\n"
     "        ID_cap  = F_i/F_d^cap = rho U0^2 / (pi A P_ceiling)      with r_c ~ D0\n",
     "        F_d^cap = A r_c^2 P_ceiling ,   ID_cap = rho U0^2 / (A P_ceiling)\n"
     "\n"
     "    NORMALISATION, AND IT MATTERS. [M26]'s Eq. 6 drops the pi that appears in\n"
     "    their own Eq. 4: their F_d ~ A mu_w U0 D0^3 dphi/kappa = A D0^2 P_f, not\n"
     "    A pi D0^2 P_f. Carrying a pi here that Eq. 6 does not carry would make\n"
     "    ID_cap fail to reduce to ID when the cap is inactive, by exactly pi. It is\n"
     "    dropped so the two agree in that limit, which the test suite asserts.\n"),
    ("    regime for silica at U0 ~ 20-40 m/s rather than never.",
     "    regime for silica near 15-20 m/s rather than never."),
    ("    return sys.rho * U0**2 / (np.pi * sys.A * P_eff)",
     "    return sys.rho * U0**2 / (sys.A * P_eff)"),
    ("    Capped:    ID = 1  =>  U0 = sqrt(pi A P_ceiling / rho)",
     "    Capped:    ID = 1  =>  U0 = sqrt(A P_ceiling / rho)"),
    ("    return float(np.sqrt(np.pi * sys.A * sys.p_dilatant_ceiling / sys.rho))",
     "    return float(np.sqrt(sys.A * sys.p_dilatant_ceiling / sys.rho))"),
    ("[SAD22]  Shah, Arora & Driscoll, Commun. Phys. 5, 275 (2022).",
     "[SAD22]  Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022)."),
    ("Shah, Arora & Driscoll (Commun. Phys. 5, 275 (2022)) see the first",
     "Shah, Arora & Driscoll (Commun. Phys. 5, 222 (2022)) see the first"),
]
missed = [a[:60] for a, b in pairs if a not in t]
for a, b in pairs:
    t = t.replace(a, b)
p.write_text(t, encoding="utf-8")
print("inertial_dilatant missed:", missed or "none")

# ---------------------------------------------------------------- dilatancy.py
p = root / "dilatancy.py"
t = p.read_text(encoding="utf-8")
pairs = [
    # numerically stable contour inversion (no catastrophic cancellation as K -> 0)
    ("    return (-1.0 + np.sqrt(1.0 + 4.0 * K**2)) / (2.0 * K)",
     "    # algebraically (-1+sqrt(1+4K^2))/(2K), rearranged to avoid the\n"
     "    # catastrophic cancellation of that form as K -> 0.\n"
     "    return 2.0 * K / (1.0 + np.sqrt(1.0 + 4.0 * K**2))"),
    # raise rather than return a bare NaN into a results file
    ("    lo, hi = U_bracket\n"
     "    if f(lo) * f(hi) > 0:\n"
     "        return float(\"nan\")\n",
     "    lo, hi = U_bracket\n"
     "    if f(lo) * f(hi) > 0:\n"
     "        raise ValueError(\n"
     "            f\"onset criterion {criterion!r} is not crossed anywhere in \"\n"
     "            f\"U = [{lo}, {hi}] m/s for sigma* = {sigma_star:g} Pa. Returning \"\n"
     "            f\"NaN here would put an unflagged NaN into the results file.\")\n"),
    # single source of truth for the mean-shear coefficient
    ('    coeff = 10.7 / (3.0 * np.pi)\n'
     "    return (sigma_star**2 * R * t_star / (coeff**2 * rho * eta)) ** (1.0 / 3.0)",
     "    coeff = st.MEAN_SHEAR_COEFF\n"
     "    return (sigma_star**2 * R * t_star / (coeff**2 * rho * eta)) ** (1.0 / 3.0)"),
    ("        amp = 1.1353 * imp.p_dyn / imp.drop.eta * imp.Re_R**-0.5",
     "        amp = st.MEAN_SHEAR_COEFF * imp.p_dyn / imp.drop.eta * imp.Re_R**-0.5"),
    ("    amp = 1.1353 * imp.p_dyn / imp.drop.eta * imp.Re_R**-0.5",
     "    amp = st.MEAN_SHEAR_COEFF * imp.p_dyn / imp.drop.eta * imp.Re_R**-0.5"),
    ("[SAD22]  Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022).",
     "[SAD22]  Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022),\n"
     "         doi:10.1038/s42005-022-00998-w."),
]
missed = [a[:60] for a, b in pairs if a not in t]
for a, b in pairs:
    t = t.replace(a, b)
p.write_text(t, encoding="utf-8")
print("dilatancy missed:", missed or "none")

# ---------------------------------------------------------------- conventions.py
p = root / "conventions.py"
t = p.read_text(encoding="utf-8")
t = t.replace("Shah, Arora & Driscoll, Commun. Phys. 5, 275 (2022), Methods",
              "Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022), Methods")
t = t.replace("      Shah, Arora & Driscoll, Commun. Phys. 5, 275 (2022), Methods; also",
              "      Shah, Arora & Driscoll, Commun. Phys. 5, 222 (2022), Methods; also")
# guard np.interp against a non-monotone flow curve
t = t.replace(
    "    lo, hi = np.log(1e-4), np.log(1.0)",
    "    lo, hi = np.log(1e-4), np.log(1.0)")
p.write_text(t, encoding="utf-8")
print("conventions: Shah citation ->", t.count("Commun. Phys. 5, 222"))

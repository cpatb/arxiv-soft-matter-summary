import numpy as np

# ---- Method: Beer-Lambert attenuation of 365 nm through a coumarin-loaded PEG resin ----
# I(z) = I0 * 10^(-eps*c*z) = I0 * exp(-alpha*z),  alpha = ln(10)*eps*c
# 1/e penetration depth delta = 1/alpha = 1/(ln(10)*eps*c)
# Dominant absorber at 365 nm = CoumAc monomer (7-alkoxy-4-methylcoumarin).
# Coumarin is a BLEACHING absorber: [2+2] dimerization destroys the C3=C4 band -> dimer transparent at 365.

# --- Resin composition: 1:9:90 = mol% PEGDA : CoumAc : PEGMA (paper Table 1) ---
x_PEGDA, x_Coum, x_PEGMA = 0.01, 0.09, 0.90
M_PEGDA, M_Coum, M_PEGMA = 700.0, 274.27, 400.0   # g/mol (paper; CoumAc C15H14O5)
rho_resin = 1.10   # g/cm^3, typical PEG-acrylate resin (assumption, +-10%)

Mbar = x_PEGDA*M_PEGDA + x_Coum*M_Coum + x_PEGMA*M_PEGMA   # mole-avg monomer MW
c_total = rho_resin*1000 / Mbar        # mol/L of monomer units (rho in g/L)
c_coum = x_Coum * c_total              # mol/L CoumAc
print(f"Mean monomer MW      = {Mbar:.1f} g/mol")
print(f"Total monomer conc   = {c_total:.3f} mol/L")
print(f"CoumAc concentration = {c_coum:.3f} mol/L\n")

L = 0.20  # cm (2 mm sample thickness)

print(f"{'eps365':>8} {'alpha(cm-1)':>12} {'delta(um)':>10} {'A(2mm)':>8} {'T_back':>10}")
for eps in [5, 20, 50, 100, 300, 1000, 3000]:
    alpha = np.log(10)*eps*c_coum          # cm^-1
    delta = 1/alpha * 1e4                    # um
    A = eps*c_coum*L                         # absorbance through 2 mm (base 10)
    T = 10**(-A)                             # transmitted fraction at back face
    print(f"{eps:>8} {alpha:>12.1f} {delta:>10.1f} {A:>8.2f} {T:>10.2e}")

# --- Bleaching-front travel time (order-of-magnitude) ---
# Once surface coumarins dimerize they bleach; front advances as each layer receives
# its "clearing dose" D* = photons needed to dimerize a layer.
# Photons per unit volume to dimerize all coumarin: N = c_coum*NA / Phi   (per cm^3, /quantum yield)
# Absorbed-dose front speed v ~ (I0 * lambda/hc * Phi) / (c_coum*NA)  [layer clears when absorbed dose = c/Phi]
print("\n--- Bleaching-front clearing-time scale (needs quantum yield & eps only via absorbed fraction) ---")
h = 6.626e-34; c_light = 3e8; NA = 6.022e23
lam = 365e-9
E_photon = h*c_light/lam                      # J per photon
I0 = 0.1e4                                     # W/cm^2 -> W/m^2? keep W/cm^2 handling below
I0_Wcm2 = 0.1                                  # W/cm^2 incident (paper & user's source)
photon_flux = I0_Wcm2 / E_photon              # photons/s/cm^2 incident
c_coum_percm3 = c_coum*1e-3*NA                 # molecules/cm^3  (mol/L * 1e-3 L/cm^3 * NA)
print(f"Incident photon flux = {photon_flux:.3e} photons/cm2/s")
print(f"CoumAc number density= {c_coum_percm3:.3e} molecules/cm3")
for Phi in [0.05, 0.1, 0.3]:
    # If ALL incident photons were used to dimerize a slab of thickness L (upper bound on speed):
    # time to dimerize full 2mm column = (molecules in column)/(photons absorbed per s * Phi)
    molecules_col = c_coum_percm3 * L          # per cm^2 footprint
    t_ideal = molecules_col/(photon_flux*Phi)  # s, if every photon absorbed & used
    print(f"Phi={Phi:>4}:  ideal min dimerization time for full 2mm column = {t_ideal/60:.1f} min "
          f"(= {t_ideal/3600:.2f} h) -- photon-supply lower bound")

# Pore relaxation and Reynolds dilatancy: 800 nm silica against cornstarch

**Written 2026-08-10. Code: `pore_relaxation.py`, `run_pore_relaxation.py`.
Tests: `test_pore_relaxation.py`, 33 passing. Numbers: `out/pore_relaxation.json`.
Physics and measurement routes: `PORE_PHYSICS.pdf`.**

The question was whether the pore relaxation and the dilatancy approximation can be worked
out for the lab's 800 nm monodisperse silica in water at high volume fraction, and compared
against cornstarch. Both have closed forms. Nothing here needed a solver; the numerical
part is only a sweep over φ.

This extends the 2026-08-06 stress brief rather than replacing it. That document already
carried the permeability, the inertial–dilatant number ID, and the Reynolds–Darcy suction.
What it did not have is the relaxation time itself, for a reason worth stating first.

---

## 1. What the earlier estimate was actually measuring

`inertial_dilatant.poroelastic_time_ratio` puts the whole pore-pressure storage into
water's own compressibility, χ_w = 4.55 × 10⁻¹⁰ Pa⁻¹. That is the rigid-skeleton limit:
grains rigid, network rigid, and the only thing that can absorb a pore-pressure increment
is the water squeezing. It answers the question Jérôme, Vandenberghe & Forterre pose,
which is whether the pore fluid may be treated as incompressible in the Darcy–Reynolds
equation, and its answers of 0.014 for cornstarch and 0.81 for silica are correct for that
question.

It is not the pore relaxation time. A saturated packing sheds an imposed pore pressure by
draining, and what resists the drainage is the stiffness of the particle network, not the
bulk modulus of water. The two storage terms add,

$$S = (1-\phi)\chi_f + \frac{1}{K_p},$$

and for every system below the second term is larger by five to nine decades. Using the
first alone understates the relaxation time by the same factor. Both terms are carried now,
and the test suite pins the new code to reproduce the old result exactly in the K_p → ∞
limit, so the two coexist without contradicting each other.

---

## 2. The analytic core

Mixture continuity plus Darcy for a saturated deformable packing gives a diffusion equation
for the pore pressure (Biot 1941; Terzaghi consolidation):

$$S\,\frac{\partial p}{\partial t} = \nabla\!\cdot\!\left[\frac{\kappa}{\mu_f}\nabla p\right] + \text{(dilation source)},
\qquad D_p = \frac{\kappa}{\mu_f S},
\qquad \tau_p = \frac{L^2}{D_p}.$$

With the network dominating the storage this is just

$$\boxed{\;\tau_p = \frac{\mu_f L^2}{\kappa K_p}\;}
\qquad
\kappa = \frac{(1-\phi)^3 a^2}{180\,\phi^2}\;\;\text{(Carman–Kozeny, }a=\text{diameter)}$$

and D_p = κK_p/μ_f is also, identically, the colloid's collective diffusion coefficient.

The dilatancy approximation is the same closure the DNS already runs behind `DIL_COEFF`
and the same as Mobaseri *et al.* Eq. 3. A packing denser than φ_c must expand to shear,
the expansion draws solvent in from a distance L_d, and Darcy resistance to that inflow
generates the suction:

$$p_{\rm dil} = \frac{\mu_f}{\kappa}\,\Delta\phi\,L_d^2\,\dot\gamma. \tag{D1}$$

Dividing (D1) by K_p ties the two halves together:

$$p_{\rm dil} = K_p\,\Delta\phi\,\bigl(\dot\gamma\,\tau_p\bigr). \tag{D2}$$

The dilatant pressure is the network modulus times the dilation strain accumulated in one
pore-relaxation time. Written that way it has a ceiling, because the strain cannot exceed
the dilation the packing has left in it:

$$p_{\rm dil}^{\max} = \frac{K_p\,\Delta\phi}{\phi}. \tag{D3}$$

The strain is Δφ/φ rather than Δφ/(1−φ), because with incompressible grains the particle
volume φV is fixed, so ε_v = −Δφ/φ, and the same φ appears in the definition
K_p = φ dP/dφ. The first version of this brief used the porosity and was wrong. The two
agree exactly at φ = 0.5, which is why the silica operating point did not catch it. Details
in `PORE_PHYSICS.pdf` §4.3.

(D3) is new here, and it settles Problem 2 of the 06 AUG brief, which was the observation
that uncapped Darcy predicts 10⁷ to 10⁸ Pa of suction for 800 nm particles, a state no pore
fluid can reach. The cavitation cap proposed then was withdrawn because the corrected
pore-throat radius made it non-binding and because it lost all particle-size dependence. A
stiffness ceiling keeps the size dependence, since K_p carries it, and it can be checked
against data rather than against another scaling argument.

**The validation.** Jérôme *et al.* measured −10 to −100 kPa of pore suction in 170 µm
glass beads. That is the one system in this comparison where the suction was measured
rather than extrapolated, so (D3) has to land on it.

The first attempt did not. With a pressure-free contact stiffness k_n ~ E·a the ceiling
came out at 1.1 × 10⁸ Pa, three decades high, and the test failed. The cause is that
k_n ~ E·a describes a contact squashed by δ ~ a, where a real pack at 10 kPa confinement
has δ/a ≈ 5 × 10⁻⁵. Replacing it with a Hertzian stiffness evaluated at the stated
confining pressure gives 4.2 × 10⁵ Pa, within a decade of the measurement and on the
expected side of it, since the affine modulus is an upper bound and real packings measure
two to five times softer near isostaticity. The pressure-free form is kept in the module,
labelled, with a test that holds it to its three-decade error so nobody reinstates it.

For colloidal silica the contacts never touch, so Hertz does not apply. An
electric-double-layer repulsion decays as F = F\* exp(−h/λ_D), giving k_n = F/λ_D, a
stiffness set by the load already carried rather than by any particle modulus. Feeding that
into the affine packing modulus, the coordination number and φ cancel exactly:

$$K_p = \frac{P\,a}{6\lambda_D}.$$

---

## 3. The numbers

φ_c = 0.42 for the silica throughout. That is the value implied by the lab's own
observation of fronts at φ = 0.42, and it is the single largest uncertainty here, because
it is unmeasured. Drainage length is the contact-patch radius at peak force, 0.49 D₀. Shear
rate is U₀/D₀. Each system sits at its own operating pressure: σ\* = 44 Pa for silica,
10⁵ Pa for cornstarch, which is the jammed-column stress Waitukaitis & Jaeger measured.

| | silica, φ=0.44 | silica, φ=0.50 | cornstarch, φ=0.42 | glass beads (JVF) |
|---|---:|---:|---:|---:|
| a (diameter) | 800 nm | 800 nm | 10 µm | 170 µm |
| Δφ | 0.02 | 0.08 | 0.05 | 0.005 |
| κ (m²) | 3.2 × 10⁻¹⁵ | 1.8 × 10⁻¹⁵ | 6.1 × 10⁻¹³ | 3.5 × 10⁻¹¹ |
| K_p (Pa) | 590 | 590 | 1.1 × 10⁶ | 4.9 × 10⁷ |
| D_p (m²/s) | 1.9 × 10⁻⁹ | 1.1 × 10⁻⁹ | 6.8 × 10⁻⁴ | 1.7 |
| **τ_p, contact scale** | **790 s** | **1430 s** | **5.7 ms** | **2.2 µs** |
| τ_p, particle scale | 0.34 ms | 0.61 ms | 0.15 µs | 17 ns |
| τ_impact = D₀/U₀ | 1.1 ms | 1.1 ms | 0.78 ms | 0.78 ms |
| **De = τ_p/τ_impact** | **7 × 10⁵** | **1 × 10⁶** | **7.2** | **0.003** |
| p_dil, uncapped Darcy | 8.6 × 10⁶ Pa | 6.2 × 10⁷ Pa | 4.0 × 10⁵ Pa | 690 Pa |
| **p_dil, ceiling (D3)** | **27 Pa** | **94 Pa** | **1.3 × 10⁵ Pa** | 4.2 × 10⁵ Pa |

Three things fall out.

**Silica is undrained, cornstarch is marginal, glass beads drain freely.** The drainage
number spans nine decades across the three systems. Cornstarch sits at De ≈ 7, which is why
velocity and drop size select its behaviour and why a re-entrant liquid regime exists for
it at all. Silica sits five to six decades higher. On impact timescales its pore pressure
does not relax over any length larger than a few particle diameters.

**The dilatant pressure for silica is 10 to 100 Pa, not 10⁷ Pa.** The ceiling binds by six
decades. At φ = 0.50 it is 94 Pa, the same order as the fitted onset σ\* = 44 Pa, which is
a ratio statement rather than a coincidence, since P_conf was set to σ\*. The content is
the ratio, and §5 gives it without that anchor.

**Drainage length has to be quoted with every relaxation time.** At L = a the time is
independent of particle size, because κ ~ a² cancels L² = a². Smaller particles do not
drain more slowly locally. They drain more slowly across a fixed geometric length. Silica's
τ_p at the particle scale is 0.3 to 0.6 ms, comparable to the impact, so pore pressure does
equilibrate between neighbouring pores while the drop is spreading. It is the contact patch
and the drop that cannot drain.

---

## 4. What the contrast is actually made of

At matched Δφ, silica's pore-pressure diffusivity is 4.8 × 10⁵ times smaller than
cornstarch's. That splits as ×257 from permeability and ×1854 from network stiffness, so on
the face of it stiffness is the larger effect and the a² permeability story told in the
06 AUG brief is not the whole picture.

That reading is wrong, and the control is in the JSON. K_p ~ P for the silica and ~P^(1/3)
for cornstarch, so the stiffness ratio is not a material property. Most of it is the
2000-fold difference in the pressures the two experiments operate at. Put both systems at
the same confining pressure and the networks turn out to be comparably stiff:

| P_conf | K_p silica | K_p cornstarch | ratio | D_p ratio | silica De |
|---:|---:|---:|---:|---:|---:|
| 10³ Pa | 1.3 × 10⁴ | 2.4 × 10⁵ | 17.9 | 6170 | 5.8 × 10⁴ |
| 10⁴ Pa | 1.3 × 10⁵ | 5.1 × 10⁵ | 3.8 | 1329 | 5.8 × 10³ |
| 10⁵ Pa | 1.3 × 10⁶ | 1.1 × 10⁶ | 0.83 | 286 | 5.8 × 10² |

At 10⁵ Pa the silica network is marginally the stiffer of the two, and the entire
diffusivity contrast collapses back onto the permeability factor of about 256. So the
honest statement is that the a² permeability contrast is the material difference, while the
stiffness contrast says where each experiment sits rather than what the particles are made
of.

The verdict that survives all of it is the last column. Silica's drainage number stays above
500 across three decades of assumed confining pressure. Undrained is the robust conclusion,
and it does not depend on getting K_p right.

---

## 5. A prediction the salt series can test

Dividing (D3) by the confining pressure removes P entirely and leaves something that
inherits none of σ\*'s uncertainty:

$$\frac{p_{\rm dil}^{\max}}{P_{\rm conf}} = \frac{a}{6\lambda_D}\cdot\frac{\Delta\phi}{\phi}.$$

Below 1, frustrated dilation cannot raise the particle pressure above what already drives
it, and Reynolds dilatancy is not available as a strengthening mechanism. Above 1, it is.
The crossover is closed-form:

$$\phi^* = \frac{\phi_c}{1 - 6\lambda_D/a}.$$

| added salt | λ_D | φ\* | ratio at φ=0.44 | ratio at φ=0.50 |
|---|---:|---:|---:|---:|
| ~0 (residual) | 30 nm | 0.542 | 0.20 | 0.71 |
| 1 mM | 10 nm | 0.454 | 0.61 | 2.1 |
| 100 mM | 1 nm | 0.423 | 6.1 | 21 |

Two consequences.

At φ = 0.44, which is where the drops actually sit, silica is on the wrong side of the
crossover at low salt. Reynolds dilatancy cannot be the arrest mechanism there unless the
salt is high. That agrees with what the DNS found independently: at φ = 0.44 the fitted
Wyart–Cates closure cannot shear-jam at any stress either, and the free-surface mechanism
(Brown & Jaeger frustrated dilatancy) is needed in any case, because a local constitutive
law has no front solution.

φ\* also moves through the accessible φ window as salt is added, from 0.54 at residual
ionic strength down to 0.42 at 100 mM. That gives a curve rather than a single threshold,
and the salt-screening series already scans the axis it lives on. The direction matches the
one measured for σ\*: more salt, shorter screening, stiffer contacts, dilatancy available at
lower φ.

Worth flagging without leaning on it. At 1 mM the prediction is φ\* = 0.454, and Shah, Arora
& Driscoll see the first localised solid patches at φ = 0.47 on the same 830 nm particles.
λ_D = 10 nm was a default rather than a measured value for their samples, so this is a
coincidence worth testing, not evidence.

---

## 6. What is uncertain, in order

1. φ_c is unmeasured for this silica. Everything with Δφ in it scales with it, φ\* moves
   with it roughly one-for-one, and the three available anchors (Shah's 0.47, the fitted
   φ_M = 0.5176, the lab's own fronts at 0.42) cannot all sit on one φ axis. This is the
   same blocker the 06 AUG brief identified and it has not moved.
2. K_p is not measured for either system. It is bracketed everywhere: λ_D from 30 down to
   1 nm for silica, E_p from 10⁷ to 10⁹ Pa for cornstarch. The undrained verdict survives
   the whole bracket. The absolute ceilings do not.
3. The affine modulus is an upper bound. Real packings relax non-affinely and measure two
   to five times softer near isostaticity, so every ceiling here is correspondingly high.
4. Δφ = 0.005 for the glass beads is an estimate rather than a reported value, and it is
   the largest single uncertainty in the one comparison against measured data.
5. Carman–Kozeny at high φ is fine over 0.40 to 0.58, since it returns 10⁻³d² at random
   close packing, which is the standard value. But it is a hydraulic-radius model and knows
   nothing about colloidal ordering.

---

## 7. What would settle it

φ_c, by DLS and ζ plus a settled φ axis. That is already the top item on the open-problems
slate and already gates three other things, and it is the same measurement.

K_p directly, by centrifugation or filtration at two or three volume fractions. Every
absolute number here depends on it, and it is a benchtop measurement rather than a
rheometer one.

The crossover itself. φ\*(salt) is a curve the existing salt series can be replotted
against at the cost of no new experiments, provided the φ axis is trusted.

---

## Figure notes

All three ordinates are logarithmic. Every quantity plotted spans four to nine decades, so
a zero-anchored linear axis would compress the entire result into the top pixel row. The
abscissas are linear and show their full physical range.

`figures/fig6_pore_relaxation.png` shows τ_p (left) and De (right) against φ at the
contact-patch drainage length. Blue is silica, orange is cornstarch. The dotted horizontals
on the left are each system's own impact time; the solid line on the right is De = 1. The
gap between the two curves is the result: five to six decades, roughly constant in φ.

`figures/fig7_dilatant_amplification.png` shows p_dil^max/P_conf against φ for the three
salt concentrations, silica only. The horizontal line is the crossover at 1. The two dotted
verticals are φ = 0.44 (the lab's drops) and φ = 0.47 (Shah's first solid patches). Where a
curve sits below the horizontal, Reynolds dilatancy cannot amplify.

`figures/fig8_contrast_decomposition.png` shows the cornstarch-to-silica ratio of
permeability, network modulus, and their product against matched Δφ. Read it with §4: the
middle curve is mostly the pressure difference between the two experiments, not a material
contrast.

---

## Reproduce

```
python -m pytest test_pore_relaxation.py -q     # 33 tests
python run_pore_relaxation.py                   # JSON + figures 6-8
```

Versions as run: numpy 2.4.2, scipy 1.17.1, matplotlib 3.10.8.

## Sources

Biot, *J. Appl. Phys.* **12**, 155 (1941) ·
Pailha & Pouliquen, *JFM* **633**, 115 (2009) ·
Jérôme, Vandenberghe & Forterre, *PRL* **117**, 098003 (2016) ·
Mobaseri *et al.* (incl. Driscoll & Cheng), *PRL* **136**, 148201 (2026) ·
Athani, Metzger, Forterre & Mari, *PRF* **10**, 043301 (2025) ·
Makse, Gland, Johnson & Schwartz, *PRL* **83**, 5070 (1999) ·
Agnolin & Roux, *PRE* **76**, 061302 (2007) ·
Waitukaitis & Jaeger, *Nature* **487**, 205 (2012) ·
Shah, Arora & Driscoll, *Commun. Phys.* **5**, 222 (2022) ·
Carnahan & Starling, *J. Chem. Phys.* **51**, 635 (1969).

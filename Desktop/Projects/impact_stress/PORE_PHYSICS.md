---
title: "Pore draining and pore relaxation in a dense suspension"
subtitle: "What the quantities mean, how to calculate them, and how to measure them, for 800 nm silica and for cornstarch"
author: "Working brief, Driscoll lab"
date: "2026-08-10"
---

# 0. What this document is for

The lab needs two numbers for its 800 nm silica at high volume fraction: how long the pore
pressure takes to relax, and how large a pressure Reynolds dilatancy can generate. Both
have closed forms. Neither can be evaluated without inputs the project does not yet
measure.

This document works the physics from the two-phase mass balances up, so that every
assumption is visible and every input can be traced to something measurable. It then gives
the measurement routes in order of cost, including one that skips the hardest input
entirely.

Companion documents:

| File | What it holds |
|---|---|
| `PORE_RELAXATION_BRIEF.md` | The answer and the comparison, four pages |
| `STRESS_AND_DILATANCY_BRIEF.pdf` | The 2026-08-06 stress ladder and the ID number |
| `pore_relaxation.py` | Every formula here, executable |
| `test_pore_relaxation.py` | 33 tests, including the validation against measured data |
| `out/pore_relaxation.json` | Every number here, machine-readable |

Notation, fixed once. $a$ is a particle **diameter**. It is a radius in the Wyart--Cates
literature and a diameter in Carman--Kozeny, and mixing them costs a factor of four in the
permeability. $\phi$ is the solid volume fraction, $\phi_c$ the dilatancy onset,
$\Delta\phi = \phi - \phi_c$. Pressures are in Pa, lengths in m, and $\mu_f$ is the pore
fluid's viscosity. A full symbol table is in §12.

\newpage

# 1. Three different things get called "pore relaxation"

Before deriving anything it is worth separating three times that all have units of seconds
and all get called the pore relaxation time. They differ by up to nine decades for the same
material, so picking the wrong one is not a rounding error.

**(i) The acoustic time, $\tau_{\rm ac} \sim L/c$.** How long a pressure signal takes to
cross the sample at the speed of sound in the mixture. For a 2.5 mm drop with
$c \approx 1630$ m/s this is 1.5 µs. It governs the water-hammer stage of impact and
nothing else here.

**(ii) The fluid-compressibility time,
$\tau_{\chi} = \mu_f\chi_f L^2/\kappa$.** How long it takes to drain a pore-pressure
increment when the *only* thing storing that increment is the pore water squeezing. This is
the rigid-skeleton limit: rigid grains, rigid network, compressible fluid. It is the time
that appears in Jérôme, Vandenberghe and Forterre's validity condition for treating the
pore fluid as incompressible, and it is what
`inertial_dilatant.poroelastic_time_ratio` computes.

**(iii) The poroelastic consolidation time, $\tau_p = \mu_f L^2/(\kappa K_p)$.** How long it
takes to drain a pore-pressure increment when the storage is provided by the *particle
network* deforming. This is the pore relaxation time in the sense that matters for
dilatancy, and it is the one this document is about.

The three differ because they use different storage terms. Water's compressibility gives an
effective stiffness of $1/\chi_f = 2.2$ GPa. The particle network of a colloidal suspension
at $\phi = 0.5$ has a bulk modulus of order $10^2$--$10^3$ Pa. The ratio is six to seven
decades, and it is the whole reason (ii) and (iii) disagree so violently.

Both (ii) and (iii) are correct answers to their own questions. Neither substitutes for the
other. §2 derives a single equation containing both.

\newpage

# 2. Deriving the consolidation equation

## 2.1 The two mass balances

Treat the suspension as two interpenetrating continua. The solid phase occupies volume
fraction $\phi$ and moves at $\mathbf v_s$; the fluid occupies $1-\phi$ and moves at
$\mathbf v_f$. With incompressible grains,

$$\frac{\partial\phi}{\partial t} + \nabla\!\cdot\!(\phi\,\mathbf v_s) = 0, \tag{2.1}$$

$$\frac{\partial}{\partial t}\bigl[(1-\phi)\rho_f\bigr] + \nabla\!\cdot\!\bigl[(1-\phi)\rho_f\mathbf v_f\bigr] = 0. \tag{2.2}$$

Let the fluid be slightly compressible, $\rho_f = \rho_{f0}(1 + \chi_f p)$, with $\chi_f$
the isothermal compressibility and $p$ the pore pressure. Expanding (2.2) and keeping terms
to first order in $\chi_f p$:

$$\frac{\partial(1-\phi)}{\partial t} + \nabla\!\cdot\!\bigl[(1-\phi)\mathbf v_f\bigr] + (1-\phi)\chi_f\frac{\partial p}{\partial t} = 0. \tag{2.3}$$

Adding (2.1) and (2.3), the $\partial\phi/\partial t$ terms cancel and we are left with a
statement about the total volume flux:

$$\nabla\!\cdot\!\mathbf v_s + \nabla\!\cdot\!\mathbf q + (1-\phi)\chi_f\frac{\partial p}{\partial t} = 0,
\qquad \mathbf q \equiv (1-\phi)(\mathbf v_f - \mathbf v_s). \tag{2.4}$$

$\mathbf q$ is the **filtration flux**, also called the superficial or Darcy velocity: the
volume of fluid crossing unit total area per unit time, measured in the frame of the moving
solid. It is not the interstitial fluid velocity, and the factor $(1-\phi)$ between them is
a standard place to lose a factor.

## 2.2 Darcy's law

Darcy's law closes $\mathbf q$:

$$\mathbf q = -\frac{\kappa}{\mu_f}\nabla p. \tag{2.5}$$

This is the Stokes-flow drag of the fluid through the pore network, coarse-grained. It
assumes: creeping flow through the pores (pore Reynolds number $\ll 1$), a pore geometry
large compared with any molecular or electrostatic length, no slip at the particle
surfaces, and a representative volume containing many pores. §5.2 shows the third of these
is the one that bites for 800 nm silica at low salt.

The pore Reynolds number for the lab's system is reassuringly small. With a pore velocity
of order $q/(1-\phi)$ and a pore length of order $a$, even a filtration flux of 1 mm/s
gives $\mathrm{Re}_{\rm pore} \sim \rho\,q\,a/\mu_f \sim 10^{-3}$. Creeping flow is safe.

## 2.3 Effective stress

The last ingredient is how the skeleton responds. Terzaghi's principle: the deformation of
a saturated granular skeleton is controlled not by the total stress $\sigma$ but by the
**effective stress** carried by the contacts,

$$\sigma' = \sigma - p. \tag{2.6}$$

This is why the distinction matters physically. Load a saturated packing fast enough that
no fluid can move and the pore pressure rises to carry essentially the whole increment, so
$\sigma'$ barely changes and the particles barely feel it. That is the undrained limit.
Wait long enough for the fluid to drain and $p$ returns to ambient, so the entire increment
transfers onto the contacts. Consolidation is the transfer, and $\tau_p$ is how long it
takes.

Define the drained bulk modulus of the skeleton by

$$K_p \equiv -V\frac{\partial\sigma'}{\partial V}\bigg|_{p} = \phi\,\frac{\partial\sigma'}{\partial\phi}, \tag{2.7}$$

where the second form uses incompressible grains: the particle volume $\phi V$ is fixed, so
$V = V_s/\phi$ and $\mathrm dV/V = -\mathrm d\phi/\phi$. **This is the same $\phi$, not
$1-\phi$, that reappears in the dilatancy ceiling in §4.3, and getting it wrong there is a
mistake this project made and corrected.**

At constant total stress, $\partial\sigma'/\partial t = -\partial p/\partial t$, so the
skeleton's volumetric strain rate is

$$\nabla\!\cdot\!\mathbf v_s = \frac{1}{V}\frac{\partial V}{\partial t} = -\frac{1}{K_p}\frac{\partial\sigma'}{\partial t} = \frac{1}{K_p}\frac{\partial p}{\partial t}. \tag{2.8}$$

## 2.4 The storage coefficient

Substituting (2.5) and (2.8) into (2.4):

$$\boxed{\;S\,\frac{\partial p}{\partial t} = \nabla\!\cdot\!\left[\frac{\kappa}{\mu_f}\nabla p\right],
\qquad S = \underbrace{\frac{1}{K_p}}_{\text{skeleton}} + \underbrace{(1-\phi)\chi_f}_{\text{pore fluid}}\;} \tag{2.9}$$

The two storage terms **add**, which means the *larger compliance wins*. This is worth
dwelling on, because the intuition runs the other way: one might expect the stiffer element
to dominate, as it does for springs in parallel. Here the elements are in series as far as
storing a pressure increment goes, and the softer one absorbs almost all of it.

For the lab's silica at $\phi = 0.50$: $1/K_p = 1.7\times10^{-3}$ Pa$^{-1}$ against
$(1-\phi)\chi_f = 2.3\times10^{-10}$ Pa$^{-1}$. The skeleton carries $99.99999\%$ of the
storage. For cornstarch it carries $99.97\%$. Only for a genuinely stiff pack such as glass
beads does the fluid term reach a percent.

So in every case of interest, (2.9) simplifies to a plain diffusion equation:

$$\frac{\partial p}{\partial t} = D_p\nabla^2 p, \qquad \boxed{\;D_p = \frac{\kappa K_p}{\mu_f}\;} \tag{2.10}$$

$D_p$ is the pore-pressure diffusivity, also called the consolidation coefficient in soil
mechanics and, in colloid science, the **collective (or gradient) diffusion coefficient**.
Those are the same quantity under three names, and §7 turns that fact into the cheapest
available measurement.

## 2.5 The relaxation time, exactly

Solve (2.10) in one dimension for a layer of thickness $L$ drained at $x=0$ and sealed at
$x=L$, starting from a uniform overpressure $p_0$. Separating variables gives a cosine
series whose slowest mode decays as $\exp(-t/\tau_1)$ with

$$\tau_1 = \frac{4L^2}{\pi^2 D_p} = 0.405\,\frac{L^2}{D_p}. \tag{2.11}$$

Throughout the brief and the code we quote the **scaling**

$$\tau_p = \frac{L^2}{D_p} = \frac{\mu_f L^2}{\kappa K_p} \tag{2.12}$$

and report the $4/\pi^2$ factor separately (`slowest_mode_time`) rather than folding it in.
The reason is arithmetic hygiene: every conclusion in this analysis is a statement about
decades, and a formula you can check on the back of an envelope is worth more than one that
is $2.5\times$ more accurate. If a number is ever needed to better than a factor of two,
use (2.11) and state the boundary conditions it assumes.

\newpage

# 3. The drainage length is a modelling choice

$L$ is the distance pore fluid must actually travel to equilibrate. It is not a material
property, and quoting $\tau_p$ without $L$ is meaningless. Three lengths matter for a drop
impact.

**The particle length, $L = a$.** Equilibration between neighbouring pores. Note what
happens to the formula: $\kappa \propto a^2$, so $L^2/\kappa$ is independent of $a$, and

$$\tau_p(L=a) = \frac{180\,\mu_f\,\phi^2}{(1-\phi)^3 K_p} \tag{3.1}$$

**contains no particle size at all.** Small particles do not drain more slowly at the pore
scale. They drain more slowly across a *fixed geometric* length, because it takes more pores
to span it. This is the single most common confusion in this area and it is worth stating
in any talk on the subject.

**The contact length, $L = 0.49\,D_0$.** The contact-patch radius at the force peak, taken
from the Wagner analysis in `stress.py`, where it is the computed $r_c$ at peak force rather
than an assumption. This is the length that governs whether the spreading lamella can
drain, and it is the length implied by Mobaseri *et al.*'s Eq. 3.

**The drop length, $L = D_0$.** Global drainage.

For the lab's silica these give $\tau_p$ = 0.6 ms, 1430 s, and 5960 s at $\phi = 0.50$. The
spread is seven decades, from a single material, purely through $L$.

The right choice depends on the question. For *whether frustrated dilatancy can raise the
particle pressure in the spreading film*, the contact length is correct: the fluid must be
drawn from somewhere in the film to fill the dilating region, and the film's own radius
bounds that distance. For *whether the packing is locally in mechanical equilibrium at the
pore scale*, the particle length is correct.

\newpage

# 4. Reynolds dilatancy

## 4.1 The dilatancy law

A packing denser than its critical state cannot shear without expanding: the particles have
to ride up over one another. Reynolds observed this in 1885. The modern statement, following
Roux and Radjai and used in the two-phase form by Pailha and Pouliquen, is that the solid
fraction relaxes toward a shear-rate-dependent equilibrium value:

$$\frac{1}{\phi}\frac{\mathrm D\phi}{\mathrm Dt} = -\tan\psi\,\dot\gamma,
\qquad \tan\psi = K_{\rm dil}\bigl(\phi - \phi_{\rm eq}(I)\bigr) \approx \alpha\,\Delta\phi. \tag{4.1}$$

$\psi$ is the dilatancy angle. The linear closure $\tan\psi = \alpha\Delta\phi$ with
$\alpha \sim 1$ is what Mobaseri *et al.* use and what this analysis inherits. It is a
one-parameter model of a genuinely more complicated constitutive response, and $\alpha$ is
not measured for either system here.

## 4.2 Coupling to Darcy

Dilation creates pore volume that must be filled by fluid drawn in from outside. Using (2.4)
in the incompressible-fluid limit, $\nabla\!\cdot\!\mathbf q = -\nabla\!\cdot\!\mathbf v_s$,
and from (4.1) the skeleton's volumetric expansion rate is $\tan\psi\,\dot\gamma$. So

$$\nabla\!\cdot\!\left[\frac{\kappa}{\mu_f}\nabla p\right] = \tan\psi\,\dot\gamma. \tag{4.2}$$

This is a Poisson equation for the pore pressure with the dilation as its source. Over a
region of size $L_d$ it scales as $(\kappa/\mu_f)(p/L_d^2) \sim \tan\psi\,\dot\gamma$, giving

$$\boxed{\;p_{\rm dil} = \frac{\mu_f}{\kappa}\,\Delta\phi\,L_d^2\,\dot\gamma\;} \tag{4.3}$$

taking $\alpha = 1$. The sign is a **suction**: the fluid pressure drops, and by Terzaghi
(2.6) the effective stress on the contacts rises by the same amount. That rise is the
mechanism. Frustrated dilation does not push the particles apart, it sucks them together.

Equation (4.3) is Mobaseri *et al.*'s Eq. 3 with $|\mathrm dh/\mathrm dt| = L_d\dot\gamma$,
and it is the same closure the DNS runs behind `DIL_COEFF` in `wc_drop_impact.c`, where it
appears as $p_{\rm dil} = (\eta_f/k)\,\epsilon_{\rm dil}L_d^2\dot\gamma$. Three
independent statements of one relation, which is a good sign.

## 4.3 The identity, and the ceiling

Divide (4.3) by $K_p$ and use $\tau_p = \mu_f L_d^2/(\kappa K_p)$:

$$p_{\rm dil} = K_p\,\Delta\phi\,\bigl(\dot\gamma\,\tau_p\bigr). \tag{4.4}$$

The dilatant pressure is the network modulus times the dilation strain accumulated in one
pore-relaxation time. This form is more useful than (4.3) because it makes the saturation
obvious. $\dot\gamma\tau_p$ is a strain, and the packing does not have unlimited strain
available: once it has expanded from $\phi$ to $\phi_c$ there is no further dilation left to
frustrate. From (2.7) that maximum volumetric strain is $\Delta\phi/\phi$, so

$$\boxed{\;p_{\rm dil}^{\max} = \frac{K_p\,\Delta\phi}{\phi}\;} \tag{4.5}$$

**A correction on the record.** The first version of this analysis wrote (4.5) with
$\Delta\phi/(1-\phi)$, reasoning from the porosity. That is wrong, and the derivation in
§2.3 is why: with incompressible grains the volumetric strain is referenced to the solid
fraction, $\mathrm dV/V = -\mathrm d\phi/\phi$, and the same $\phi$ appears in the
definition of $K_p$ itself. Using the porosity in one and the solid fraction in the other
mixes conventions. The two forms coincide *exactly* at $\phi = 0.5$, which is precisely the
operating point the analysis was first evaluated at, so nothing caught it until the
derivation was written out longhand. The error is $27\%$ at $\phi = 0.44$ and $38\%$ at
$\phi = 0.58$. It changes no conclusion and it moves every printed ceiling.
`test_pore_relaxation.py` now pins the corrected form, including the coincidence at
$\phi = 0.5$ that hid it.

## 4.4 Why the ceiling matters

Uncapped, (4.3) predicts $6\times10^7$ Pa of pore suction for 800 nm silica at
$\phi = 0.50$. That is 600 atmospheres of tension in water. The 2026-08-06 brief flagged
this as an open problem, proposed a cavitation cap, and then withdrew the cap for four
reasons, the strongest being that once the pore-throat radius is computed correctly the
capillary ceiling never binds and the capped expression loses all particle-size dependence.

(4.5) is a better answer because it is a statement about what the *material* can do rather
than about what the *fluid* can survive, and because it keeps a size dependence through
$K_p$. It also has the right structure: it says the Darcy expression (4.3) describes a
*transient build-up rate*, valid until the packing runs out of dilation, at which point the
pressure plateaus. The physical picture is a spring being wound up, not a pressure vessel
being burst.

**It can be checked against data.** Jérôme, Vandenberghe and Forterre measured pore
suctions of $-10$ to $-100$ kPa in 170 µm glass beads under impact. That is the one system
in this comparison where the suction was measured rather than extrapolated. §6.4 and §9.3
show what (4.5) gives there, and what had to be fixed to make it agree.

\newpage

# 5. How to know the permeability

## 5.1 Carman--Kozeny

The standard estimate treats the pore space as a bundle of tortuous channels and applies
Poiseuille flow. The hydraulic radius of the pore space is

$$R_h = \frac{\text{void volume}}{\text{wetted area}} = \frac{1-\phi}{S_v\,\phi},
\qquad S_v = \frac{6}{a}\ \text{for spheres}, \tag{5.1}$$

and folding the Kozeny constant and the tortuosity into a single empirical factor of 180
gives

$$\boxed{\;\kappa = \frac{(1-\phi)^3\,a^2}{180\,\phi^2}\;} \tag{5.2}$$

**Validity.** (5.2) works well for packings of roughly equal spheres over
$0.25 \lesssim \phi \lesssim \phi_{\rm RCP}$. It returns $\kappa = 1.0\times10^{-3}a^2$ at
random close packing, which is the standard measured value, and the test suite checks that.
It degrades for high-porosity fibrous media, for strongly non-spherical particles, and for
broad size distributions, where the Sauter mean diameter should be used instead of the
arithmetic mean. For a monodisperse silica this is the right formula.

**The size dependence is the whole story.** $\kappa \propto a^2$, so the lab's 800 nm silica
is $(10\,\mu\mathrm m/800\,\mathrm{nm})^2 = 156$ times less permeable than cornstarch at
matched $\phi$, and about 256 times at the volume fractions the two experiments actually
use. Everything downstream inherits this factor.

## 5.2 Where it breaks for 800 nm silica: the double layer in the pore throat

Carman--Kozeny treats the pore space as inert plumbing. That requires the electrical double
layer to be thin compared with the pore throat. The throat radius of a dense sphere pack is
the largest sphere that fits between three touching particles,
$r_{\rm throat} = 0.155\,R = 0.0775\,a$, which is **62 nm** for 800 nm silica.

Comparing that with the Debye length across the lab's salt range
(`debye_to_throat_ratio`):

| added salt | $\lambda_D$ | $\lambda_D/r_{\rm throat}$ | continuum Darcy? |
|---|---:|---:|---|
| 100 mM | 1.0 nm | 0.016 | fine |
| 1 mM | 9.6 nm | 0.16 | fine |
| nominal 0 (0.01 mM residual) | 95 nm | **1.5** | **no** |

At the low-salt end the double layers of opposing pore walls overlap across the entire
throat. Two things follow. The counterion cloud is dragged along by the flow and generates
a streaming potential whose back-conduction opposes the flow, which is the **electroviscous
effect** and raises the apparent viscosity. And the geometric pore area available to bulk
flow is reduced.

Both reduce the *effective* permeability below the geometric one, so the low-salt end of the
salt series is **less** permeable than (5.2) says, not more. That strengthens the undrained
conclusion rather than threatening it, which is why no correction is applied here. But it
does mean that a measured $\kappa$ at 0 mM should not be expected to match Carman--Kozeny,
and that the salt axis moves the permeability as well as the contact stiffness. Anyone
fitting $\kappa$ across the salt series should expect two effects, not one.

A caution on the residual ionic strength itself. The 95 nm figure assumes $10^{-5}$ M from
dissolved CO$_2$. A concentrated silica suspension also releases its own counterions, which
at $\phi = 0.5$ can push the ionic strength considerably higher and $\lambda_D$
correspondingly lower. The bracket used in the calculations (30 nm) is a compromise between
these. **The ionic strength of the nominally salt-free samples is not measured, and it
should be**: a conductivity probe reading takes minutes and removes a factor that currently
spans a decade.

## 5.3 Three ways to measure $\kappa$

**Constant-flux filtration (a Darcy cell).** Push fluid through a packed bed at a known
volumetric flux and measure the pressure drop; $\kappa = \mu_f q L_{\rm bed}/\Delta p$. This
is the definitional measurement. For a colloid it is awkward: the required pressures are
high, the bed consolidates during the measurement (so $\phi$ is not the $\phi$ you set), and
you need a membrane fine enough to retain 800 nm particles without dominating the
resistance.

**Hindered settling.** Let the suspension sediment and track the interface. At steady state
the buoyant weight of the particles is balanced by Darcy drag, giving

$$\kappa = \frac{\mu_f\,(1-\phi)\,U_{\rm settle}}{\phi\,\Delta\rho\,g}. \tag{5.3}$$

For the lab's silica at $\phi = 0.50$ this predicts $U_{\rm settle} \approx 1.5$ mm per day.
Slow, but a cathetometer or a photograph a day for a week measures it, and it needs no
equipment the lab does not have. The caveat is that at high $\phi$ the suspension may not
settle as a clean interface, and the network can support itself against gravity entirely if
it has yielded.

**Infer it from the collective diffusivity.** See §7. This is the cheapest route and it
measures the product $\kappa K_p$, which is what $\tau_p$ actually needs.

\newpage

# 6. How to know the network modulus

$K_p$ is the input with the widest uncertainty and the least literature support. This
section gives three physical regimes and four measurement routes.

## 6.1 What it is

From (2.7), $K_p = \phi\,\partial\sigma'/\partial\phi$: the change in effective particle
pressure per unit fractional compression, at fixed pore pressure. In colloid science the
same quantity, measured under quiescent conditions, is the **osmotic bulk modulus**
$K = \phi\,\partial\Pi/\partial\phi$. The two coincide when the particle pressure is
thermodynamic in origin, and diverge once contacts form, because contact forces are not in
equilibrium with a reservoir.

That distinction sets up the three regimes.

## 6.2 The Brownian regime, $\phi < \phi_c$

Below the contact onset the particles are held apart by thermal pressure. For hard spheres
the Carnahan--Starling equation of state gives

$$\Pi = n k_B T\,Z(\phi), \qquad Z = \frac{1+\phi+\phi^2-\phi^3}{(1-\phi)^3},
\qquad n = \frac{6\phi}{\pi a^3}. \tag{6.1}$$

Its worth as an anchor is that it is exact to fourth order in $\phi$ and reproduces the
hard-sphere freezing pressure to a few percent, which the test suite checks.

For 800 nm silica at $\phi = 0.50$ this gives $K = 0.43$ Pa. That number deserves a moment.
It is *tiny*: the entire osmotic stiffness of a concentrated colloidal suspension is well
under a pascal, because $n k_BT$ for 800 nm particles is only $7.6\times10^{-3}$ Pa. This is
why colloidal consolidation takes days, why these suspensions do not settle, and why any
network stiffness above a few pascals must come from contacts rather than from thermal
pressure.

For cornstarch, $K_{\rm osm} = 9\times10^{-5}$ Pa. Thermal pressure is irrelevant. The shear
Péclet number confirms it from the other side: $\mathrm{Pe} = 7\times10^{5}$ for cornstarch
under impact against $\mathrm{Pe} = 274$ for the silica. Both are athermal in flow, but the
silica only by two decades, and its *equilibrium* modulus is Brownian.

## 6.3 The double-layer regime: the silica case

Once particles are pressed into contact the stiffness is set by whatever repulsion holds
them apart. For electrostatically stabilised silica that is the double layer, and a
screened-Coulomb repulsion decays exponentially with surface separation $h$:

$$F(h) = F^\ast e^{-h/\lambda_D} \quad\Longrightarrow\quad k_n = -\frac{\mathrm dF}{\mathrm dh} = \frac{F}{\lambda_D}. \tag{6.2}$$

The stiffness is set by the load the contact already carries, divided by the screening
length. There is no particle modulus in it, correctly, since these surfaces never touch.

The mean contact force in a packing at isotropic particle pressure $P$ follows from summing
$f\ell$ over the $Z n/2$ contacts per unit volume:

$$f = \frac{\pi P a^2}{Z\phi}. \tag{6.3}$$

The affine (Voigt) estimate of the packing's bulk modulus with normal springs is
$K_p = Z\phi k_n/(6\pi a)$. Substituting (6.2) and (6.3), $Z$ and $\phi$ **cancel exactly**:

$$\boxed{\;K_p = \frac{P\,a}{6\lambda_D}\;} \tag{6.4}$$

This is a satisfying result. It has one free parameter, $\lambda_D$, which is set by a
measurable salt concentration. It says the network's stiffness *relative to its own
confining pressure* is fixed by a pure geometric ratio, $a/6\lambda_D \approx 13$ at
$\lambda_D = 10$ nm. And it predicts that **adding salt stiffens the network**, in the same
direction as the measured rise of $\sigma^\ast$ with salt.

Two caveats. The affine estimate is an upper bound; real packings relax non-affinely and
measure two to five times softer near isostaticity. And (6.2) presumes a pure
screened-Coulomb repulsion, whereas Ishibashi *et al.* fit these particles with a
double-layer term *plus* a thermal term and find the double-layer contribution vanishing
above a critical coagulation concentration. Above roughly 100 mM, (6.4) should not be
trusted.

## 6.4 The Hertzian regime: cornstarch and glass beads

For solid grains in true contact the stiffness is Hertzian:

$$E^\ast = \frac{E}{2(1-\nu^2)},\quad R^\ast = \frac{a}{4},\quad
\delta = \left(\frac{3f}{4E^\ast\sqrt{R^\ast}}\right)^{2/3},\quad
k_n = \frac{\mathrm df}{\mathrm d\delta} = 2E^\ast\sqrt{R^\ast\delta}. \tag{6.5}$$

The essential feature is $k_n \propto P^{1/3}$: **a granular network has no
pressure-independent modulus**, and quoting $K_p$ for a granular pack without stating the
confining pressure is meaningless.

This is not a fine point. The first version of this analysis used the pressure-free scale
$k_n \sim E a$, which corresponds to a contact squashed by $\delta \sim a$. For glass beads
at 10 kPa confinement the real overlap is $\delta/a \approx 5\times10^{-5}$, and the
pressure-free estimate is **three decades too stiff**. Fed into (4.5) it predicted
$1.1\times10^8$ Pa of dilatant suction in Jérôme's beads against a measured $10^4$ to
$10^5$ Pa, and the validation test failed. With (6.5) at the stated confinement the
prediction is $4.2\times10^5$ Pa: within a decade, and on the expected side, since the
affine modulus is an upper bound. The pressure-free form is retained in the code with a
test that holds it to its three-decade error, so it cannot be quietly reinstated.

## 6.5 Four ways to measure $K_p$

**Osmotic stress / dialysis.** Equilibrate the suspension against a polymer solution of
known osmotic pressure through a semipermeable membrane, and measure the resulting $\phi$.
Repeat at several pressures to get $\Pi(\phi)$, then differentiate. This is the gold
standard for colloids and it directly returns the quantity in (2.7). It is slow: equilibration
at $\phi \sim 0.5$ takes days to weeks, for exactly the reason this document is about.

**Centrifugation to equilibrium.** Spin to a known $g$, let the sediment reach equilibrium,
and read the $\phi(z)$ profile. The profile *is* $\Pi(\phi)$, since at equilibrium
$\mathrm d\Pi/\mathrm dz = -\phi\Delta\rho\,g(z)$. One tube gives the whole curve. Reading
$\phi(z)$ non-destructively needs X-ray or NMR profiling; destructively, it needs careful
sectioning.

**Compressional rheology / filtration.** Apply a step pressure to a suspension above a
membrane and record the height against time. This is the most direct route because it
returns $K_p$ *and* the consolidation coefficient from one experiment, and the Buscall--White
and Landman--White framework for extracting both is standard in the literature.

**Confined oscillatory compression.** Compress a suspension in a rheometer between a plate
and a porous frit, oscillating at frequency $\omega$, and measure the normal force. Below the
drainage frequency the response is drained and returns $K_p$; above it the response is
undrained. The crossover frequency is $1/\tau_p$, which makes this a measurement of both.
The practical obstacle is that most rheometers do not have a normal-force resolution good
enough for a network modulus of order $10^2$ Pa.

\newpage

# 7. The shortcut: you do not need $\kappa$ and $K_p$ separately

Look again at (2.10). The relaxation time depends on $\kappa$ and $K_p$ only through their
product:

$$\tau_p = \frac{L^2}{D_p}, \qquad D_p = \frac{\kappa K_p}{\mu_f}.$$

$D_p$ is a single measurable quantity. Splitting it into a permeability and a modulus is a
modelling convenience, not a requirement. **If the goal is the relaxation time, measure
$D_p$ and stop.** Three routes, cheapest first.

**Dynamic light scattering at small $q$.** The collective diffusion coefficient is exactly
$D_p$, and DLS measures it directly from the decay rate $\Gamma = D_c q^2$ in the limit
$qa \ll 1$. The obstacle at $\phi = 0.5$ is multiple scattering: the sample is opaque. The
standard fixes are refractive-index matching (feasible for silica in a
water/glycerol/DMSO mix, though that changes $\mu_f$ and therefore $\tau_p$, so the
measurement must be corrected back), 3D cross-correlation DLS, or diffusing-wave
spectroscopy. **This is the measurement to ask a collaborator for**, and it is a
half-day on an instrument the lab may already have access to.

**A consolidation (filtration) experiment.** Load the suspension into a cell above a
membrane, apply a step pressure, record the piston height $h(t)$. The approach to
equilibrium is exponential with the time constant of (2.11), so a single fit returns $D_p$
for a known cell height. This measures the actual physical process rather than a proxy, and
it needs no optics. It is the most defensible number available.

**A centrifuge consolidation run.** The same experiment with centrifugal rather than
mechanical loading. Faster, and it reaches higher effective pressures, at the cost of a
non-uniform body force that has to be modelled.

All three measure the same $D_p$, and any two agreeing would be a genuine cross-check of the
kind this project currently lacks.

\newpage

# 8. The dimensionless groups, and what each is for

Four groups appear in this problem and they are routinely confused.

**The drainage number, $\mathrm{De} = \tau_p/\tau_{\rm imp}$**, with
$\tau_{\rm imp} = D_0/U_0$. Asks whether the pore pressure can relax during the event.
$\mathrm{De} \gg 1$ is undrained, so an imposed dilation is fully frustrated and generates
the largest particle pressure the network can carry. $\mathrm{De} \ll 1$ is drained, and
suction bleeds away as fast as dilation makes it. This is the group this document is
built on.

**The inertial--dilatant number,
$\mathrm{ID} = \rho\kappa U_0/(A\mu_f\Delta\phi D_0)$**, from Mobaseri *et al.*, Eq. 6.
Asks whether inertia or dilatancy wins the *force* balance. Note it contains no $K_p$: it
is a Reynolds--Darcy construct, and it compares $F_i \sim \rho U_0^2D_0^2$ with
$F_d \sim A\mu_fU_0D_0^3\Delta\phi/\kappa$. It answers a different question from De and the
two must not be conflated. The 2026-08-06 brief documents two problems with Eq. 6 as
printed, including that it does not evaluate to 1 at the regime boundary it is drawn at.

**The fluid-compressibility ratio, $\tau_\chi/\tau_{\rm imp} = \mu_f\chi_f D_0U_0/\kappa$.**
Asks only whether the pore fluid may be treated as incompressible in the Darcy--Reynolds
equation. Jérôme *et al.*'s own validity condition. 0.014 for cornstarch, 0.81 for the
lab's silica, which is marginal, and 17.6 with 0.1\% entrained gas, which is not.

**The Péclet number, $\mathrm{Pe} = 6\pi\mu_f(a/2)^3\dot\gamma/k_BT$.** Asks whether the
suspension is athermal. 274 for the silica, $7\times10^5$ for cornstarch.

A table of what each one decides:

| group | decides | silica | cornstarch |
|---|---|---:|---:|
| De | can pore pressure relax? | $1\times10^6$ | 7.2 |
| ID | inertia or dilatancy? | $6\times10^{-5}$ | $1\times10^{-2}$ |
| $\tau_\chi/\tau_{\rm imp}$ | is the pore fluid incompressible? | 0.81 | 0.014 |
| Pe | is it athermal? | 274 | $7\times10^5$ |

\newpage

# 9. Worked calculations

## 9.1 800 nm silica at $\phi = 0.50$

Inputs: $a = 800$ nm, $\phi = 0.50$, $\phi_c = 0.42$ (assumed, see §10),
$\mu_f = 10^{-3}$ Pa s, $\lambda_D = 10$ nm, $P_{\rm conf} = \sigma^\ast = 44.2$ Pa,
$D_0 = 2.5$ mm, $U_0 = 2.3$ m/s.

Permeability, (5.2):

$$\kappa = \frac{(0.50)^3(8\times10^{-7})^2}{180(0.50)^2} = 1.78\times10^{-15}\ \mathrm{m^2}$$

Network modulus, (6.4):

$$K_p = \frac{44.2\times8\times10^{-7}}{6\times10^{-8}} = 590\ \mathrm{Pa}$$

Storage, (2.9): $1/K_p = 1.70\times10^{-3}$ against $(1-\phi)\chi_f = 2.3\times10^{-10}$.
The skeleton carries all of it, to seven figures.

Diffusivity, (2.10):

$$D_p = \frac{1.78\times10^{-15}\times590}{10^{-3}} = 1.05\times10^{-9}\ \mathrm{m^2/s}$$

Relaxation, (2.12), at $L = 0.49D_0 = 1.23$ mm:

$$\tau_p = \frac{(1.23\times10^{-3})^2}{1.05\times10^{-9}} = 1.43\times10^{3}\ \mathrm{s}$$

against an impact time $\tau_{\rm imp} = D_0/U_0 = 1.09$ ms. So

$$\mathrm{De} = 1.3\times10^{6}.$$

Uncapped dilatant pressure, (4.3), with $\dot\gamma = U_0/D_0 = 920$ s$^{-1}$:

$$p_{\rm dil} = \frac{10^{-3}}{1.78\times10^{-15}}\times0.08\times(1.23\times10^{-3})^2\times920 = 6.2\times10^{7}\ \mathrm{Pa}$$

Ceiling, (4.5):

$$p_{\rm dil}^{\max} = \frac{590\times0.08}{0.50} = 94\ \mathrm{Pa}$$

The ceiling binds by six decades. Everything the material can do is bounded by 94 Pa, of
the same order as $\sigma^\ast = 44$ Pa. That coincidence is not independent, since
$P_{\rm conf}$ was set to $\sigma^\ast$; §9.4 removes the anchor.

## 9.2 Cornstarch at $\phi = 0.42$

Inputs: $a = 10$ µm, $\phi_c = 0.37$, $E_p = 10^8$ Pa, $\nu = 0.3$,
$P_{\rm conf} = 10^5$ Pa, $D_0 = 4$ mm, $U_0 = 5.1$ m/s.

$\kappa = 6.14\times10^{-13}$ m$^2$. Hertz, (6.5), gives $k_n = 83$ N/m at
$\delta/a = 2.3\times10^{-2}$, and the affine modulus gives $K_p = 1.10\times10^6$ Pa.
Then $D_p = 6.8\times10^{-4}$ m$^2$/s, and at $L = 1.96$ mm,

$$\tau_p = 5.7\ \mathrm{ms}, \qquad \mathrm{De} = 7.2.$$

Ceiling: $p_{\rm dil}^{\max} = 1.10\times10^6\times0.05/0.42 = 1.3\times10^5$ Pa, against an
uncapped Darcy value of $4.0\times10^5$ Pa. The ceiling binds, but only by a factor of 3,
and the resulting number sits in the $10^5$--$10^6$ Pa band Waitukaitis and Jaeger measured
for jammed cornstarch columns under impact.

**The contrast.** Silica's $D_p$ is $4.8\times10^5$ times smaller. That decomposes as
$\times257$ from permeability and $\times1854$ from network stiffness, but the second factor
is mostly an artefact of the two experiments operating at pressures 2000-fold apart. Put
both at $P_{\rm conf} = 10^5$ Pa and the moduli come out within 20\% of each other, and the
whole contrast collapses back onto the $a^2$ permeability factor. **The permeability
contrast is the material difference; the stiffness contrast says where each experiment
sits.**

## 9.3 The validation: Jérôme's glass beads

Inputs: $a = 170$ µm, $\phi = 0.58$, $\Delta\phi = 0.005$ (estimated),
$E = 70$ GPa, $\nu = 0.25$, $P_{\rm conf} = 10^4$ Pa, the low end of their measured suction.

$\kappa = 3.5\times10^{-11}$ m$^2$. Hertz gives $k_n = 4.5\times10^4$ N/m at
$\delta/a = 5.1\times10^{-5}$, so $K_p = 4.9\times10^7$ Pa, and

$$p_{\rm dil}^{\max} = \frac{4.9\times10^7\times0.005}{0.58} = 4.2\times10^{5}\ \mathrm{Pa}$$

against a measured $10^4$ to $10^5$ Pa. Within a decade, and above rather than below, which
is the expected direction for an affine upper bound. The dominant remaining uncertainty is
$\Delta\phi$, which is not reported in their paper: $\Delta\phi = 0.001$ would put the
prediction at $8\times10^4$ Pa, dead centre of the measurement.

Their drainage number is $\mathrm{De} = 0.003$. Glass beads drain freely, which is why their
suction is set by the Darcy build-up rate rather than by the ceiling.

## 9.4 The one prediction that is anchor-free

Divide (4.5) by the confining pressure. Using (6.4), $P$ cancels:

$$\boxed{\;\frac{p_{\rm dil}^{\max}}{P_{\rm conf}} = \frac{a}{6\lambda_D}\cdot\frac{\Delta\phi}{\phi}\;} \tag{9.1}$$

No $\sigma^\ast$, no $K_p$, no permeability. Only a particle size, a Debye length, and the
distance above the dilatancy onset. Below 1, frustrated dilation cannot raise the particle
pressure above what already drives it, so Reynolds dilatancy is not available as a
strengthening mechanism. Above 1, it is. Setting (9.1) to unity,

$$\phi^\ast = \frac{\phi_c}{1 - 6\lambda_D/a}. \tag{9.2}$$

| added salt | $\lambda_D$ | $\phi^\ast$ | ratio at $\phi=0.44$ | at $\phi=0.50$ |
|---|---:|---:|---:|---:|
| nominal 0 (residual) | 30 nm | 0.542 | 0.20 | 0.71 |
| 1 mM | 10 nm | 0.454 | 0.61 | 2.13 |
| 100 mM | 1 nm | 0.423 | 6.06 | 21.3 |

At the volume fraction the drops are actually made at, $\phi = 0.44$, the silica is on the
wrong side of the crossover at low salt. **Reynolds dilatancy cannot be the arrest
mechanism there unless the salt is high.** That is consistent with what the DNS found by an
entirely different route: at $\phi = 0.44$ the fitted Wyart--Cates closure cannot shear-jam
at any stress either, and a free-surface mechanism is required in any case because a local
constitutive law has no propagating-front solution.

$\phi^\ast$ moves from 0.54 to 0.42 across the salt range, which is a curve the existing
salt-screening series can be replotted against at no experimental cost, provided the $\phi$
axis is trusted. Worth flagging without leaning on it: at 1 mM, (9.2) gives
$\phi^\ast = 0.454$, and Shah, Arora and Driscoll report the first localised solid patches
at $\phi = 0.47$ on the same 830 nm particles. $\lambda_D = 10$ nm was a default rather
than a measurement on their samples, so that is a coincidence worth testing rather than
evidence.

\newpage

# 10. How each assumption fails

In rough order of how much damage it does.

**$\phi_c$ is unmeasured.** Every quantity containing $\Delta\phi$ scales with it, and
$\phi^\ast$ moves with it roughly one-for-one. Three anchors exist and they cannot all sit
on one $\phi$ axis: Shah's 0.47 for the first solid patches, the fitted $\phi_M = 0.5176$,
and the lab's own observation of fronts at 0.42. Until the $\phi$ axis is settled, every
absolute number here carries that spread. Nothing else on this list matters as much.

**$K_p$ is unmeasured for either system.** Bracketed here by $\lambda_D$ from 30 to 1 nm
for the silica and $E_p$ from $10^7$ to $10^9$ Pa for cornstarch. The undrained verdict
survives the whole bracket by at least two decades; the absolute ceilings do not survive it
at all.

**The affine modulus is an upper bound.** Real packings relax non-affinely, softening by
two to five times near isostaticity. Every ceiling quoted here is correspondingly high.

**The single-length scaling of (4.3).** Replacing the Poisson problem (4.2) by
$p \sim (\mu/\kappa)\Delta\phi L_d^2\dot\gamma$ is a scaling, not a solution. Mobaseri
*et al.*'s own $F_d$ carries the same substitution and it is where their factor-of-150
discrepancy most likely lives. If a factor of 3 ever matters, solve (4.2) on the actual
contact geometry.

**Carman--Kozeny at low salt.** §5.2: the double layer fills the pore throat at the
salt-free end, so the effective permeability is below the geometric one. Direction known,
magnitude not.

**Steady-state constitutive assumptions.** (4.1) presumes the dilatancy law responds
instantaneously. Under a 1 ms impact with $\mathrm{De} \sim 10^6$ that is not obviously
true, and a rate-dependent $\phi_{\rm eq}$ would change the transient.

**Monodispersity.** Assumed throughout, and true for this silica, which is why the analysis
is cleaner for it than it would be for cornstarch. Cornstarch granules are polydisperse and
faceted, so its $\kappa$ carries an unquantified shape factor.

\newpage

# 11. What to do about it, cheapest first

**A conductivity reading on the salt-free samples.** Minutes. Removes a factor spanning a
decade in $\lambda_D$, which propagates directly into $K_p$, $\phi^\ast$ and the
electroviscous correction. There is no reason not to have this number.

**$\phi_c$ and the $\phi$ axis, by DLS and $\zeta$.** Half a day, already the top item on the
open-problems slate, and already gating three other things. It is the same measurement in
all four cases.

**$D_p$ by DLS at small $q$, or by a consolidation cell.** §7. Returns the pore relaxation
time directly, without needing $\kappa$ and $K_p$ separately. If only one new measurement
is possible, this is the one, because it turns the central quantity of this document from an
estimate into data.

**Replot the salt series against $\phi^\ast(\lambda_D)$.** No new experiments. Tests (9.2)
against data the lab already owns.

**A settling column.** A week of passive observation gives $\kappa$ from (5.3) and, with
$D_p$, closes the system: $K_p = \mu_f D_p/\kappa$ with no micromechanical model at all.
That combination would make every number in §9 a measurement rather than an estimate.

\newpage

# 12. Symbols

| symbol | meaning | units |
|---|---|---|
| $a$ | particle **diameter** | m |
| $\phi$, $\phi_c$, $\Delta\phi$ | solid fraction, dilatancy onset, difference | -- |
| $\kappa$ | permeability, $(1-\phi)^3a^2/180\phi^2$ | m$^2$ |
| $\mu_f$, $\chi_f$ | pore-fluid viscosity, compressibility | Pa s, Pa$^{-1}$ |
| $K_p$ | drained bulk modulus of the particle network | Pa |
| $S$ | storage coefficient, $1/K_p + (1-\phi)\chi_f$ | Pa$^{-1}$ |
| $D_p$ | pore-pressure diffusivity $=\kappa K_p/\mu_f$; also the collective diffusivity | m$^2$/s |
| $L$, $L_d$ | drainage length, dilation length | m |
| $\tau_p$, $\tau_1$ | relaxation scaling $L^2/D_p$; exact slowest mode $0.405\,\tau_p$ | s |
| $\tau_{\rm imp}$ | impact time $D_0/U_0$ | s |
| $\mathrm{De}$ | drainage number $\tau_p/\tau_{\rm imp}$ | -- |
| $p$, $\sigma'$ | pore pressure, effective stress | Pa |
| $\psi$, $\alpha$ | dilatancy angle, $\tan\psi = \alpha\Delta\phi$ | -- |
| $\lambda_D$ | Debye screening length | m |
| $r_{\rm throat}$ | pore-throat radius, $0.0775\,a$ | m |
| $k_n$, $Z$ | contact normal stiffness, coordination number | N/m, -- |
| $\sigma^\ast$ | Wyart--Cates onset stress, fitted 44.2 Pa | Pa |

## Formula card

$$\kappa = \frac{(1-\phi)^3a^2}{180\phi^2}
\qquad
S = \frac{1}{K_p} + (1-\phi)\chi_f
\qquad
D_p = \frac{\kappa}{\mu_f S} \simeq \frac{\kappa K_p}{\mu_f}$$

$$\tau_p = \frac{\mu_f L^2}{\kappa K_p}
\qquad
\tau_1 = \frac{4}{\pi^2}\tau_p
\qquad
\mathrm{De} = \frac{\tau_p U_0}{D_0}$$

$$p_{\rm dil} = \frac{\mu_f}{\kappa}\Delta\phi L_d^2\dot\gamma = K_p\Delta\phi\,\dot\gamma\tau_p
\qquad
p_{\rm dil}^{\max} = \frac{K_p\Delta\phi}{\phi}$$

$$K_p^{\rm EDL} = \frac{Pa}{6\lambda_D}
\qquad
\frac{p_{\rm dil}^{\max}}{P} = \frac{a}{6\lambda_D}\frac{\Delta\phi}{\phi}
\qquad
\phi^\ast = \frac{\phi_c}{1-6\lambda_D/a}$$

# 13. Sources

Biot, *J. Appl. Phys.* **12**, 155 (1941).
Terzaghi, *Erdbaumechanik auf bodenphysikalischer Grundlage* (1925).
Carman, *Trans. Inst. Chem. Eng.* **15**, 150 (1937).
Carnahan & Starling, *J. Chem. Phys.* **51**, 635 (1969).
Reynolds, *Phil. Mag.* **20**, 469 (1885).
Pailha & Pouliquen, *J. Fluid Mech.* **633**, 115 (2009).
Jérôme, Vandenberghe & Forterre, *Phys. Rev. Lett.* **117**, 098003 (2016).
Athani, Metzger, Forterre & Mari, *Phys. Rev. Fluids* **10**, 043301 (2025).
Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng, *Phys. Rev. Lett.* **136**, 148201 (2026).
Makse, Gland, Johnson & Schwartz, *Phys. Rev. Lett.* **83**, 5070 (1999).
Agnolin & Roux, *Phys. Rev. E* **76**, 061302 (2007).
Buscall & White, *J. Chem. Soc. Faraday Trans.* **83**, 873 (1987).
Waitukaitis & Jaeger, *Nature* **487**, 205 (2012).
Shah, Arora & Driscoll, *Commun. Phys.* **5**, 222 (2022).
Ishibashi *et al.*, *Colloids Surf. A* **732**, 139148 (2026).
Wyart & Cates, *Phys. Rev. Lett.* **112**, 098302 (2014).

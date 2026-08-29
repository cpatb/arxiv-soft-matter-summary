# Solver landscape: axisymmetric VOF + Wyart–Cates shear-thickening suspension solver

Compiled 2026-08-09. Sources cited inline. Where a claim could not be verified from a
primary source it is marked **UNVERIFIED**. Statements are labelled [source states] vs
[inference] where the distinction matters.

Status: COMPLETE. Sections: 1 competitive landscape, 2 commercial CFD, 3 industrial
applications, 4 funding, then an **EXECUTIVE VERDICT at the end of the file** — read that
first if short on time.

Three things in here contradict the assumptions the task was written with; they are in
§1.9, §2.3 and §3 (rank 3). Read those before drafting any positioning claim.

---

## TASK 1 — COMPETITIVE LANDSCAPE

### 1.1 Baumgarten & Kamrin, PNAS 116, 20828 (2019); J. Rheol. 64, 367 (2020)

**What it is.** A 3-D two-phase (mixture-theory) continuum constitutive model for dense
fine-particle suspensions, implemented in an adaptation of the **material point method
(MPM)**. Validated across many geometries.
DOI 10.1073/pnas.1908065116; open full text at
https://pmc.ncbi.nlm.nih.gov/articles/PMC6800318/

**Critical correction to the assumed positioning.** Their microstructural variable *f*
("the fraction of frictional contacts within the granular skeleton") is **not algebraic** —
it carries its own evolution equation of the form

    f_dot / (K0 * gammadot) = H(f_m - f) - S f

with hardening term H and softening term S, and *f* modulates both the jamming fraction
phi_m = phi_hat_m(f) and the dilation parameter a = a_hat(f).
[source states — quoted from the PMC full text]

**Therefore: "a frictional-contact-fraction state variable with memory in a continuum
solver" is NOT novel by itself.** Baumgarten & Kamrin got there in 2019. What differs:

- Their evolution is **strain-clock** driven (the rate term is K0*gammadot), so it has a
  characteristic *strain*, not a characteristic *time*. The solver's
  `-(gammadot/gamma0 + 1/t_r)` form carries **both** a strain clock and an absolute-time
  clock t_r, so it relaxes even at gammadot -> 0. That is the physically load-bearing
  difference for un-jamming after impact, when the flow has stopped. [inference, but a
  direct structural comparison of the two equations]
- Their constants are fitted per-material; the solver's gamma0 and t_r are taken from
  independent experiments (Richards PRL 123 038004; Han PRF 3 073301). [inference from
  the task brief — I did not re-derive their fitting procedure]
- **Surface tension is absent.** It is not mentioned anywhere in the PNAS paper.
  Free surfaces are handled only as zero-traction MPM boundaries: "The upper surface of
  the simulated mixture is modeled as 'free' (zero traction)." That is a *stress-free
  boundary*, not a capillary interface — no Laplace pressure, no Weber number, no
  breakup/splash physics. [source states]
- Stated open items include fabric-tensor contributions for "rheologically chaotic flows",
  material-property dependence of tau*, eta_B and K0, and numerical pathologies
  ("kinematic locking", "ringing instability"). [source states]

**Verdict on this competitor:** the closest prior art on the *rheology* side, and it must
be cited and distinguished honestly. The differentiator is not "advected relaxing f" — it
is (a) absolute-time relaxation in addition to strain relaxation, and (b) a genuine
sharp-interface capillary free surface.

### 1.2 Mobaseri, Kumar & Cheng, PNAS 122, e2500163122 (2025)

arXiv:2503.23555. Directly adjacent: **drop impact, Basilisk, VOF, non-Newtonian.**

- Method [source states]: "We perform direct numerical simulations of drop impact based on
  the volume-of-liquid solver Basilisk, which solves the time-dependent mass and momentum
  conservation equations with second-order accuracy in both space and time."
- Constitutive model [source states]: "we employ the Carreau model, a constitutive model
  for generalized Newtonian liquids." No memory, no thixotropy, no viscoelasticity, no
  jamming variable. Viscosity ratio beta = mu_inf/mu_0 spanned 0.01–100.
- **The explicit gap, in their own words** [source states]: the liquid-to-solid transition
  via discontinuous shear thickening and shear jamming "is not captured by the current
  study, which focuses on continuous shear thickening at relatively low cornstarch mass
  fractions." Experiments were held at 10–40 wt% specifically to stay continuous.
- Also excluded by their own statement: viscoelasticity, viscoplastic/yield-stress drops
  ("beyond the scope of this study").

**Verdict:** same platform (Basilisk), same problem (drop impact), same material
(cornstarch), and the authors state in print that the shear-jamming transition is the
thing they do not model. This is the single strongest citation for the solver's premise.

### 1.3 Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng, arXiv:2601.12642 (Jan 2026)

"Inertia-Dilatancy Interplay Governs Shear-Thickening Drop Impact." Note this is
**experiment + a reduced analytical model**, not a solver: high-speed photography plus
direct force measurement, unified by "the classic drop-impact theory with the
Reynolds-Darcy mechanism for dilatancy." [source states, abstract]

Its key result is directly load-bearing for the solver's premise: they report a regime in
which high-concentration mixtures are **liquid-like at impact onset (high shear rate) and
only transition to solid-like at later times as shear rate falls**. A time-delayed,
memory-carrying jamming variable is exactly the object needed to reproduce that; an
instantaneous eta(gammadot) closure cannot produce it by construction. [inference — the
paper does not claim a solver is needed, but the regime is non-instantaneous]


### 1.4 Cui, Wang, Ye et al., Rheologica Acta 62, 225–238 (2023)

"Numerical simulation of shear jamming in a shear thickening fluid under impact",
DOI 10.1007/s00397-023-01391-2.

Bibliographic details and method summary confirmed from the publisher landing page and
NASA ADS (2023AcRhe..62..225C). Their own abstract-level description: "The STF was defined
with strain rate-dependent viscosity and compressibility, with the air-STF interface
modeled using the volume of fluid method and impact penetration reproduced by a dynamic
mesh method." [source states]

The limitation quoted in the brief — that "the phase transition between liquid and jamming
state is not modelled and the relaxation process of the shear jamming after impact cannot
be simulated" — is **taken as previously verified by the requester from the full text.**
I could **not independently re-verify it**: the Springer full text is paywalled and
redirects to an IdP login. Treat the quote as verified-by-user, not verified-by-me.
Before it goes into a pitch or a paper, re-check the exact wording and page number against
the PDF. **UNVERIFIED (by me).**

Note the structural point that survives regardless: their viscosity is an algebraic
function of strain rate. A strain-rate-algebraic closure cannot represent hysteresis or
post-impact relaxation, because the state is a pure function of the instantaneous
kinematics. [inference, but a mathematical one]

### 1.5 Angerman, Sandnes, Seto & Ellero — SPH with a Wyart–Cates scalar

Two papers: Phys. Fluids 36, 033103 (2024), arXiv:2311.13759 ("Microstructural SPH model
and simulations of DST fluids"); and J. Fluid Mech. 1016, A33 (2025),
doi:10.1017/jfm.2025.10178, arXiv:2408.08285.

**This is the closest competitor on the scalar-field formulation, and it is close.**
Their JFM model [source states, Eq. 6]:

    Df(x,t)/Dt = k_f * gammadot * ( f_hat(x,t) - f(x,t) ) + alpha * grad^2 f(x,t)

with f "the fraction of particles in frictional contact at a given position in space",
f_hat the Wyart–Cates steady-state, and D/Dt the material derivative — i.e. the field is
**advected**, and relaxes toward the WC steady state. So "advected, relaxing WC contact
fraction in a continuum solver" is **occupied ground**, in SPH, since 2024.

What they do **not** have [source states]:
- Geometry is **planar simple shear** between walls: "We present numerical simulations of a
  scalar model for DST fluids in a planar simple shear using the Smoothed Particle
  Hydrodynamics (SPH) approach." Box of gap h = 0.01 m, periodic in flow direction.
  **No free surface, no interface, no surface tension anywhere.**
- "The primary limitation of our implementation is the lack of vorticity direction."
- 2-D planar only; no tensorial microstructure; vorticity-direction instabilities left to
  future work.

Structural difference that remains: their relaxation rate is **k_f * gammadot only** — a
pure strain clock. At gammadot -> 0 the microstructure freezes. The solver's extra
`1/t_r` term is what lets a jammed region un-jam after the flow has stopped, which is
precisely the post-impact relaxation problem. They instead add a **spatial diffusion**
term alpha*grad^2 f (non-locality) that the solver does not have. [inference from
side-by-side comparison of the two equations]

### 1.6 Athani, Metzger, Forterre & Mari, PRF 10, 043301 (2025)

arXiv:2408.15130. Particle-based **DEM under pressure-imposed conditions**, plus a
continuum model built on a Reynolds-like dilatancy law coupled to Wyart–Cates.
No free surface, no air–liquid interface, no surface tension anywhere in the paper.
[source states / absence confirmed]

Their closing claim is the clearest statement in the literature that the target ground is
open [source states, abstract]: "Beyond globally jammed suspensions, our results enable the
modeling of inhomogeneous flows where shear jamming is local, e.g. under impact, which
eludes usual shear thickening rheological laws."

They **do not simulate an impact.** "our results enable the modeling" is a stated
capability for future work, not a result in the paper. [source states + inference]

Their physical mechanism — Reynolds dilation coupled to Darcy backflow setting the finite
stress — is **not** in the solver as described. That is a real competitive exposure: see
§1.9.

### 1.7 Field-level review: Ayar, Siddani, Srivastava & Singh, arXiv:2603.12410 (June 2026)

"Recent Modeling Advances in Dense Suspension Rheology" — a current review. Useful because
it states the gap independently [source states]:

- "Many state-of-the-art models assume idealized simple shear conditions without explicitly
  resolving wall slip, localization, free surfaces, or complex boundary conditions."
- "In the Wyart–Cates framework, the fraction of frictional contacts f(sigma) is a
  stress-activated order parameter." — the review treats WC only in its **static** form;
  it does not survey advected/relaxing implementations.
- The review's method table covers Stokesian Dynamics, ASD, FSD, DEM+hydrodynamics, Force
  Coupling, Fictitious Domain, Immersed Boundary, two-fluid continuum, and the Suspension
  Balance Model. **None of these is a free-surface method.**

Caveat: the "no method integrates free-surface tracking with frictional contact mechanics"
formulation is my summariser's synthesis of the review's gap statements, not a verbatim
sentence. Quote only the two verbatim lines above. **Partially UNVERIFIED.**

### 1.8 Also checked, not competitors

- Shi, Hu & Zhao, arXiv:2602.05669 (Feb 2026), "Boundary compliance selects heterogeneous
  dynamics in shear-thickening suspensions" — **experimental** (cornstarch under a silicone
  oil layer, rheometry + visualisation). Relevant as validation data for a two-fluid
  interface problem; not a solver. [source states, abstract]
- Slutzky, Pelosse, van der Naald & Jaeger, arXiv:2511.17742, Phys. Rev. E (2026),
  "Fracture and failure of shear-jammed dense suspensions under impact" — **experimental**.
  Directly relevant: they vary **surface tension** and find "this likelihood [of radial
  fracturing] ... increases when the solvent's viscosity or surface tension are reduced."
  [source states] This is the strongest published evidence that capillarity is not a
  decoration in impacted dense suspensions. It is also a warning: at high stress the
  jammed state **fractures**, which a viscous continuum closure cannot represent.
- Angerman, Kumar, Sandnes, Seto & Ellero, Phys. Fluids 36, 033103 (2024), abstract
  [source states]: "Despite the recent interest in the discontinuous shear-thickening (DST)
  behaviour, few computational works tackle the rich hydrodynamics of these fluids."

### 1.9 VERDICT — is the ground unoccupied?

**Partly. The ground is more crowded than the framing implies, and the honest claim is
narrower than "first solver for the jamming transition with memory."**

Occupied already:
1. A **frictional-contact-fraction state variable with its own evolution equation** in a
   continuum solver — Baumgarten & Kamrin 2019/2020 (MPM), and Angerman et al. 2024/2025
   (SPH), the latter explicitly advected and explicitly Wyart–Cates.
2. **Wyart–Cates + Reynolds dilatancy as a continuum closure aimed at impact** — Athani
   et al. 2025 state exactly that ambition.
3. **VOF + drop impact + cornstarch on Basilisk** — Mobaseri, Kumar & Cheng 2025.

Genuinely unoccupied, as far as this search reaches:
- **The intersection.** No published work combines (a) a Wyart–Cates contact-fraction
  field with its own dynamics, (b) a sharp two-phase interface with surface tension, and
  (c) impact/free-surface flow. Each competitor holds two of the three at most:
  Baumgarten–Kamrin has (a) + free surfaces but no capillarity; Angerman has (a) with no
  interface at all; Mobaseri has (b)+(c) with no (a) and says so in print.
- **Absolute-time relaxation.** Both Baumgarten–Kamrin and Angerman relax f on a *strain*
  clock (rate proportional to gammadot). Neither can un-jam a region that has stopped
  moving. The `1/t_r` term is the distinguishing piece, and its experimental provenance
  (Richards PRL 123 038004; Han PRF 3 073301) is the defensible part of the claim.
- **Un-jamming / hysteresis as a *simulated* outcome.** Every source above either
  models the onset only, or states the reverse transition is out of scope (Cui;
  Mobaseri).

Claims to **avoid** in a pitch, because a referee will hold them against you:
- "First to treat the frictional contact fraction as a dynamical variable" — false,
  Baumgarten & Kamrin 2019.
- "First advected Wyart–Cates field" — false, Angerman et al. 2024.
- "No one has modelled shear jamming under impact" — Cui et al. 2023 attempted it, and
  Athani et al. 2025 built the continuum closure for it.

The safe, defensible claim: **the first sharp-interface, surface-tension-resolved
free-surface solver in which the Wyart–Cates contact fraction is advected and relaxes on
an experimentally measured absolute timescale, so that jamming and un-jamming are both
outcomes rather than inputs.** Every clause in that sentence is doing work.

### 1.10 Competitive exposure to fix or acknowledge

1. **No Reynolds dilatancy / Darcy backflow.** Athani et al. (PRF 2025) show the stress in
   a locally shear-jammed region is set by dilation + pore-pressure backflow, not by the
   steady-state flow rule. A single-phase eta(phi, f) closure with fixed phi cannot
   produce that. Expect this as referee objection #1. [inference from their result]
2. **No fracture.** Slutzky et al. (PRE 2026) show the jammed state fails by cracking above
   a stress threshold. A viscous closure will keep thickening instead.
3. **Non-locality.** Angerman et al. found the local model gives "stress-splitting
   instability" without a non-local term; they add alpha*grad^2 f. A purely local f may be
   ill-posed in the DST branch. Worth checking whether the solver's S-curve region is
   regularised by anything (mesh, surface tension, the t_r term) or is silently
   mesh-dependent. [inference — this is a genuine numerical risk, not a criticism of a
   result I have seen]
4. **Axisymmetry.** Jamming rays, radial cracks, and banding are non-axisymmetric.
   Axisymmetry is a real ceiling on the phenomena reachable.


---

## TASK 2 â€” WHAT COMMERCIAL / MAINSTREAM CFD CAN ACTUALLY DO

**Headline: none of the mainstream packages ships a jamming model. All ship
generalized-Newtonian closures; at least one ships viscoelastic constitutive equations; and
exactly one ships a time-evolving microstructure scalar â€” which is a thixotropy (breakdown)
model, not a jamming model.**

Verification caveat up front: ANSYS and Siemens keep current documentation behind a
customer login. Where I could only reach third-party or legacy mirrors I say so. Do not put
an unverified capability claim about a competitor's product into a pitch.

### 2.1 ANSYS Fluent

Native non-Newtonian viscosity models (material panel): **power law, Carreau /
Carreau-Yasuda, Cross, Casson, Herschel-Bulkley, Bingham**. All are algebraic functions of
the local strain-rate magnitude â€” no history, no state variable.
Sources: legacy Fluent User's Guide section "Viscosity for Non-Newtonian Fluids"
(mirror: https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node297.htm â€” I
could not fetch this over TLS, only see it summarised in search results) and the
third-party summary https://cfdland.com/non-newtonian-fluids-in-ansys-fluent/ .
**The current official Fluent Theory Guide is login-gated and I could not fetch it â€” treat
the exact model list as approximately, not exactly, verified. UNVERIFIED against current
docs.**

Viscoelastic models (PTT, Giesekus, FENE): the third-party source states plainly, "to
simulate these advanced models, you would need to write UDFs (User Defined Functions) ...
they are not directly available in the standard model library." I could **not** confirm
from Ansys's own release notes whether recent releases added native viscoelasticity to
Fluent; searches of Ansys release-highlight pages returned nothing on this.
**UNVERIFIED â€” check before asserting.** (ANSYS **Polyflow**, the separate
polymer-processing product, has long carried differential viscoelastic models; it is a
distinct product from Fluent and is not a free-surface impact solver.)

**Verdict:** anything beyond generalized Newtonian is a UDF you write yourself. A
Wyart-Cates field with relaxation is implementable as a user-defined scalar plus a UDF
viscosity â€” no published example found. No jamming model exists in the box.

### 2.2 Siemens Simcenter STAR-CCM+

Generalized-Newtonian models are standard and Siemens markets viscoelastic material
modelling. I could **not** obtain the exact constitutive-model list from Siemens'
documentation â€” the docs portal is gated, and the release-note blogs reachable publicly did
not enumerate Giesekus / PTT / FENE-P. **UNVERIFIED.** No evidence of any jamming or
contact-fraction model; that is absence of evidence, and should be stated as such.

### 2.3 COMSOL Multiphysics â€” the closest commercial analogue, take it seriously

The **Polymer Flow Module** product page lists, verbatim:

- Inelastic: "Power law, Carreau, Carreau-Yasuda, Cross, Cross-Williamson, Sisko, Ellis,
  Bingham-Papanastasiou (Viscoplastic), Casson-Papanastasiou (Viscoplastic),
  Herschel-Bulkley-Papanastasiou, Robertson-Stiff-Papanastasiou,
  DeKee-Turcotte-Papanastasiou, **Houska thixotropy (Thixotropic)**"
- Viscoelastic: "Oldroyd-B, Giesekus, FENE-P, FENE-CR, Linear Phan-Thien-Tanner (LPTT),
  Exponential Phan-Thien-Tanner (EPTT), Rolie-Poly"
- Interface methods: "Level Set", "Phase Field", "Moving Mesh"

Sources: https://www.comsol.com/polymer-flow-module and
https://doc.comsol.com/6.3/doc/com.comsol.help.polymer/polymer_ug_fluidflow.06.52.html

**The Houska thixotropy model is the one commercial model with memory.** Per the COMSOL 6.3
docs it carries a structure field xi obeying its own differential equation with rebuild and
breakdown coefficients, and the consistency and yield stress depend linearly on xi.
[source states]

This cuts both ways:
- *Against you:* an industrial CFD lead can say "COMSOL already couples a time-evolving
  structure parameter to a level-set free surface." Structurally it is the same class of
  model, and that objection is fair.
- *For you:* Houska is a **thixotropy** model â€” structure breaks down under shear and
  rebuilds at rest, so viscosity **falls** with stress. Wyart-Cates has the opposite sign:
  frictional contacts **build** under stress, and viscosity diverges at a stress-dependent
  jamming point phiJ(f). Houska has no phiJ, no divergence, and no volume fraction. It
  cannot represent shear jamming. [inference from the two model forms â€” a model-structure
  argument, not a quotation]

COMSOL has no suspension/jamming model: the Polymer Flow Module page names "colloidal
suspensions" only as an application example, with no particle-migration or jamming physics.

**Consulting-house evidence that industry wants this:** Veryst Engineering publishes a
COMSOL case study, "Shear Jamming in Dense Suspensions", using "a phase transport and
mixture model", distinguishing "lubricated" from "frictional" particles and tracking "the
fraction of frictional contacts". Their named applications: "the blockage of inkjet printer
nozzles in additive manufacturing processes, the clogging in medical injection devices by
protein aggregates, and the occlusion of pores by organic contaminants in water filtration
membranes." [source states â€”
https://www.veryst.com/case-studies/shear-jamming-dense-suspensions]

Read that carefully. It is simultaneously (a) proof of commercial demand, from a firm that
sells engineering simulation to industry, and (b) proof that a Wyart-Cates-flavoured model
has already been hand-built inside commercial CFD by a consultancy. It is a custom build,
not a product feature, and there is no indication of a free surface in it â€” but do **not**
claim "no one has done this in commercial CFD."

### 2.4 OpenFOAM

Core `viscosityModel` classes: **Newtonian, powerLaw, BirdCarreau, CrossPowerLaw, Casson,
HerschelBulkley, strainRateFunction**. All algebraic in strain rate; `strainRateFunction`
accepts an arbitrary eta(gammadot) but still carries no history.
Sources: OpenFOAM v13 User Guide section 8.3 "Transport/rheology models",
https://doc.cfd.direct/openfoam/user-guide-v13/transport-rheology ; API guide
https://www.openfoam.com/documentation/guides/latest/api/classFoam_1_1viscosityModel.html

**rheoTool** (Pimenta & Alves, https://github.com/fppimenta/rheoTool) is the serious
open-source rheology extension: GNF plus differential viscoelastic models (Oldroyd-B,
Giesekus, PTT variants, FENE-P/FENE-CR), Papanastasiou regularisation for yield-stress
models, and a two-phase solver `rheoInterFoam` in which "any GNF or viscoelastic model can
be used for each phase", described as under development. Thixotropic models exist in this
broader ecosystem (the HiGFlow package is described as carrying "viscoelastic models that
incorporate non-constant viscosity and time-dependent phenomena (such as thixotropy)").
**No jamming / contact-fraction model in any of them.**

### 2.5 Particle-based commercial routes

CFD-DEM couplings (Rocky DEM + Fluent, Altair EDEM + Fluent, CPFD Barracuda) can in
principle resolve frictional contacts directly, but at dense-suspension volume fractions
with lubrication films the cost is prohibitive, and they are not used for free-surface
impact of cornstarch at engineering scale. **I did not verify any vendor claim here â€” this
is inference from the cost discussion in the Ayar et al. review (arXiv:2603.12410) and must
be checked before being asserted. UNVERIFIED.**

### 2.6 Commercial verdict, one sentence

Every mainstream package gives you eta(gammadot); COMSOL alone gives you a structure
variable with its own ODE (Houska thixotropy â€” wrong sign for jamming, no phiJ, no volume
fraction); nobody ships a jamming transition; and the one commercial-adjacent Wyart-Cates
implementation found (Veryst, in COMSOL) is a bespoke consulting build with no free surface.


---

## TASK 3 â€” INDUSTRIAL APPLICATIONS

Scored on (fit to what the solver already does) x (size of the pain point). "Fit" is
judged against the actual capability: **axisymmetric**, two-phase VOF with **surface
tension**, dense STF with a jamming/un-jamming memory variable, Herschel-Bulkley and
Oldroyd-B/FENE-P branches, textured substrate. Axisymmetry and single-phase-suspension
(no separate particle phase) are the two limits that decide most of the ranking.

A note on market figures: the only public numbers for STF-specific sectors come from
commercial market-research publishers (GrowthMarketReports, DataIntelo, GMInsights,
FutureMarketInsights). **These are not peer-reviewed and their methodologies are not
published. Cite them as "vendor estimates" or not at all.** Where I have no defensible
number I say so rather than inventing one.

### Rank 1 â€” Direct ink writing / dense-paste extrusion additive manufacturing

**Fit: very high. Pain: high and growing.**

- *Process:* extrusion of highly loaded ceramic, metal, battery, bioink or energetic pastes
  through a fine nozzle, deposited as a free filament onto a substrate or a previous bead.
- *What goes wrong:* nozzle clogging, filament breakup, die swell, filament-shape error,
  phase separation (liquid extruded preferentially, leaving the solids behind).
  A 2025-26 ceramics DIW review states that "agglomerates that act as stress concentrators
  induce nozzle clogging, disrupt filament continuity, and generate pores or microcracks
  during sintering", and that "High particle loading and insufficient particle dispersion
  or stability can be responsible for shear-thickening effects at high shear, which may
  disrupt ink flow during printing." [source states â€” via search summary of
  https://www.mdpi.com/2571-6131/9/2/16 and
  https://www.sciencedirect.com/science/article/pii/S2666539525000756 ; I did not read
  these in full, so treat as indicative. PARTIALLY UNVERIFIED]
- *Authoritative statement of the gap:* Sauret, Ray & Compton, **Annu. Rev. Fluid Mech. 58,
  413-442 (2026)**, "Fluid Mechanics Challenges in Direct-Ink-Writing Additive
  Manufacturing" (abstract verbatim, via NSF PAR record 10680843): "The rheology of fluids
  and suspensions governs the flow behavior through narrow nozzles, posing questions about
  extrudability, confined flow dynamics, and clogging mechanisms. Downstream, the formation
  and deposition of extruded filaments involve extensional flows and potential
  instabilities, while postdeposition dynamics introduces complexities related to yield
  stress and structural stability. ... developing fundamental fluid mechanics frameworks is
  essential to replace trial-and-error approaches with predictive design methodologies".
  Full text is embargoed to Jan 2027; only the abstract was reachable.
- *Why the fit is good:* a nozzle is natively axisymmetric; the extrudate has a genuine
  capillary free surface; the failure mode is literally "the material jams, then must
  un-jam"; the Herschel-Bulkley branch already covers the post-deposition
  shape-retention requirement; the textured-substrate capability maps onto printing onto a
  previously deposited bead.
- **Physics the solver would need to add, in priority order:**
  1. **Wall slip.** Paste extrusion is slip-dominated at the die wall; a no-slip
     axisymmetric pipe flow of a jamming suspension will over-predict pressure drop badly.
     This is the single highest-value addition for this market. [inference, standard paste
     rheology]
  2. **Shear-induced particle migration** (a transported phi field / suspension balance).
     Phase separation and the liquid-first extrusion failure mode are migration, and the
     current model holds phi fixed.
  3. **Yield stress AND jamming simultaneously**, not as alternative branches.
  4. 3-D, eventually, for cross-sectional filament shape and bead-to-bead coalescence;
     axisymmetric gets you the nozzle and the die swell but not the deposited road.
  5. Non-isothermal / solvent loss only for specific inks.

### Rank 2 â€” High-solids coating and atomization (battery electrode slurries first)

**Fit: high on the free-surface physics, medium on geometry. Pain: very large.**

- *Processes:* slot-die and roll coating of Li-ion electrode slurries; spray coating and
  atomization of high-solids paints and functional coatings; agricultural spray.
- *What goes wrong:* shear thickening at coating/nozzle shear rates destabilises the
  coating bead. From the coating literature: "shear thickening has the opposite effect of
  shear thinning in promoting air entrainment", and "shear thickening particle dispersions
  [cause] coating failure through ribbing instability, where critical web velocity for
  defect onset decreases with increasing particle concentration and fluid viscosity."
  [source states â€” via search summary of the slot-die/battery-slurry coating literature,
  incl. Reynolds et al., Energy Technology 10, 2200545 (2022),
  https://onlinelibrary.wiley.com/doi/full/10.1002/ente.202200545 . I did not read these
  in full. PARTIALLY UNVERIFIED â€” verify before quoting.]
- *Shear-rate context (agricultural spray, and a good general number):* drift-control
  adjuvants see "continuous shear forces of between 50 s-1 during recirculation in the
  spray tank, reaching up to 10^5 s-1 at the nozzle outlet" â€” Hewitt et al., Pest Manag.
  Sci. (2024), https://scijournals.onlinelibrary.wiley.com/doi/10.1002/ps.8255 .
  10^5 s-1 is far above the DST onset for most dense suspensions, so the thickened branch
  is the operating branch, not an edge case. [source states the shear rates; the inference
  about which branch is mine]
- *Why the fit is only partly good:* the two headline defects â€” **ribbing** in slot-die
  coating and **ligament/satellite formation** in atomization â€” are intrinsically
  three-dimensional and non-axisymmetric. An axisymmetric solver can do a single spray
  orifice, a jet, and drop-scale breakup, but cannot do ribbing at all.
- **Physics the solver would need:**
  1. **3-D.** Non-negotiable for ribbing and for realistic atomization. This is the
     largest single lift in this whole document.
  2. Moving contact line / dynamic wetting for coating-bead problems.
  3. Air-phase inertia and density ratio at spray velocities (the VOF is already
     two-phase, so this is parameter range rather than new physics).
  4. Solvent evaporation / drying for the coating-to-film step; non-isothermal.
  5. Particle migration (the coated film's through-thickness composition is a migration
     outcome and is what actually sets electrode performance).

### Rank 3 â€” STF impact protection: body armour, sports padding, industrial PPE

**Fit: high on physics, partial on geometry and on whether a free surface is even present.
Pain: real, and the current state of engineering practice is weak.**

- *What goes wrong today:* design is empirical. The current engineering-simulation approach
  is to replace the fluid entirely: "Engineering-level numerical modeling of shear-thickening
  fluid layers under low-energy impact using **solid-equivalent surrogate models**",
  Mech. Adv. Mater. Struct. (2026), doi:10.1080/15376494.2026.2646291 â€” the title alone
  makes the point. **The article is paywalled (403) and I could not read it; I am relying
  on the title and the search index. UNVERIFIED.**
- *Scientific pull is strong and current:* Slutzky, Pelosse, van der Naald & Jaeger,
  Phys. Rev. E (2026) / arXiv:2511.17742 impact dense suspensions with a rod, vary impact
  speed over orders of magnitude, and find failure by fracture whose likelihood "increases
  when the solvent's viscosity or **surface tension** are reduced." That is a
  free-surface impact experiment with capillarity as a control parameter â€” an almost exact
  validation target for this solver.
- *Market:* vendor estimates put "shear-thickening fluid body armor" at roughly
  USD 1.1-1.2 billion (2024) with ~8.7% CAGR, and the narrower "liquid armor materials"
  segment at USD 82.8 million (2025) growing to USD 222.7 million (2035).
  **These two numbers are mutually inconsistent by an order of magnitude, which tells you
  what these reports are worth. Do not put either in a pitch without caveating it.**
  Sources: https://growthmarketreports.com/report/shear-thickening-fluid-body-armor-market ,
  https://www.futuremarketinsights.com/reports/liquid-armor-materials-market
- *The honest caveat on fit:* in fielded armour the STF is impregnated into a fabric and
  is **confined** â€” there is no free surface during the ballistic event. The solver's
  capillary interface is load-bearing for (a) laboratory rod-into-open-bath impact, which
  is where the physics is validated, and (b) the impregnation/coating manufacturing step.
  For the actual ballistic problem it is not the governing feature.
- **Physics the solver would need:**
  1. **Coupling to a fibre/fabric network** (porous-medium drag, or an FSI coupling to a
     textile model). Without it you are modelling the bath experiment, not the product.
  2. **3-D**, and high strain rate: real threats are non-axisymmetric and involve
     compressibility of the entrained air and of the packing.
  3. **Fracture / failure of the jammed state.** Slutzky et al. show the jammed solid
     cracks. A purely viscous closure will keep thickening where the material actually
     breaks, so the solver will over-predict protection. This is a credibility risk, not
     just a missing feature.
  4. Reynolds dilatancy + Darcy backflow (per Athani et al., PRF 2025) to set the stress
     level correctly in a locally jammed region.

### Rank 4 â€” Injectable/dispensed high-concentration formulations and small-bore clogging

**Fit: good (axisymmetric, free surface at the orifice, jamming is the failure). Pain:
high value per unit, modest total volume.**

- Veryst Engineering's own COMSOL case study names exactly this cluster: "the blockage of
  inkjet printer nozzles in additive manufacturing processes, the clogging in medical
  injection devices by protein aggregates, and the occlusion of pores by organic
  contaminants in water filtration membranes." [source states â€”
  https://www.veryst.com/case-studies/shear-jamming-dense-suspensions] That is a
  consultancy telling you which industries pay for this analysis.
- *What goes wrong:* high-concentration biologic formulations and dense particulate
  dispersions jam in needles, nozzles and print heads; the failure is intermittent and
  therefore very expensive to characterise experimentally.
- **Physics the solver would need:** Brownian/colloidal stability and aggregation kinetics
  (protein aggregation is not frictional contact â€” the analogy is loose and a referee from
  that industry will say so); wall slip; migration; and a constriction/geometry sweep
  capability. **Honest assessment: the mechanism is only partly the same physics.** Silica
  and cornstarch map well; protein aggregates do not.

### Explicitly ranked LOW, with the reason (state these to pre-empt the question)

- **Drilling fluids and proppant transport / screen-out.** The flow is fully confined
  inside a fracture or wellbore â€” **no free surface** â€” and the failure mechanism is
  proppant *bridging* across a slot, a granular geometric-jamming problem, not
  stress-induced frictional shear jamming. Near-wellbore friction is cited as responsible
  for "90-95% of screenouts". [source states â€” SPE/industry sources via
  https://onepetro.org/SPEHFTC/proceedings-abstract/22HFTC/22HFTC/D021S004R005/476943 and
  https://shalemag.com/detecting-avoiding-screenouts-hydraulic-fracturing/ ]. Poor fit;
  would need a full two-phase granular (mu(J)) model and 3-D.
- **STF dampers / adaptive damping devices.** Confined, no free surface, and the design
  question is a lumped force-displacement law rather than a field solution. The solver
  would be over-specified for it.
- **Food processing.** Genuine free surfaces (depositing, filling, coating) and genuine
  starch/chocolate jamming, but I found no quantified simulation-driven pain point in this
  search, and R&D simulation budgets in food processing are small. Plausible but
  unevidenced â€” do not lead with it.


---

## TASK 4 â€” FUNDING ROUTES (US)

Named programs only, with sources. I did not verify current deadlines or budgets â€” check
each before writing anything.

### NSF â€” the natural academic home

- **CBET / Transport Phenomena cluster â€” Particulate and Multiphase Processes, PD 23-1415.**
  Supports "physico-chemical phenomena that govern particulate and multiphase systems,
  including flow of suspensions, drops and bubbles, granular and granular-fluid flows,
  behavior of micro- and nanostructured fluids", and states "Collaborative and
  interdisciplinary proposals are encouraged, especially those that involve a combination
  of experiment with theory and/or modeling." [source states â€”
  https://www.nsf.gov/funding/opportunities/particulate-multiphase-processes ]
  **This is the single best-matched program.** Dense suspensions are explicitly named in
  NSF's own routing guidance.
- **CBET / Transport Phenomena â€” Fluid Dynamics, PD 23-1443** (supersedes PD 17-1443).
  Emphasis on "the behavior of new fluid materials and innovative uses of fluids in
  manufacturing". The free-surface/impact framing fits here; the suspension framing fits
  1415 better. Do not submit to both. [source: NSF program pages via
  https://www.nsf.gov/funding/opportunities/fluid-dynamics/505698 ]
- **DMR â€” Condensed Matter Physics / Polymers**, if the pitch is the jamming transition as
  physics rather than as engineering. Weaker fit for a solver.
- **CMMI â€” Advanced Manufacturing**, if the pitch is direct ink writing / paste extrusion.
  NSF's manufacturing focus area explicitly covers additive and digital manufacturing and
  process modeling. [https://www.nsf.gov/focus-areas/manufacturing]
  Note that the Annual Review that best states the DIW gap (Sauret, Ray & Compton 2026) is
  itself indexed in the **NSF Public Access Repository**, i.e. NSF already funds this
  problem statement.
- **NSF SBIR/STTR (America's Seed Fund)** if a company is the vehicle. STTR is the right
  instrument for a university-origin solver.

### DoD â€” the strongest fit for impact protection

- **Army Research Office (ARO) Core BAA for fundamental research**, Mechanical Sciences
  and Materials Science divisions. The current core BAA (W911NF-23-S-0001) is open with
  rolling submission through Nov 2027. [https://arl.devcom.army.mil/collaborate-with-us/opportunity/arl-baa/]
  ARO's stated disciplines include mechanical sciences and materials science.
- **DEVCOM Soldier Center BAA for Basic and Applied Research (W911QY-25-R-0023)**, Soldier
  Protection and Survivability topics; concept papers accepted on a rolling basis until
  28 Feb 2030. [source states, via
  https://defencescienceinstitute.com/funding-opportunity/us-devcom-soldier-center-baa-for-basic-and-applied-research-w911qy-25-r-0023/ ]
  A rolling-basis concept paper is the cheapest possible shot on goal here.
- **ONR / AFOSR** have complex-fluids and materials portfolios but I did not verify a
  matching program element. **UNVERIFIED.**
- **DoD SBIR/STTR** â€” the Army runs open topics (e.g. a "Battery Focused Open Topic",
  https://armysbir.army.mil/topics/battery-focused-open-topic/ ), which are a route in
  without waiting for a bespoke solicitation.

### DOE

- **DOE SBIR/STTR** â€” management was consolidated under the Office of Technology
  Commercialization; FY26 Phase I awards are up to ~$200k for ~6 months, with Phase II
  follow-on. [https://www.energy.gov/technologycommercialization/doe-small-business-innovation-research-sbir-and-small-business]
  I could **not** find a 2026 topic that names slurry/suspension modelling specifically.
  **UNVERIFIED â€” the topic list must be checked when the annual solicitation opens.**
- **Vehicle Technologies Office / battery manufacturing R&D** is the obvious home for the
  electrode-slurry coating application, but I did not verify an open, matching FOA.
  **UNVERIFIED.**
- **Office of Science / BES** funds dense-suspension physics but is a poor fit for an
  applied solver.

### Industry consortia and direct routes

- **Consultancies as channel partners rather than competitors.** Veryst Engineering
  already sells shear-jamming analysis into inkjet/AM, medical injection and filtration
  (see Â§2.3). A firm that has hand-built this model in COMSOL is a plausible partner or
  acquirer of the physics, and a source of named end customers.
- I did **not** find a US consortium specifically for dense-suspension processing. The
  nearest analogues are the coating/converting community (AIMCAL) and battery-manufacturing
  consortia. **UNVERIFIED â€” do not assert a consortium exists without checking.**

---

## EXECUTIVE VERDICT

**Competitive.** The ground is **partly occupied**, and more so than the framing assumed.
A frictional-contact-fraction field with its own evolution equation inside a continuum
solver already exists twice: Baumgarten & Kamrin (PNAS 2019, MPM) and Angerman, Sandnes,
Seto & Ellero (Phys. Fluids 2024 / JFM 2025, SPH, explicitly Wyart-Cates and explicitly
advected). Athani, Metzger, Forterre & Mari (PRF 2025) have built the Wyart-Cates +
Reynolds-dilatancy closure and named impact as its target. Mobaseri, Kumar & Cheng (PNAS
2025) already run cornstarch drop impact on Basilisk with VOF. What **nobody** has is the
intersection: a Wyart-Cates contact fraction that is advected **and** relaxes on an
absolute timescale, inside a **sharp-interface, surface-tension-resolved** free-surface
flow. Every competitor holds at most two of those three, and two of them say so in print.
The `1/t_r` term is the sharpest technical differentiator â€” both prior implementations
relax f on a pure strain clock, so neither can un-jam a region that has stopped moving.

**Commercial CFD.** Nothing off the shelf does this. Fluent, STAR-CCM+ and OpenFOAM offer
only algebraic eta(gammadot) (power law, Carreau/Carreau-Yasuda, Cross, Casson,
Herschel-Bulkley, Bingham). COMSOL is the exception worth respecting: its Polymer Flow
Module carries seven differential viscoelastic models and, uniquely, the **Houska
thixotropy** model with a structure field obeying its own ODE â€” but thixotropy has the
wrong sign and no jamming point, so it cannot represent shear jamming. rheoTool adds
viscoelasticity and a two-phase solver to OpenFOAM, still with no jamming model. The one
commercial-adjacent Wyart-Cates implementation found is a bespoke **Veryst Engineering**
COMSOL build that tracks "the fraction of frictional contacts" â€” proof of demand, and a
reason not to claim commercial CFD has never touched this.

**Top 4 applications (fit x pain), with the physics gap for each.**
1. **Direct ink writing / dense-paste extrusion AM.** Axisymmetric nozzle, capillary
   free filament, clogging-then-unclogging is literally the failure mode; the field's own
   review (Annu. Rev. Fluid Mech. 58, 2026) calls for replacing "trial-and-error
   approaches with predictive design methodologies". *Needs:* wall slip (highest value),
   particle migration, yield stress and jamming together, eventually 3-D.
2. **High-solids coating and atomization, battery electrode slurries first.** Nozzle shear
   rates reach 10^5 s^-1, well past DST onset; thickening drives ribbing and air
   entrainment. *Needs:* **3-D** (ribbing and ligaments are non-axisymmetric â€” the biggest
   lift in this document), moving contact line, drying/non-isothermal, migration.
3. **STF impact protection.** Best physics match, and current engineering practice is to
   substitute "solid-equivalent surrogate models" for the fluid. *Needs:* coupling to the
   fibre/fabric network, 3-D, fracture of the jammed state (Slutzky et al., PRE 2026, show
   it cracks â€” a viscous closure will over-predict protection), Reynolds dilatancy + Darcy
   backflow.
4. **Small-bore clogging: injectables, inkjet, dispensing.** Veryst names these as paying
   customers. *Needs:* wall slip, migration, and honesty that protein aggregation is not
   frictional contact.
   *Poor fits, say so first:* proppant screen-out (confined, no free surface, granular
   bridging) and STF dampers (confined, lumped force law).

**Funding.** NSF **CBET Particulate and Multiphase Processes (PD 23-1415)** is the
best-matched program and names suspensions explicitly; **CBET Fluid Dynamics (PD 23-1443)**
is the alternative framing; **CMMI Advanced Manufacturing** if pitched as DIW. On the
defence side, the **DEVCOM Soldier Center BAA (W911QY-25-R-0023)** takes concept papers on
a rolling basis to Feb 2030 and the **ARO core BAA (W911NF-23-S-0001)** is open to Nov
2027 â€” both are cheap shots on goal. DOE SBIR/STTR is plausible but I found no matching
2026 topic and could not verify one. Veryst Engineering is a credible channel partner
rather than only a competitor.


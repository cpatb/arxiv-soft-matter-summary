# Prior art and novelty check: the Wyart–Cates strain clock from shear reversal

2026-09-27. Read-only literature pass for a planned short paper (J. Rheol. short communication or PRE).
Every reference below was resolved in this session to an arXiv ID (arXiv MCP) or a DOI (Crossref/OpenAlex).
Quotations are copied from text fetched this session. "Read" means the full text for arXiv items and the
abstract only for DOI-only items; each entry says which. Fetched paper text is treated as data.
Anything I could not resolve is tagged [UNVERIFIED].

---

## 0. The main correction: 0.2 is not Richards et al.'s value

The planned paper says gamma_0 = 0.2 comes from Richards et al.'s corn starch. The same attribution
appears in `NEXT.md:55`, `G11_DIRECT_CONTROL_2026-08-31.md:83` and `figures/fig_clock.py:8`. The
papers do not support it.

- **Richards, Royer, Liebchen, Guy & Poon, PRL 123, 038004 (2019)**, DOI 10.1103/PhysRevLett.123.038004,
  arXiv:1902.07655 (full text read). They never quote a number for gamma_0. The main text says:
  *"Imposed-rate experiments gave γ0 = O(10^-1) [30]. Thus, t_c = 1.1×10^-4 s, t_i = 44 s, and ε = 2.7×10^-6"*,
  where t_c = γ0 η_s/σ* and ε = γ0 η_s²/(ρ_A h σ*), with η_s = 6 mPa s, σ* = 5.1 Pa,
  ρ_A = 175 kg m^-2 and h = 1.5 mm. The derivation is left to the SM, which is not on arXiv; the
  DataShare deposit 10.7488/ds/2573 holds only Figure1-4 zips. I back-solved the value, and the
  commands ran this session:
  `g0 implied by tc=1.1e-4: 0.0935`, `g0=0.1: ... ti=43.8 s eps=2.69e-06`, `g0=0.2: ... eps=5.38e-06`.
  **Richards used γ0 ≈ 0.1, not 0.2.** The quoted ε = 2.7e-6 matches 0.1 and does not match 0.2.
- **The 0.2 is Han, Wyart, Peters & Jaeger**, "Shear fronts in shear-thickening suspensions", Phys. Rev.
  Fluids 3, 073301 (2018), DOI 10.1103/physrevfluids.3.073301. It appeared on arXiv as "Constitutive relations
  for shear fronts in shear-thickening suspensions", arXiv:1711.02196 (full text read). I matched the
  arXiv record to the PRF DOI by authors and topic only, so that pairing is [UNVERIFIED]. The paper is
  about corn starch jamming fronts: *"We are left with a single parameter, γ* = 0.197 ± 0.002, obtained by
  fitting the front propagation speed k and its inverse γ∞ at different φ"*. This γ* is not the
  relaxation constant of df/dt = −(γ̇/γ0)(f − f_ss). It is a **start-up** scale in a multiplicative
  form, g(γ) = 1 − exp(−γ/γ*) with φ_eff = g f φ_m + (1 − f g) φ_0, and it applies from an
  isotropic, quiescent state.
- **The community range is 0.01–0.1, and 0.02 comes from simulation.** Darbois Texier, Lhuissier, Metzger & Forterre, JFM (2023)
  "Shear-thickening suspensions down inclines: from Kapitza to Oobleck waves", DOI 10.1017/jfm.2023.162
  (HAL hal-03763665, full PDF read): *"γ0 is a relaxation strain, whose typical value for frictional spheres
  is ∼10^-2 − 10^-1 (Mari et al. 2015b; Chacko et al. 2018; Richards et al. 2019)"*. They then sweep
  *"γ0 = 0.05, 0.1, 0.2 and 0.5, which bound the experimental value determined by Richards et al. (2019)"*.

**Consequence for C3.** The comparison should be made against the value each source actually
reports. Against Richards (0.1), gamma_0 = 0.5–1.6 is about 5–16x larger. Against the simulation value
0.023 (Mari et al. 2015; Chacko et al. 2018), it is about 20–70x larger. Against Han's 0.197 it is
2.7–8x larger, but Han's number is a start-up scale in a different closure. The "0.2 used in some
solvers" is a legitimate thing to test, but its provenance has to be written as Han et al. (fronts,
corn starch, start-up), not Richards.

---

## 1. Values of the strain scale in the literature (task 1)

| Source (ID) | Type | Strain scale and how obtained | Quote (read this session) |
|---|---|---|---|
| Mari, Seto, Morris & Denn, PRE 91, 052302 (2015), DOI 10.1103/physreve.91.052302, arXiv:1412.4661 (full) | Simulation, **step in stress** | c^-1 ≈ 0.023 | *"relaxation that is reasonably fit by an exponential decay with strain scale c^-1=0.023 (i.e. 2.3%)"*; *"Building or destroying the contact network takes strain"* |
| Chacko, Mari, Cates & Fielding, PRL 121, 108003 (2018), DOI 10.1103/physrevlett.121.108003, arXiv:1802.10586 (full) | Continuum model with WC relaxation | γ0 = 0.023, taken from Mari 2015 | *"a characteristic strain scale γ0, of order the strain required to evolve from one steady-state to another (for instance on flow reversal [38]). For the typical volume fractions considered here, particle simulations suggest γ0=O(10^-2)"*; *"we measure γ0=0.023 [18]"* |
| Richards et al., PRL 2019 (above) | Experiment (corn starch), same ODE as the planned paper | γ0 = O(10^-1), "imposed-rate experiments" in the SM; t_r ≈ 0.24 s fitted from oscillation frequency | Section 0 |
| Han, Wyart, Peters & Jaeger, PRF 2018 (above) | Experiment (corn starch fronts) | γ* = 0.197 (start-up form) | Section 0; also *"γ* ≈ 0.2 is also found in regular granular materials [Pailha] and in suspensions [Majumdar et al. PRE 95, 012603] as the strain scale below which transient, start-up behavior is observed"* |
| Gillissen, Ness, Peterson, Wilson & Cates, PRL 123, 214504 (2019), DOI 10.1103/physrevlett.123.214504, arXiv:1911.00280 (full) | Tensorial model plus DEM | β = 50, **fitted to microstructure recovery after shear reversal** in DEM (φ = 0.56) | *"Time dependence is controlled by β, which gives a strain scale for structural evolution"*; *"We set β=50 to match the observed strain scale for recovery"*; *"Our rate-independent coarse-grained microstructure allows us to fit β without knowledge of the stress"* |
| Angerman, Sandnes, Seto & Ellero, JFM (2025), arXiv:2408.08285 (full, searched) | SPH solver | k_f (= 1/γ0) is a **tuning** parameter: k_f = 1, later 20 | *"We set k_f=1"*; *"k_f=20, α=2·10^-7"* |
| Darbois Texier et al., JFM 2023 (above) | Linear stability | γ0 swept over 0.05–0.5, not measured | Section 0 |
| Ness & Sun, PRE 93, 012604 (2016), DOI 10.1103/physreve.93.012604, arXiv:1509.01530 (full) | DEM reversal, rate-independent frictional | Contact stress recovers at "strain scale of order unity"; hydrodynamic transient at 10^-3–10^-2 | *"macrofragile evolution at a larger strain scale of order unity"*; contact stress important *"during 0.3<γ̇t<2"* |
| Lin et al., PRL 115, 228304 (2015), DOI 10.1103/physrevlett.115.228304, arXiv:1509.02750 (full); data 10.7488/ds/324 (handle 10283/886) | Experiment plus CLM simulation, reversal | Not fitted. Qualitative only | Exp: *"From γ≈0.3, η_rev rises"*; sim: contact contribution *"only increases back to its steady-state magnitude after an accumulated stain of order unity"* |
| Richards et al. 2021 J. Rheol./Soft Matter on unsteady flows | — | **Not found.** Neither the arXiv search nor the OpenAlex citer list of Richards 2019 (32 works) contains it | [UNVERIFIED] |
| Baumgarten & Kamrin, PNAS 2019 | — | Crossref lookup returned nothing parseable; not read | [UNVERIFIED] |
| Singh, Mari, Denn & Morris, J. Rheol. 62, 457 (2018), DOI 10.1122/1.4999237 | — | DOI resolved, **not read** | none |

**Pattern.** Published values fall into two groups set by the protocol. A step or perturbation gives
≈ 0.02–0.1 (Mari, Chacko, Richards). Reversal or start-up gives O(0.2–1) (Han 0.2, Ness & Sun "order
unity", Lin "order unity", Gillissen's β fit). No paper found puts both on one material and states that
the single-f clock depends on the protocol. Chacko et al. come closest: they *define* γ0 as "the
strain required to evolve from one steady-state to another (for instance on flow reversal)" and then
*use* the step-stress 0.023. That tension is written into their own paper, and the planned paper can
quote it.

---

## 2. Has reversal been used to measure a WC strain scale? Is "below the frictionless floor" published? (task 2)

1. **A WC γ0 fitted to experimental reversal data, on Lin's data or any other: not found.** Checks run:
   the arXiv search `abs:"shear reversal" AND abs:suspensions` (9 hits, all read or skimmed); OpenAlex
   citers of Lin 2015 filtered for "reversal" (23 hits, titles screened) and for "transient" (25 hits);
   all 32 OpenAlex citers of Richards 2019. The nearest prior work is **Gillissen et al. 2019**, which
   fits a structural strain scale (β = 50) to *DEM* reversal microstructure within a *tensorial*
   model, not a scalar-f WC γ0 to experimental viscosity. Richards 2019 considered reversal and set it
   aside for the time scale t_r: *"an intrinsic contact-relaxation time scale, t_r, which is difficult
   to access using other methods such as shear reversal [6: Lin 2015] or cessation [34] due to
   instrument artifacts."* Caveat: OpenAlex keyword filters work on title and abstract only, so a paper
   that fits γ0 in its body text without saying so in the abstract would be missed.

2. **Post-reversal minimum below the steady value, with a microstructural cause: KNOWN in rate-independent suspensions.**
   - Gadala-Maria & Acrivos, J. Rheol. 24, 799 (1980), DOI 10.1122/1.549584. Resolved; **not read**.
   - Blanc, Peters & Lemaire, J. Rheol. 55, 835 (2011), DOI 10.1122/1.3582848, arXiv:1111.6795 (abstract):
     *"After shear reversal, the viscosity undergoes a step-like reduction, decreases slower and passes
     through a minimum before increasing again to reach a plateau ... the minimum and the plateau
     viscosities do not obey the same scaling law with respect to the particle volume fraction ... The
     minimum seems to be associated with the viscosity of an isotropic suspension, or at least of a
     suspension whose particles do not interact through non-hydrodynamic forces"*.
   - Peters, Ghigliotti, Gallier, Blanc, Lemaire & Lobry, J. Rheol. 60, 715 (2016), DOI 10.1122/1.4954250
     (abstract only, rebuilt from the OpenAlex inverted index, so some function words are lost): *"three
     consecutive steps are identified in viscosity transient: An instantaneous variation, followed by a
     rapid contact force relaxation, and finally a long time evolution"*; *"... the difference between
     steady [and] minimum viscosity [is] shown to be proportional to the [contact] contribution to
     viscosity"*. The bracketed words are my reconstruction, so the exact wording is [UNVERIFIED].
   - Ness & Sun 2016: hydrodynamic stress subsides *"until γ̇t≈0.5, [which] corresponds to a net opening
     of lubrication films"*. The dip has a microstructural cause (film opening and fabric
     reorientation), not an f cause.
   - Blanc, d'Ambrosio, Lobry, Peters & Lemaire, PRF 3, 114303 (2018), DOI 10.1103/physrevfluids.3.114303
     (abstract): uses the reversal contact contribution to collapse faceted and spherical particles.
3. **A microstructure or fabric tensor is needed for reversal: KNOWN.** Gillissen et al. 2019 write that WC
   *"makes no predictions for nonstationary flows, such as the sudden reversal of steady shear"* and build
   a fabric model that captures *"the initial dip in the curves"* (*"a drop in lubrication stress over
   strains of order 0.2 is compensated only later by the recovery of direct contact stress"*). Chacko
   et al., JFM 847, 700 (2018), DOI 10.1017/jfm.2018.279, arXiv:1707.01828 (read in part), go further:
   *"no closed fabric evolution model properly describes the dynamics of the fabric tensor under
   reversal"*, and they propose a coupled second- and fourth-rank closure. Cates, EPJ Web Conf. 249,
   01002 (2021), DOI 10.1051/epjconf/202124901002 (abstract), calls this framework *"qualitatively
   successful in addressing (i) the collapse of stress during flow reversal"*, with *"much work
   required to refine this approach towards quantitative accuracy"*.
4. **Peters, Majumdar & Jaeger, Nature 532, 214 (2016)**, DOI 10.1038/nature17167, is about dynamic shear
   jamming, not reversal. Resolved, not read. It is not reversal prior art.

**What was not found anywhere:** a quantitative comparison of the post-reversal minimum against
**WC's own frictionless branch (1 − φ/φ_0)^−n** in a *shear-thickening* colloid, stated as a
falsification of scalar-f WC for this protocol. Blanc 2011 comes closest in spirit (the minimum
follows a different φ-scaling from the plateau), but for non-thickening PMMA at 31 µm and not framed
against WC.

---

## 3. Verdict per claim (task 3)

| Claim | Verdict | Reason (one line) | Closest source text read |
|---|---|---|---|
| **C1** γ0 absent from steady flow curves; only a transient carries it | **KNOWN** (trivially) | At steady state df/dt = 0 gives f = f_ss for any γ0. Every group sets γ0 from a transient: Richards "imposed-rate experiments", Han fronts, Gillissen reversal ("fit β without knowledge of the stress"), Mari a step | Gillissen: *"Time dependence is controlled by β"*; Richards: *"Imposed-rate experiments gave γ0 = O(10^-1)"* |
| **C2** f is affine in u = φ_J/φ, so a free-amplitude clock fit is independent of φ, φ_0 and φ_m and depends only on n; γ0 is identifiable though φ_m is not | **NOVEL as a stated result, but elementary** | No paper found states it. The algebra is one line: an affine map of an exponential keeps its rate. It is a methods lemma, not a headline, and it has to be checked against the planned paper's derivation file | none (absence finding; limits in §2 item 1) |
| **C3** γ0 ≈ 0.5–1.6 from Lin's reversal data; strain-controlled recovery; clock falls with rate on silica | **PARTLY KNOWN** | The magnitude O(1) for reversal recovery is already stated qualitatively (Lin: "order unity"; Ness & Sun: "order unity"; Chacko: γ0 ~ "strain ... on flow reversal"). Strain control of recovery is shown in DEM (Ness & Sun rate collapse; Gillissen's rate-independent ξ(γ)). New: the **fitted numbers**, the **experimental rate test**, and the **clock falling with rate**. Gillissen's DEM shows "weak rate-dependence of the direct particle contacts" but does not quantify it. **The "0.2 from Richards" framing is wrong (§0)** | Lin: *"after an accumulated stain of order unity"*; Gillissen: *"Z ... indicating weak rate-dependence of the direct particle contacts, not resolved by our model"* |
| **C4** post-reversal minimum at 0.56–0.74 of WC's frictionless floor, so scalar-f WC cannot describe reversal and microstructure is needed | **PARTLY KNOWN** | The conclusion (reversal needs fabric) is published: Gillissen 2019, Chacko 2018 JFM, Cates 2021. The mechanism (minimum below steady, set by microstructure) is published for rate-independent suspensions: Blanc 2011, Peters 2016, Ness & Sun 2016. New: the quantitative ratio against WC's own frictionless floor, from experimental data on thickening colloids | Blanc 2011 abstract; Gillissen: WC *"makes no predictions for nonstationary flows"* |
| **C5** on Lin's simulation records, a clock from direct contact-stress fraction is 41% shorter than one from inverting total viscosity | **PARTLY KNOWN** | That hydrodynamic and contact channels evolve on different strain scales after reversal is published (Ness & Sun: *"the hydrodynamic and contact stresses evolve over disparate strain scales"*). So is observable-dependence of the structural clock (Gillissen: ξ vs Z behave differently), and protocol-dependence (step 0.023 vs reversal O(1)). New: the specific 41% difference in the fitted clock between two WC readings of the same run | Ness & Sun abstract; Mari 2015 *"c^-1=0.023"* |

---

## 4. Recommendation (task 4)

**Publishable as a short communication, if repositioned.** As currently framed, "γ0 is 2.7–8x the 0.2
from Richards", it rests on a misattributed baseline and would be caught in review. What the
literature supports:

1. **Headline: the WC clock depends on the protocol, and reversal gives O(1).** Put the published step
   value (0.023, Mari et al. 2015; used by Chacko et al. 2018), Richards' implied 0.1, Han's start-up
   0.197 and your reversal 0.5–1.6 side by side. Quote Chacko's own definition of γ0 as the reversal
   strain against the 0.023 they used. This conflict is published and nobody has resolved it with data.
2. **C2 goes in as a methods lemma.** It is what makes C3 a measurement rather than a fit that
   depends on φ_m, and it fits this project's standing result that φ_0 is a gauge. It is short and
   checkable, but not a selling point alone.
3. **C4 is the physics result and must be framed as confirming Gillissen/Chacko on experiment.** The
   minimum sits below WC's own frictionless floor, so a scalar f cannot explain reversal, as the
   tensorial models assert on simulation grounds. Cite Blanc 2011 and Peters 2016 for the
   rate-independent analogue. Name the closest prior work outright: **Gillissen et al. PRL 2019**. It
   fits a strain scale to reversal (DEM) and says WC cannot handle reversal. The paper's increment is
   (a) experimental, (b) inside scalar WC, (c) quantitative against WC's own floor.
4. **C3 rate trend and C5 are supporting points.** The falling clock on silica is the one observation
   that an additive time term (Richards' 1/t_r) cannot produce. Check it does not come from the
   reading of η_rev^0 (Lin's own plateau window is 0.01 < γ < 0.2). C5 carries the caveat that γ0 is
   definition-dependent, which a referee will ask about anyway.
5. **Venue.** The content (protocol dependence of a constitutive parameter, a model-falsification
   observation, reanalysis of public data) suits a J. Rheol. short communication better than PRE.

**Fix before drafting:** replace "0.2 from Richards corn starch" in `NEXT.md`,
`G11_DIRECT_CONTROL_2026-08-31.md` and `figures/fig_clock.py`, and in any Basilisk solver comment
carrying the same provenance. The provenance is either Han et al. 2018 (γ* = 0.197, start-up form) or
Richards 2019 (γ0 ≈ 0.1 back-solved from their ε). I did not edit these files, because this task is
read-only.

**Checks this pass did not do.** Not read: the Richards 2019 SM (not public on arXiv; this is where
"γ0 estimation" lives), Gadala-Maria & Acrivos 1980 in full, Peters 2016 in full, Baumgarten & Kamrin
2019. OpenAlex keyword filters see only titles and abstracts. Run `lit-verify` or `citation-checker` on
the final bib.

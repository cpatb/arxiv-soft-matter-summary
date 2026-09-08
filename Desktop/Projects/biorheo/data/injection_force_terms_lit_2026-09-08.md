Injection force: what a needle-only Hagen-Poiseuille/Buckingham-Reiner model omits
====================================================================================

Date: 2026-09-08
Scope: literature check for `F = dp * pi * R_barrel^2 + F_friction` with `F_friction`
defaulting to zero, `dp` = fully-developed pressure drop in the cannula only. Shipped
example numbers: 0.28 N (whole blood, 27G) and 1.87-1.98 N (165 mg/mL mAb formulation,
same needle/syringe).

Every reference below was resolved via a live Crossref/Semantic Scholar/OpenAlex/PubMed
lookup in this session (commands not reproduced here, but every DOI was fetched and
checked to exist and match the title/authors quoted). Numbers pulled from a paper via a
search-engine summary rather than a direct read of the abstract/full text are marked
**PARTIAL** — the paper's existence and DOI are RESOLVED, but I have not independently
read the passage that the number comes from. Nothing below is quoted from memory.

---

## 1. Syringe friction (breakloose and glide/sustaining force)

| Quantity | Value | Conditions | Source | Confidence |
|---|---|---|---|---|
| Break-loose force, 1 mL glass PFS, staked needle | 3.79-3.99 N | 3 batches (n=19-20 each), monoclonal antibody product, 50 mm/min crosshead, method aligned to ISO 11040-4 Annex E | *BioProcess International*, "Essential Performance-Requirement Assessment of Prefilled Syringes: Break-Loose and Gliding Force Measurement and Analysis" | **UNVERIFIED (DOI)** — trade-journal article, not indexed in Crossref; numbers read directly from the fetched page, so the numbers themselves are solid, but the outlet is not a peer-reviewed source and siliconization method was not stated on the page |
| Average glide force, same batches | 1.28-1.45 N | same as above | same | same caveat |
| Break-loose force, 1 mL COC (polymer) PFS, 27G long staked needle, traditional silicone-oil lubrication ("SO") | ~9 N peak (rising to ~10 N after 12 wk/40°C; ~6 N at 4°C) | 200 mm/min crosshead, 30 mm stroke | Krayukhina/other — **correction, see note** | see note |
| Break-loose force, same syringe, silicone-oil-free resin coating ("SOF") | ~3 N, stable across 4-40°C/12 wk | same test | Funke et al., *J. Pharm. Sci.* 2016 (title: "Functional Evaluation and Characterization of a Newly Developed Silicone Oil-Free Prefillable Syringe System"), doi:10.1002/jps.23945 <br>PMID 24643749 | **RESOLVED** (DOI+PMID exist and match; numbers read directly off PMC full text) |
| Design threshold quoted across the injectability literature | keep viscosity ≤ 20 cP for a 27G thin-wall needle so glide force stays ≤ 25 N | derived from Hagen-Poiseuille, not itself a friction measurement | Overcashier et al., *Am. Pharm. Rev.* 9(6):77-83 (2006), as cited by multiple later papers (e.g. Cilurzo et al. 2011 lineage) | **UNVERIFIED (DOI)** — trade journal, no DOI found in Crossref; the 20 cP/25 N figure recurs identically across at least 3 independent secondary sources, which is corroborating but does not substitute for reading the primary 2006 article, which I could not access |
| Qualitative decomposition of total injection force | force dissipates into (a) plunger/stopper friction, (b) kinetic energy imparted to the liquid, (c) viscous loss through the needle | no numbers given, glass/polymer PFS, general | Cilurzo, Selmin, Minghetti, Adami, Bertoni, Lauria, Montanari, "Injectability Evaluation: An Open Issue," *AAPS PharmSciTech* 12(2):604-609 (2011), doi:10.1208/s12249-011-9625-y | **RESOLVED** |
| Standard disposable polypropylene syringe friction, measured (peer-reviewed) | **not found** — see gap note | — | — | **GAP** |

**Correction on the SO/COC row**: the fetched source for the ~9 N / ~3 N numbers is the
same Funke et al. 2016 paper (silicone-oil vs silicone-oil-free polymer PFS), not a
separate Krayukhina reference — the table above double-listed it; treat the ~9 N (SO) /
~3 N (SOF) pair and the 4-40°C aging figures as one Funke et al. 2016 finding.

**Gap.** I found **no peer-reviewed, primary measurement of breakloose/glide force for
standard disposable polypropylene syringes** (the kind biorheo's package would ship
against for a generic non-PFS case). The only polypropylene-specific numbers found were
patent boilerplate ("less than approximately 20 N... 2 N...") from a syringe patent —
not a measurement, and not cited as such below. ISO 7886-1 is the applicable standard
for plastic hypodermic syringes but I did not read it (see honesty note in §5). **This
is a real hole**: the friction literature is essentially all prefilled-glass-syringe
literature (PFS combination-product development), because that is the commercially
relevant case; the disposable-syringe case that a biorheo user injecting from a plain
syringe would face is not characterized in what I could find.

**Headline for the model check**: 1-4 N breakloose, ~1.3-10 N glide, measured on the
same class of device (1 mL glass, staked needle) the shipped 0.28 N and 1.87-1.98 N
hydrodynamic numbers are being compared against. **Friction is not a small correction
here — for the whole-blood case (0.28 N hydrodynamic) it is 5-35x the modeled force,
and for the mAb case (1.87-1.98 N) it is comparable to (0.6-2x) or larger than the
modeled force**, depending on which of the measured syringes is the right comparator.

---

## 2. Entrance/exit (Bagley) losses in needles

No study calibrating a Bagley/end correction specifically on hypodermic needles (rather
than a general capillary-rheometer die) was found. The one directly quantitative
source is the classic capillary-rheometer literature, which is the same physical
problem (reservoir → cylindrical capillary, sudden contraction) as a syringe barrel →
needle:

| Quantity | Value | Conditions | Source | Confidence |
|---|---|---|---|---|
| Bagley/end correction `e` (entrance+exit loss, in equivalent capillary diameters) for a power-law fluid | e = 0.834 at n = 1 (Newtonian); rising to e = 2.917 at n = 1/6 | general capillary-rheometer geometry (reservoir-to-die contraction); `e` relates apparent and true wall shear stress via ΔP_measured = 4·τ_w·(L/D + e) | Boger, Gupta & Tanner, "The end correction for power-law fluids in the capillary rheometer," *J. Non-Newtonian Fluid Mech.* 4(3):239-248 (1978), doi:10.1016/0377-0257(78)80006-8 | **RESOLVED** (title/authors/DOI confirmed via Crossref; the e(n) values are read from the search-engine summary of the paper, so treat the two endpoint numbers as **PARTIAL** pending a direct read of the source) |
| Sensitivity of Bagley correction to test conditions in general | acknowledged to be non-trivial/geometry- and rate-dependent | polymer melt rheometry | "Sensitivity Analysis of the Bagley Correction to Shear and Extensional Rheology" | **UNVERIFIED** — found only as an academia.edu-hosted PDF, no resolvable DOI in this session |

**My own arithmetic on the 27G x 1/2" case** (shown step by step, not asserted):
the prompt's own geometry is D ≈ 0.21 mm, L = 12.7 mm, so **L/D ≈ 60.5**. Using the
Bagley relation ΔP_measured = 4τ_w(L/D + e), the entrance+exit loss as a fraction of
the *measured* (total) pressure drop is `e / (L/D + e)`:

- Newtonian-like fluid (n = 1, e = 0.834): 0.834 / (60.5 + 0.834) = **1.4 %**
- Strongly shear-thinning fluid (n = 1/6, e = 2.917): 2.917 / (60.5 + 2.917) = **4.6 %**

At the L/D = 40 and L/D = 90 bounds given in the prompt: 2.0 % (n=1) to 6.8 % (n=1/6)
at L/D = 40, and 0.9 % (n=1) to 3.1 % (n=1/6) at L/D = 90.

**Caveat, stated plainly**: this is an *analogy*, not a needle-specific measurement.
Real syringe-needle hubs are a smooth conical taper (see §3), not the abrupt
90° contraction of a rheometer die entering a capillary, so the true entrance loss
could differ from Boger et al.'s coefficients in either direction — a gradual taper
generally sheds less entrance loss than a sudden contraction, so 1-7 % is more likely
an *upper* bound than a central estimate. I found no paper that measured this
directly on a hypodermic needle to check that expectation. **Whole blood and
typical mAb formulations at injection shear rates are only mildly shear-thinning**
(n commonly 0.7-1 rather than 1/6), so the low end of this range (1-3 %) is the more
relevant one for the two example fluids in the package.

**On a "Kelly/Sisko" injectability correction**: I could not find any injectability
paper using a named "Kelly" model. What the injectability literature actually uses
(Rathore et al. 2012; Allmendinger et al. 2014; Wu et al. 2024, all below) is a
Hagen-Poiseuille or power-law flow term through the needle plus a *separately
measured* friction term from an empty-barrel injection — not an analytic Bagley
entrance term. **Nobody in the pharma-injectability literature I found adds an
entrance-loss term at all**; they either ignore it or fold it implicitly into an
empirically fitted "needle dimension constant" (Wu et al. 2024's phrase). That is
itself a finding: the omission in biorheo's model (dp from fully-developed flow only)
matches the field's own common practice, not a biorheo-specific shortcut — but the
field's practice has not been shown to be *right*, only common.

---

## 3. The hub/taper contraction

| Quantity | Value | Conditions | Source | Confidence |
|---|---|---|---|---|
| Standard Luer taper geometry | 6% convergence angle, 7.5-10.5 mm length | generic Luer/Luer-lock hub | patent literature (not a measurement) | **not citable as a measurement** |
| Effect of replacing the abrupt barrel-to-needle transition with a continuous taper ("tapered needle", TNN) | lower injection force than a standard 27G thin-wall needle for shear-thinning (non-Newtonian) protein solutions; roughly comparable to 27G TW for Newtonian solutions | 27G/29G comparison, multiple model and protein solutions, 3 and 6 mL/min | Krayukhina, Fukuhara & Uchiyama, "Assessment of the Injection Performance of a Tapered Needle for Use in Prefilled Biopharmaceutical Products," *J. Pharm. Sci.* 109(1):515-523 (2020), doi:10.1016/j.xphs.2019.10.033 | **RESOLVED** (DOI/authors/journal/year confirmed via Crossref); the *comparative* finding (TNN < TW for shear-thinning fluids) is corroborated by two independent search summaries, but I was blocked (HTTP 403) from reading the full text directly, so treat as **PARTIAL** |

**Gap, stated plainly**: I found **no paper that isolates what fraction of total
pressure drop the hub/taper region itself carries** (as opposed to comparing two whole
needle designs). Krayukhina et al. compare a fully tapered needle against a
conventional thin-wall needle as complete systems and report the outcome (force
reduction) but do not decompose "how much of a standard needle's pressure drop is the
few-mm conical hub versus the straight bore." This is a real and — as far as I can
tell — unaddressed question in the published injectability literature; the closest
indirect evidence is that a shear-thinning fluid benefits disproportionately from
removing the abrupt contraction, consistent with the hub being a locally
high-shear-rate region, but no number quantifies "X% of dp is the hub."

---

## 4. Peer-reviewed decomposition of total injection force into components

No single peer-reviewed paper was found that decomposes total force into all four
terms named in the prompt (friction + needle + hub + entrance) for a non-Newtonian or
high-concentration protein formulation. The closest results, each covering a subset:

| Component isolated | Value | Conditions | Source | Confidence |
|---|---|---|---|---|
| Friction, isolated by measuring empty-barrel injection separately from filled-syringe injection, then subtracted from total to get the hydrodynamic term | friction force explicitly modeled as its own additive term, "derived from empty barrel injection measurements"; combined model also corrects for actual injection temperature | prefilled syringes, Newtonian and shear-thinning solutions | Wu, Li, Wang, Liu, Zhao, Zhuang, Chen, Zhou & Guo, "Advancing injection force modeling and viscosity-dependent injectability evaluation for prefilled syringes," *Eur. J. Pharm. Biopharm.* (2024), doi:10.1016/j.ejpb.2024.114221 | **RESOLVED** (DOI/authors confirmed, abstract read directly via Semantic Scholar) — but the abstract does not give the friction:hydrodynamic *ratio* as a number, so the split's magnitude is **not quantified** here, only its existence as a separately-measured, non-zero term |
| Subcutaneous tissue back-pressure, isolated in vivo and added as its own term to a needle+friction force model | 0.6 ± 0.5 N to 1.0 ± 0.4 N at ~1 mPa·s; 0.7 ± 0.2 N to 2.4 ± 1.9 N at ~10 mPa·s; 1.8 ± 0.6 N to 4.7 ± 3.3 N at ~20 mPa·s (ranges span the injection rates tested) | in vivo (minipig/porcine model), subcutaneous injection, viscosity-matched Newtonian surrogates | Allmendinger, Mueller, Schwarb, Chipperfield, Huwyler, Mahler & Fischer, "Measuring Tissue Back-Pressure – In Vivo Injection Forces During Subcutaneous Injection," *Pharm. Res.* 32:2229-2240 (2015), doi:10.1007/s11095-014-1611-0 | **PARTIAL** — DOI/authors/journal confirmed; the numeric ranges come from a search-engine summary, not a direct read of the paper's tables |
| Hydrodynamic model for concentrated (non-Newtonian, shear-thinning) protein solutions, with syringe/needle geometry and rheology as inputs | model exists and is validated against measured injection forces; abstract states injection force reduction from shear-thinning is significant relative to a Newtonian-viscosity prediction, but no single friction:hydrodynamic ratio is quoted in what I could access | glass PFS + staked needle, high-concentration mAb | Allmendinger, Fischer, Huwyler, Mahler, Schwarb, Zarraga & Mueller, "Rheological characterization and injection forces of concentrated protein formulations: an alternative predictive model for non-Newtonian solutions," *Eur. J. Pharm. Biopharm.* 87(2):318-328 (2014), doi:10.1016/j.ejpb.2014.01.009 | **RESOLVED** (existence/DOI), but abstract was paywalled/elided in every route tried — I could not read the numbers myself, so no quantitative claim from this paper is reported above beyond what the qualitative summary supports |
| Combination-product force budget (needle + friction + spring/autoinjector force) framework, no single fluid's numeric split extracted | framework exists: barrel size, needle size, autoinjector spring force and "frictional forces" all named as separate, independently characterized terms | high-concentration antibody product, PFS, wide operating range (temperature, concentration, shear) | Rathore, Pranay, Bernacki, Eu, Ji & Walls, "Characterization of protein rheology and delivery forces for combination products," *J. Pharm. Sci.* 101(12):4472-4480 (2012), doi:10.1002/jps.23297 | **RESOLVED** (existence/DOI); numeric split not independently read (abstract elided by publisher in every access route tried) |

**Bottom line for item 4**: the field's best practice (Wu et al. 2024; Rathore et al.
2012) treats friction as a *necessary, separately-measured, non-zero* term — never as
a term you can default to zero — and at least one group (Allmendinger et al. 2015)
treats subcutaneous back-pressure the same way. I could not extract a clean
percentage split (e.g., "friction is 40% of total force for this formulation") from
any paper without full-text access I did not obtain in this session; that specific
number should be flagged as **not found**, not approximated.

---

## 5. Human-factors force limits

| Quantity | Value | Basis | Source | Confidence |
|---|---|---|---|---|
| Average maximum thumb-injection force | 79.8 N (males 95.4 N, females 64.1 N) | **Not an injection-specific measurement.** This is a general finger/pinch force-capability figure, from an unpublished MS thesis, quoted by a review article without re-measurement | Original: Astin, A.D., *Finger Force Capability: Measurement and Prediction Using Anthropometric and Myoelectric Measures*, MS thesis, Virginia Polytechnic Institute and State University, 1999 (no DOI; persistent handle `hdl.handle.net/10919/30923`). Quoted by: Vo, Vardanyan, Cassiano & ..., "The Biomechanics and Optimization of the Needle-Syringe System for Injecting Triamcinolone Acetonide into Keloids," *J. Med. Eng.* 2016:5162394, doi:10.1155/2016/5162394 | **UNVERIFIED as an injection-force limit.** The review's DOI is RESOLVED and I confirmed directly (fetched the paper) that it cites Astin 1999 as the source of this number rather than measuring it itself. This is exactly the "review quoting a review/thesis" pattern to flag: it is a real, locatable primary source (not a phantom citation), but it measures general finger force capability, not thumb force actually delivered while pushing a syringe plunger, and it is a single unpublished thesis, not a peer-reviewed measurement or consensus standard. |
| Isometric pinch/injection force, large (30 mL) syringe, healthy subjects, direct measurement | thenar pinch 103.6 ± 22.9 N; two-handed 104.7 ± 17.1 N (highest of 4 grip types tested) | 20 healthy female hospital workers, isometric MVC against a fixed 30 cc syringe plunger | MacDonald, Wilson, Sonne & Keir, "Grip Type Alters Maximal Pinch Forces in Syringe Use," *Hum. Factors* (2017), doi:10.1177/0018720817720339 | **RESOLVED** (DOI confirmed) — but this is a *maximum voluntary force capability* measurement (how hard someone *can* push), not a comfort/acceptability limit for what they *should* be asked to push during a real injection; and it used a 30 mL syringe, a different geometry from a 1 mL PFS |
| Empirically assessed maximum acceptable injection force for device design | ~80 N found difficult by most subjects; ~40 N assumed as an acceptable design ceiling (taken as roughly half the difficult value) | reported from a study on n≈18 adults using real syringes/injection devices | Verwulgen, Beyers, Van Mulder, Peeters, Truijen, Dams & Vankerckhoven, "Assessment of Forces in Intradermal Injection Devices: Hydrodynamic Versus Human Factors," *Pharm. Res.* 35:113 (2018), doi:10.1007/s11095-018-2397-2 | **PARTIAL** — DOI/authors/journal confirmed to exist and match; every route to the actual abstract/full text was paywalled or blocked in this session (PubMed cookie wall, Springer auth redirect, Semantic Scholar elided-abstract), so the 80 N/40 N figures come from a third-party search-engine summary, not from my own reading of the paper. Do not treat the "40 N is half of 80 N" reasoning as verified methodology — it may be the paper's own stated logic or a paraphrase; I cannot currently tell which. |
| ISO 11608 series ("Needle-based injection systems for medical use — Requirements and test methods") actuation-force requirements | **not resolved** | the standard defines activation-force testing and a "disabled device" force multiplier (≥2x normal max), but I could not access the standard's own numeric force limit | ISO 11608-1 / -5 | **NOT READ.** Per instruction, I am explicitly not quoting any ISO numeric limit — I do not have access to the standard text in this session, only secondary marketing pages (Instron, West Pharma, Mecmesin) that describe the *existence* of force-testing requirements without giving the pass/fail newton value. If a specific ISO force limit is needed, it must be read from a purchased/licensed copy of ISO 11608-1 or -5, not inferred from any of the pages found here. |

**Overall honesty note for §5**: there is no clean, independently-replicated,
peer-reviewed consensus number for "the force a patient/nurse can apply." What exists
is (a) one 1999 unpublished thesis on general finger force, repeatedly re-quoted by
downstream reviews without new measurement (the 79.8 N figure), (b) at least one
direct biomechanics measurement using MVC pinch tasks on a large syringe (component
102-105 N depending on grip, item not injection-specific and not a comfort limit),
and (c) one device-design study that appears to have actually tested injection
difficulty on real subjects and derived an "acceptable" design target by halving an
empirically "difficult" value (Verwulgen et al. 2018) — but I could not read that
paper directly this session, so its number is reported at PARTIAL confidence, not
RESOLVED. **Anyone using a hard number here (40 N, 80 N, or the ISO limit) should
re-derive it from a primary source before it goes into a client deliverable.**

---

## What this means for a model that omits friction and entrance losses

Friction is the dominant omission, not a rounding error: measured breakloose/glide
forces on the same class of device (1 mL glass PFS, staked needle) run roughly 1-10 N
across the studies found here, against the package's shipped hydrodynamic terms of
0.28 N (blood) and 1.87-1.98 N (mAb) — so for the blood case friction alone is
plausibly 5-35x the entire modeled force, and for the mAb case it is comparable to
(and could exceed) the hydrodynamic term; defaulting `F_friction` to zero is not a
second-order simplification for either shipped example, it is dropping the term most
likely to set the answer's order of magnitude. By contrast, entrance/exit (Bagley)
losses, estimated by analogy to capillary-rheometer end corrections rather than any
needle-specific measurement I could find, come out to roughly 1-7% of the needle
pressure drop at the stated L/D range (with the low end, 1-3%, being the more
relevant one for the mildly shear-thinning fluids the package actually targets) —
small relative to friction, and arguably small enough that the field's common
practice of ignoring it (which every injectability paper found here also does) is
defensible, though unproven. The hub/taper contraction and the friction:hydrodynamic
split for a specific formulation are both real gaps in the published literature, not
just gaps in this search: no paper isolates the hub's own pressure-drop share, and no
paper I could fully read gives a clean numeric friction-to-hydrodynamic ratio for a
protein formulation, so the package cannot currently be validated or corrected
against a published split — only against the qualitative fact, well supported across
multiple sources, that friction is real, additive, and of order a few newtons.

# Literature novelty check on Claim 10 (sigma* vs polydispersity at fixed repulsion) - 2026-09-27

Checks Claim 10 of CLAIMS_SIGMA_GEOMETRY_2026-09-24.md: No published measurement of
the shear-thickening onset stress sigma* vs polydispersity at fixed repulsive force
exists.

## Table

| Candidate | Status | DOI / arXiv | Claim it bears on | Note |
|---|---|---|---|---|
| Guy, Ness, Hermes, Sawiak, Sun, Poon (2020) Soft Matter 16, 229 | OK | 10.1039/c9sm00041k / arXiv:1901.02066 | (a) partial pre-empt, (b) pre-empts | Full text read (PDF). Fixed size-independent critical load F star; bidisperse xi sweep; explicit large-first activation order and Delta_11 much greater than Delta_12 much greater than Delta_22 quantification. Simulation plus PMMA experiment. |
| Mari, Seto, Morris, Denn (2014) J Rheol 58, 1693 | OK | 10.1122/1.4890747 | (a) does not pre-empt | Single fixed bidisperse PSD (1 to 1.4, equal volumes); no polydispersity sweep. Confirms the projects own reading. |
| Singh, Mari, Denn, Morris (2018) J Rheol 62, 457 | OK | 10.1122/1.4999237 | (a) does not pre-empt | Constitutive-model paper; DEM at fixed or unspecified PSD, not a polydispersity sweep. Confirms the projects own reading. |
| Pednekar, Chun, Morris (2018) J Rheol 62, 513 | OK | 10.1122/1.5011353 | (a) does not pre-empt | Bidisperse and polydisperse viscosity and phi_J at large solid fraction; no repulsive force and no sigma* onset scan -- a jamming-fraction paper, not an onset-stress paper. Confirms the projects own reading. |
| Guy, Hermes, Poon (2015) PRL 115, 088304 | OK | 10.1103/PhysRevLett.115.088304 | (c) supports | Establishes sigma* proportional to F star over d squared (fixed force, one size at a time) for PMMA -- the scaling claim (c) rests on, not a polydispersity study. |
| Ishibashi, Sugimoto, Seto, Kobayashi (2026) Colloids Surf A 732, 139148 | OK | 10.1016/j.colsurfa.2025.139148 | (c) supports; (a) and (b) do not pre-empt | Experimental, open access (CC BY 4.0). Silica suspensions, separate monodisperse batches (2a = 100 to 2500 nm, not a size distribution within one sample), salt-tuned EDL repulsion; measures critical stress sigma_c vs size and repulsion strength -- a modern, quantitative confirmation and refinement of the sigma* proportional to F over a squared scaling, with a DLVO force law. No polydispersity axis. Co-authored by Ryohei Seto, the LF-DEM electrostatic-repulsion code author this project intends to use -- worth citing for the force-law calibration even though it does not touch CV. |
| Malbranche, Chakraborty, Morris (2023) J Rheol 67, 91 | OK | 10.1122/8.0000495 | (a) cannot tell, abstract only; (c) supports | DEM, lubrication plus Coulomb friction plus EDL repulsion with F_R independent of size (sigma proportional to F_R times a to the minus 2); bidisperse composition and size-ratio sweep. Abstract (AIP page blocked with a 403; read via Crossref plus the CUNY thesis abstract) reports viscosity and phi_J vs composition, a Farris-type non-monotonicity, not sigma* vs CV directly -- full text would be needed to settle whether an onset-stress reading is buried in it. Flag as needs-full-text-check before citing further. |
| Pelosse and Jaeger (2025 or 2026), Softening of shear thickening in suspensions by the addition of large particles, Soft Matter 22, 2379 / arXiv:2501.03529 | OK | PubMed 41757923 (DOI page returned 403, journal record confirmed by title and metadata search); arXiv:2501.03529 | (a) and (b) do not pre-empt as scoped, but this is the closest experimental analogue | Experimental, fumed-silica and nanosilica plus granules, size ratio 20 to 120 (extreme, not moderate CV), no controlled fixed repulsive force. Reports DST to CST softening and onset shift with granule fraction -- corroborates that bidispersity moves the onset, but at a size ratio and control variable well outside the regime this project studies. |
| Singh, Tyagi, Santra (2025 or 2026), Rigidity transition in polydisperse shear-thickening suspensions, Phys Rev Fluids (in press) / arXiv:2510.16464 | OK | arXiv:2510.16464 | does not pre-empt the sigma* claim; adjacent | Full abstract read directly. 2D LF-DEM, high-phi and high-sigma rigidity percolation via the pebble game; phi_c and phi_J^mu non-monotonic in the polydispersity index. This is the jamming-fraction side of the same question, not the onset-stress side -- it does not measure sigma* as a function of CV. Already in the projects own LITERATURE_2026-08-31.md (full text read there on 2026-09-01; the PSDs used are zero-skew normals plus one lognormal). No new information beyond what the project already has on file; listed here for completeness. |

## Verdict, up front

Claim 10 does not stand as worded. Guy, Ness, Hermes, Sawiak, Sun and Poon, Soft
Matter 16, 229 (2020), doi:10.1039/c9sm00041k (arXiv:1901.02066) is already in the
projects own paper/refs.bib under the key guy2020testing, verified there by a citation
agent, but cited only for the WC-model-failure point. It is a DEM simulation (Ness and
Sun lubrication plus Coulomb-friction-above-a-threshold code, same lineage as LF-DEM)
that:

1. holds the critical contact force F* fixed and independent of particle size (the
   fixed repulsive force condition Claim 10 asks about, realised there as a size-blind
   critical load rather than a continuous Derjaguin force law);
2. systematically varies bidispersity (small-particle volume fraction xi at fixed size
   ratio alpha = 0.25) and reports how the onset of frictional engagement shifts with it
   (their Fig 4a and Fig 6a);
3. explicitly measures and states the size-selected activation order that Claim 10 pitches
   as the papers own mechanism: quote, "frictional contacts first form amongst the large
   species, at sigma-star-1 about 0.06 sigma-zero; ... The latest frictional contacts to
   form are the small-small ones ... large-large contacts dominate despite the smallness
   of f11, and many small-large and small-small contacts seem to contribute little to the
   shift in phi_J" (their section 4, page 5). They then quantify it as a per-contact
   weight Delta_ij = kappa_ij divided by f_ij at large stress, with Delta_11 much greater
   than Delta_12 much greater than Delta_22 for every bidisperse composition tested (their
   Fig 8) -- in their words, small particles are largely redundant for stress
   transmission.

This is a bidisperse, critical-load-model measurement, not the continuous lognormal-CV,
Derjaguin-electrostatic, mu = 0.3 system this project runs, so a narrower, correctly
scoped novelty claim survives -- see What remains defensible, below. But the literal
sentence in Claim 10, and the mechanism pitched as the papers contribution (the first
contacts to activate are large-large pairs, so the onset network is size-selected), is
already published, in a paper this project has already read and cited for something
else. This has the same shape as this projects own LITERATURE_2026-08-31.md, which
flagged Guy et al 2020 as something that must be cited by anything here that touches WC
and size distribution together -- a caveat that did not survive into the 2026-09-24
claims file.

## What remains defensible after this check

The literal sentence in Claim 10 is false: a published paper, already cited in the
projects own bib for a different point, measures onset-stress shift with bidispersity
at a size-independent repulsive threshold, and explicitly states and quantifies
large-first size-selected contact activation. Claim 10, and the framing of the mechanism
as novel, need to be narrowed, not deleted. A defensible restatement, given everything
found here:

No published measurement exists of sigma* and the activation sharpness beta as a
function of continuous, lognormal polydispersity (CV = 0.2, 0.35, 0.5, moderate, not the
1-to-4 size-ratio bidisperse extreme of Guy et al 2020, nor the 20-to-120 granule extreme
of Pelosse and Jaeger) under a continuous Derjaguin-scaled electrostatic repulsion and
Coulomb friction mu = 0.3, rather than the size-blind critical-load model at mu = 1 used
by Guy et al. The qualitative mechanism -- large-large contacts activate first and carry
disproportionate weight -- is already published for bidisperse mixtures (Guy et al
2020); the contribution this project could still make is a quantitative test of whether
the same size-selection and its magnitude survive, and how they scale, under a
continuous size distribution and a continuous rather than step-function repulsive force
law, and whether sigma* and beta have a closed-form CV-dependence analogous to the
projects own sigma*-over-kappa contact-geometry identity, Claim 1.

That is a real but substantially smaller claim than Claim 10 as written, and the
size-selected-activation result specifically cannot be presented as a novel finding of
this project -- at most as an independent confirmation, in a different continuous-PSD,
continuous-force regime, of a mechanism Guy et al (2020) already established and this
project has already read.

## Action for the project

1. Rewrite Claim 10 in CLAIMS_SIGMA_GEOMETRY_2026-09-24.md to name Guy et al (2020) as
   prior art for the mechanism, not just for the WC-model-failure point, and to scope the
   remaining novelty to continuous PSD, continuous repulsion, mu = 0.3, per the
   restatement above.
2. If a short paper is drafted, its introduction must cite Guy et al (2020) for the
   large-first activation order and the Delta_11 much greater than Delta_12 much greater
   than Delta_22 quantification, not only for WC fails under bidispersity.
3. Flagged needs-full-text-check: Malbranche, Chakraborty and Morris (2023) -- get past
   the AIP 403 (institutional access or inter-library loan) before citing further; it may
   already contain an onset-stress-vs-composition reading that would narrow the surviving
   claim even further.
4. No claim in this document should be read as approval of anything beyond the
   literature question asked; the DEM and derivation claims 1 through 9 and 11 in
   CLAIMS_SIGMA_GEOMETRY_2026-09-24.md were not in scope here and were not re-checked.

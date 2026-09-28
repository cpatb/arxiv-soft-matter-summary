# How the polydispersity / particle-tuning work becomes papers (framing proposal, 2026-09-27)

Written overnight while `lfsig2` runs. A proposal for the author, not a decision. Every number below is copied from
the file named beside it; nothing here is new data. Status words: **standing** (graded, refereed), **scoped**
(standing with a stated limit), **pending** (tonight), **withdrawn** (must not appear).

## 0. Recommendation in one paragraph

Two papers, with a possible third short one decided by tonight's σ\* result.
**Paper A** is the packing paper already drafted (`paper/main.tex`): the shape of the size distribution, not its width,
sets φ_J. It is nearly ready and should go out first, strengthened by one addition that explains its own rattler
finding. **Paper B** is the paper the programme was built for: *for a drop driven past onset, the size distribution and aggregation act
through two packing ratios, x = φ/φ₀ and χ = φ_m/φ₀* (drop-tested on those two axes; friction and shape are predicted
by the same reduction from DEM χ but were never run through the drop solver). The DNS shows the jam-to-spread switch goes entirely through
χ, and it gives experimentalists a hard warning: the footprint cannot tell a jammed drop from a thickened one.
Paper B needs one bench anchor before submission. The σ\* work becomes either a section of B or a short letter (C).

## 1. Paper A — packing (drafted, 6 pp RevTeX, three referee passes)

**Claim.** At fixed number-mean size and fixed number CV, three PSD families (lognormal, uniform, bidisperse) shift φ_J by
amounts that differ 3–5×. So width is not the control parameter. Two binaries matched in size ratio and number CV jam
0.0635 apart (30× the protocol scatter), and a Furnas-type filling model fixed in advance predicts the pair without a
fitted parameter. Above jamming, polydispersity removes grains from the backbone: the rattler fraction rises ~11×,
while the backbone's coordination and orientation stay family-independent.
(`paper/main.tex` abstract; status: standing, refereed ×3.)

**One addition worth its cost: rattling follows surface area.** Poisson contact counts with mean ∝ a² (no continuous
parameter) give the rattler fraction to rms 0.016 on 702 isotropic records, against 0.128 size-blind. ⚠ The exponent
p = 2 was SELECTED from a sweep of six forms on that same set, so the 702-record number is not out of sample, and since
the law takes the measured Z as input while Z_nr is flat, the total f_r largely follows from Z anyway. The content, and
the clean out-of-sample test, is SIZE-RESOLVED: per species on 135 binaries the sweep never saw, rms 0.017 / 0.018
against 0.133 / 0.346 size-blind (`CLAIMS_SIGMA_GEOMETRY_2026-09-24.md` claim 5; `derivations/rattler_closure` §4d).
That is the *why* behind Paper A's elevenfold rattler rise: polydispersity sorts grains into the backbone by surface
area. Draft section: `paper/draft_rattler_section.tex` (not merged). Cost: one panel (from G2), two paragraphs.
Scope it honestly: multi-class packs spanning ratio 4 do worse (rms 0.11–0.12, under-predicting by ~0.09).

**Leave out of A:** the filling factor w in φ_J = Z/(Z + 2√3 w) as a law (its registered test R5 FAILED, 83 % against a
90 % bar; it is a realistic-friction description only); the K-axis additivity saga; shape/facets.

**Title options.** "The shape of the size distribution, not its width, sets random jamming of frictional spheres";
"Polydisperse spheres jam by filling, not by width".
**Venue.** PRE (or Soft Matter). **Blocking item:** the author line (`paper/main.tex:31`), a standing author call.

## 2. Paper B — particle tuning of drop impact reduces to two packing ratios

**The one-sentence claim.** Once an impact is past the thickening onset, the size distribution and one-primary
aggregation change a dense-suspension drop's fate through two packing ratios, x and χ, with the activation channel
carrying under 0.05 % (CV) and 0.08 % (doublets, U₀ 4). Those are the only two particle axes run through the drop solver.
Friction, elongation, concavity and facets are shown to move χ in dry DEM, which makes them *predicted* by the reduction,
not tested in the drop (`REGIME_MAP_2026-09-26.md`: "Whether friction, shape and potential also collapse onto (x, χ) in
the drop was not tested directly"). The onset stress σ\*, which carries size and surface potential, decides *whether*
the drop reaches that regime.

### Proposed structure and the evidence under each section

| § | content | key evidence | status |
|---|---|---|---|
| 1 | **Reduction.** Wyart–Cates is exactly covariant in (φ, φ₀, φ_m) → (x, χ, σ\*, δ, n); exact phase diagram: x < χ never jams, χ ≤ x < 1 is the shear-jamming window | `PROGRAMME_2026-09-22.md` §1; `derivations/thickening_type_2026-09-25.md` | standing (algebra) |
| 2 | **The drop sees two classes, not three.** CST vs DST is invisible: over 20 drops, peak spread is strictly monotone in Θ = [(1−x)/(1−x/χ)]ⁿ with no jump at x_DST. Drops are either jammable (x ≥ χ) or stiffened by Θ | `PREREG_DSTTYPE_2026-09-26` RESULT (D4 cleanest; D1′ UNDETERMINED) | scoped (one closure, n = 2) |
| 3 | **The knob atlas (DEM).** What each particle property does to χ, with sign: CV narrows the window slightly and lowers x (w 1.00 → 0.76, φ₀ +9.9 %); family shape moves x as much as CV; doublets widen it ×1.49 (χ 0.908 → 0.863, μ_p 0.3, within ±0.005 of ideal mixing); elongation and concavity widen it; friction widens it (0.09 → 0.11); rolling resistance reproduces the doublet's φ_m by the opposite geometry (R3; R7 ρ 0.97 at p̂ 1e-3, NOT INFORMATIVE at p̂ 3e-5), a friction–shape degeneracy worth one sentence. Convex facets RAISE φ with no simulated narrowing | `ANSWER_2026-09-23.md` §2; `REGIME_MAP_2026-09-26.md` §2; `grade_kmer.txt`; `grade_rollchi.txt`; `FACET_RESULT_2026-09-27.md` | CV, family, doublets, friction, rolling: standing. Facets: **descriptive only** (F6 one registered RISE at ε 0.5; F7/F8 INCOMPLETE; node check FAIL) |
| 4 | **The DNS test (the heart).** At fixed weighed φ, CV 0 → 0.5 turns a drop that jams (J 0.29–0.91) into one that cannot, and it spreads 35–67 % further at peak. The window channel carries 99.5–101.8 % of the change and activation < 0.05 %. Doublets jam where singlets cannot (x 0.88), and χ alone reproduces it (J 0.0591 vs 0.0579) | `grade_cvmapE03.txt` (C2a′/C2b′/C3′ PASS); `grade_kmermap.txt` R1; `PREREG_KMERDECOMP` | standing; classes cap-robust (`kmermapE03` 18/18) |
| 5 | **What the bench can and cannot read.** No DNS observable tested separates a jammed drop from a near-edge thickened one: peak spread, terminal footprint and apex recoil all read stiffness. The terminal footprint is sign-mixed (−8.2 … +14.4 %). The largest footprint cut sits where nothing can jam | `PREREG_CVMAP` C4/C2t′; `grade_kmermap.txt` (x 0.85, −31 % Rg) | standing negatives; this is the paper's most useful message to experimentalists |
| 6 | **Onset.** σ\* carries size (∝ a^(m−2)) and potential (∝ κ); the dry-DEM chain gives CV ×0.768 through Zφ. The first lubricated check (`lfsig`, one configuration per arm) was already non-confirmatory: registered readout undefined, post hoc 0.96–1.17, never near 0.768. Tonight's `lfsig2` is the 3-seed redo with a mechanism test. Either way, §4 shows activation is negligible for a drop past onset (G_act ≈ 4–70 for the lab silica), and `kmerdecomp` E1 shows it grows as U₀ falls (1.8 % at U₀ 2, 2.8 % at U₀ 1): that is the scope boundary | `CLAIMS_SIGMA_GEOMETRY` 1–9; `kmerdecomp` Am. 1 | pending |
| 7 | **Predictions.** (i) The collapse test: two suspensions matched in (x, χ, σ_imp/σ\*) but mismatched in size or PSD must impact identically (bench only; trivial in the DNS). ⚠ Not yet well-posed: the impact-stress gauge σ_imp must be named before it runs, because ρkU₀² and the deviatoric wall shear give opposite verdicts on saturation (`PROGRAMME` §3); design the pair so the verdict holds under both. (ii) Fractionated cornstarch has a wider window (lower χ) than Stöber silica at matched CV. (iii) Doublet counts at 0 and 100 mM salt decide whether the salt effect is screening or aggregation | `PROGRAMME` §3; card §3 | predictions |

**What makes it a paper rather than a report.** The thesis is falsifiable (§7i), the reduction is exact, and the DNS
isolates channels by construction, which no experiment can. §5 corrects a reading experimentalists routinely make:
"it spread less, so it jammed."

### Weaknesses a referee will find first, and the answer to each

1. **No bench data.** This is the decisive weakness. The founding observation (CV at similar size → spreading) is the
   author's verbal report and is *not on disk*. **Minimum fix: put that observation on disk as data (drops, PSDs,
   φ), even a small set, and make it figure 1.** Better: the fractionated-cornstarch χ measurement.
2. **The jammed state in the DNS is the viscosity cap** (EPS_MIN; J moves ×5.9 per decade of it). Answer: every claim
   is a CLASS claim, and classes are cap-robust (`epssig`, `kmermapE03` 18/18). Plot classes, never J magnitudes.
3. **One closure, n = 2, L9.** The L10 check passed (R1, F1). Varying n is a registered one-night campaign if needed.
4. **Dry DEM → suspension transfer.** χ comes from frictional dry packs, and the lab φ 0.44 sits below every DEM window.
   Answer: only ratios transfer, and they are stated as such. One transferred quantity, σ\*(CV), has already been checked
   lubricated once (`lfsig`) and did NOT reproduce the dry ratio; `lfsig2` is the redo. The critique's call to calibrate μ_p against the lab rheometer is the real fix, and it needs rheometry the
   project does not yet trust (`si_fit` withdrawn).
5. **"Stress-controlled drop."** Rests on one informative pair (`jamwin` P2). It is stated as a hypothesis in the paper.
6. **Only two particle axes were run through the drop.** CV and doublets. Friction and shape enter only through DEM χ.
   Either scope the thesis to those two (as the claim above now does) or run a friction and an elongation arm through
   the solver at matched χ before submission; that is a one-night, cvmap-sized campaign.

**Leave out of B (withdrawn or superseded):** apex recoil as a jamming signature (refuted out of sample, C4 17/19);
`kmerdrop` (fed the μ_m 3 *predicted* χ); the level-set facet pilot's "lower φ₀"; μ_p ≥ 0.6 clump numbers
(interpenetration); "+81 % spread" (arrested state against a splash peak); the rate-driven jamwin amendment.

**Title options.** "Particle tuning of dense-suspension drop impact reduces to two packing ratios";
"Why polydispersity stops a shear-thickening drop from jamming".
**Venue.** JFM or Phys. Rev. Fluids (a DNS + theory paper with bench predictions); Soft Matter if the knob atlas leads.

## 3. Paper C (conditional) — the onset stress of polydisperse suspensions

Decided by `lfsig2` (graded by ~03:00). The dry-DEM σ\* geometry is refereed: the Love–Weber identity, configuration-model pairing
and rattler law are standing; the closed chain to 2 % is PROVISIONAL (its mean-degree exponent 1.625 comes from two
isotropic cells and is not universal: 1.45 at ratio ≤ 1.5, 2.4–2.9 at 4–5). But on its own it is a dry-packing result about a quantity defined by
lubricated rheology.
- **If P1 = GEO** (lubricated σ\* falls ×0.77): the geometry transfers. That makes a clean short paper, "σ\* is contact
  geometry", and a quantitative input to Paper B §6.
- **If P1 = INERT or NEITHER, with M1–M3 passing:** the more interesting paper. The dry chain describes the *plateau*
  network, while the onset is set by a size-selected sub-network whose geometric factor is CV-neutral. It is a
  mechanism for why dry-DEM σ\* misleads, and it fits a letter.
- **If NOT INFORMATIVE:** fold the geometry into Paper B §6 as a scoped section and do not write C.

⚠ **Prior art (found 2026-09-27, quote verified):** large-first activation in bidisperse suspensions is published —
Guy et al., Soft Matter 16, 229 (2020), doi:10.1039/c9sm00041k, already in `paper/refs.bib`. Paper C must not present the
mechanism as new. What it can own: (i) a continuous PSD swept in CV rather than bidisperse compositions; (ii) a
Derjaguin (electrostatic) threshold rather than a size-independent critical load; and (iii) **the activation exponent**.
Constant force gives σ_act ∝ a⁻² (Guy's ~16× split at size ratio 4), Derjaguin predicts ∝ a⁻¹. If `PREREG_LFSIG2`
Am. 6 S4 returns DERJAGUIN, the letter's headline becomes "the size selectivity of shear-thickening activation is set by
how the repulsion scales with size", which is directly relevant to the lab's charge-stabilised silica. If S4 is MIXED,
Paper C shrinks to a section of B.

## 4. The two-paper picture as one programme

A is the packing law (what φ₀ and φ_m do). B is the drop consequence (the window decides, the footprint cannot tell).
C, if it exists, is the onset. The experimental Paper 1 (800 nm silica, salt → σ\*) sits on B's §6 axis: salt moves the
onset, not the arrested state, which B predicts. Each cites the others, and none depends on an unpublished number from
another.

## 5. Calls only the author can make

1. The author line for A (unchanged since 09-17).
2. Whether to add the surface-area rattler law to A (recommended) or keep A as is and submit this week.
3. What bench anchor B waits for: the founding CV observation on disk (cheapest), fractionated-cornstarch χ, or none
   (submit as theory + DNS with predictions).
4. Whether B leads with the DNS (JFM/PRFluids) or the knob atlas (Soft Matter).

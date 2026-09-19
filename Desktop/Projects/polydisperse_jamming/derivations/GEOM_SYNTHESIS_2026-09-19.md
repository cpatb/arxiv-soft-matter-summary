# The geometry of the polydispersity shift — synthesis, 2026-09-19

> # 🔴 CORRECTED 2026-09-19 by a physics-referee pass — read `geom_protocol_floor_2026-09-19.md` with this page
>
> Every significance on this page was computed from the **seed-paired standard error
> alone**, and this repository documents in its own artifact that that is the wrong bar.
> `mixture_model.py` prints: *"grading scale PHI_SCALE = 2.1e-3: the protocol's own
> run-to-run non-determinism, NOT the seed sem (~2e-04), which is 10x smaller than the
> reproducibility of the thing being modelled"*, and `test_mixture.py` enforces it for a
> structurally identical claim. `\phiScale` is already a macro in the manuscript.
>
> Re-graded against `max(seed s.e., sqrt(L)·PHI_SCALE/sqrt(n))`, with `L` the number of
> cells that do **not** cancel:
>
> | statistic | t on the seed s.e. | t on the protocol floor | verdict |
> |---|---|---|---|
> | item A, additive (lognormal) | 27.70 | **1.67** | 🔴 REFUTED → **UNDETERMINED** |
> | item A, multiplicative (lognormal) | 56.97 | **4.03** | 🔴 REFUTED → **UNDETERMINED** |
> | item C2, zero-skewness pair | resolved 3 of 5 CV | resolved **1 of 5** | conclusion holds, robustness does not |
>
> **What survives:** the additive form fits *better* than the multiplicative one on every
> family — an ordering of two one-parameter forms, which needs no significance test. What
> does not survive is "REFUTED", i.e. the claim that either form is excluded.
>
> ⚠ This is a **stronger** correction than the referee proposed: it computed no separate
> floor for the multiplicative (ratio) statistic and left that one resolved at t = 4.45.
>
> Numbers printed by `geom_protocol_floor_2026-09-19.py`.


Five derivations, each with its own `.md` + `.py` + `.txt` in this directory. Every number
below is printed by one of those scripts; nothing is transcribed. Shared loader:
`geom_common_20260919.py` (excludes `out/phase2_wcbidi*.jsonl` by name — that campaign is
mid-flight).

| item | file | verdict |
|---|---|---|
| A. Separability in (CV, μ_p) | `geom_separability_2026-09-19` | **REFUTED** for multiplicative; **REFUTED** for additive too, but 2.3× less badly. Honest statement: *additive to ~6 % of the shift.* |
| B. The rattler decomposition | `geom_rattler_decomposition_2026-09-19` | **UNDETERMINED** — not impossible (0 of 60 cells refute it seed-paired), not confirmable (the rattler VOLUME fraction is on no shard). |
| C1. Which functional class the binary pair excludes | `geom_psd_functional_class_2026-09-19` | **SUPPORTED** — exact moment algebra; excludes *width*, excludes nothing that reads the third moment. |
| C2. Is it a skewness effect? | `geom_skewness_closure_2026-09-19` | **PARTIALLY SUPPORTED, REFUTED AS COMPLETE** — an external closure predicts the headline gap to 2 se, and is refuted by uniform-vs-bidisperse at t = +6.02. |
| D. Constraint counting | `geom_isostatic_backbone_2026-09-19` | **SUPPORTED** as the reason Z_nr is flat; **REFUTED** as exact — a resolved residual drift, and only at finite friction. |
| E. Shape | `geom_shape_2026-09-19` | **UNDETERMINED** on magnitude, **SUPPORTED** on sign; two of three geometric accounts **REFUTED**. |

---

## 1. The geometric story that survives all five

**φ_J is set by two nearly independent things: a frictional skeleton whose CONNECTIVITY
is size-blind by constraint counting, and a size-distribution effect that acts through
the third moment and is carried by the VOLUME of the grains that fall out of the
skeleton, not by their number.**

The chain, in the order it is established:

1. **The connectivity is size-blind by derivation, not by measurement.** `z_iso =
   d(d+1)/(d − x_m)` contains no length (item D §1). Polydispersity has exactly one
   channel into it — the mobilised contact fraction `x_m` — and that channel exists only
   at finite friction. Measured: at μ_p = 0 the CV drift in z_iso is +0.013 ± 0.030 and
   −0.020 ± 0.034 (both unresolved); at μ_p = 0.3 it is −0.149 ± 0.037 and −0.076 ± 0.011
   (the second resolved at t = −7.25). *The flatness of Z_nr is forced; the residual is
   exactly where the theory says a residual can live.*

2. **The decomposition is an identity, and it names the moment.** `φ = φ_B/(1 − ψ)` with
   `ψ = q·⟨d³⟩_R/⟨d³⟩` (item B §1). The functional that matters is the **third moment
   conditioned on the rattler set**. That is neither d₃₂ nor d₄₃ nor any whole-pack
   moment, which is why four weeks of moment-hunting did not find it.

3. **The number of rattlers is the wrong variable, twice over.** Perturb friction: the
   rattler-count CV response changes by 1.3–2.6× while φ_J's changes by 1.07×
   (item A §6). Perturb shape: the rattler-count CV response changes by 0.996 while
   φ_J's changes by 0.62 (item E3). **Two independent perturbations, and in both the
   count channel and the density channel fail to scale together.** Under the
   decomposition, the carrier must therefore be the conditional rattler *volume*.

4. **The size-distribution effect is, to leading order, a SKEWNESS effect.** Desmond &
   Weeks' published `φ_RCP = 0.634 + 0.0658δ + 0.0857Sδ²` — three constants fitted in
   2013 to 81 frictionless packings that do not include any of this project's cells —
   reproduces the lognormal/uniform contrast at CV = 0.3, 0.4 and 0.5 (all |t| < 1.9),
   and predicts the manuscript's headline matched-binary gap as **+0.0571 against a
   measured +0.0635 ± 0.0032**, i.e. 2.0 se (item C2 §3, §5). For scale, this project's
   own purpose-built two-exponent filling closure predicts +0.0604.

5. **But (δ, S) is not sufficient, and the residual variable is the SUPPORT.** Uniform
   and bidisperse have S = 0 *exactly* and the same δ, so D&W must predict them
   identical; seed-paired they separate at t = +4.90, +7.94, +6.02 at CV = 0.1, 0.4, 0.5
   (item C2 §4). What differs between them is the support: at CV = 0.5 the uniform runs
   on [0.134, 1.866] (size ratio 13.9) while the bidisperse has two atoms at ratio 3.0.
   That is the axis the project's own filling functional reads and D&W does not.

6. **Friction and the distribution add in the DENSITY, not in the void.** Item A's
   transform analysis: additivity in `(1−φ)/φ` or `1/φ` — the free-volume reading —
   requires `R = 0.825` where R is measured at `1.0696 ± 0.0055`, i.e. 44 se away. Of the
   four candidate transforms, the free-volume one is the most strongly refuted. The
   near-additivity that licenses comparing a frictionless published increment against a
   frictional measured one (step 4) is additivity in φ itself.

---

## 2. Derived vs fitted vs assumed

**DERIVED (SymPy-checked, limits asserted, no data):**

* `φ = φ_B/(1−ψ)` and `ψ = q⟨d³⟩_R/⟨d³⟩` — identities.
* `Ψ(q) ≤ ψ ≤ 1 − Ψ(1−q)` — rearrangement, with `Ψ` in closed form for all three
  families (Monte-Carlo checked to ≤ 3e−3).
* `z_iso = d(d+1)/(d − x_m)`, reproducing `d+1 ≤ z_c ≤ 2d` at both ends and in 2D.
* `z_iso = 2n_f/n_c`, and the DOF table that makes Donev's singularity a frictionless
  statement (6 → 10 → 12 frictionless; 6 → 6 frictional).
* The level set `{mean, CV, support}` on two-point PSDs = the two roots of
  `(1+c²)(r−1)²q² + (r−1)[2(1+c²)−r−1]q + c² = 0`, which at r = 3, CV = ½ is
  `5q² − 3q + ¼ = 0` with roots 1/10 and 1/2 — the manuscript's pair, as an output.
* `γ₁ = (1−2q)/√(q(1−q))` for any two-point law, free of the atom values; strictly
  decreasing, hence injective; and the same for every standardized moment.
* `φ_max − β = β(1−β)` in the Furnas ideal limit.
* The transform-separability condition `D(c,m1)/D(c,m0) = h′(φ₀)/h′(φ₁)` ∀c.

**FITTED (someone else's, or this project's):**

* Desmond & Weeks' `(0.634, 0.0658, 0.0857)` — external, published, no uncertainties
  given, and used here only for *increments*.
* The linear packing model's `(PA, PB) = (1.02, 1.50)`, used only in item E4.
* The `√p̂` extrapolation of z_nr — justified by the linear contact law, not fitted.

**ASSUMED (and each is flagged in its own file):**

* **H1**, backbone invariance in CV. Motivated by fact 3, *not implied by it* — a
  structure can keep its contact number and fabric while changing its density.
* Marginal rigidity, which item D needs to read z_iso as a count at all.
* The map `u → u·d_eq/d_s` in item E4 — asserted, not derived, and labelled as the
  reason E4 is "sign only".
* That D&W are frictionless (inferred from their cited algorithm; the word *friction*
  does not appear in their paper).

---

## 3. The single cheapest measurement that would most advance this

> **Record `rattler_vol_frac_k2/k3/k4` — the rattler set's share of the solid volume —
> beside the existing `rattler_frac_k2/k3/k4` in `packing.measure`.**

Why it is the best buy:

* **It costs nothing.** The rattler set is already constructed
  (`packing.iterated_rattlers`), the radii are already in memory, the k = 2/3/4
  convention already exists, and no new DEM cell is needed for any cell that is re-run.
  It is one sum over a list that is already built.
* **It closes item B outright.** `φ·(1−ψ)` becomes directly comparable across CV, and the
  zero-free-parameter test of backbone invariance goes from UNDETERMINED to a verdict.
* **It is the variable both other items point at.** Items A§6 and E3 independently show
  the rattler *count* is the wrong carrier; the decomposition says the right one is the
  conditional volume; this is that quantity.
* **The literature says it matters and has not plotted it against width.** Monti *et al.*
  2022 (doi:10.1103/PhysRevE.106.034901) state the disparity outright — *"rattlers may
  constitute a large fraction of the total NUMBER of particles [but] only a small fraction
  of the total particle VOLUME"* — and Biazzo *et al.* 2010 (doi:10.1063/1.3428387) plot
  both channels for a binary. But **every paper with a distribution-WIDTH x-axis reports
  the number fraction only** (Hoy 2025, Zhang & Scheffold 2015, Cantor 2018,
  Göncü/Durán/Luding 2010). ⚠ The defensible claim is "ψ against CV for a continuous PSD
  at finite friction was not found published", not "nobody has distinguished volume from
  number" — that stronger version would be wrong. Three full texts are unreachable and
  leave the gap **UNRESOLVED** (doi:10.1007/s10035-020-01043-9,
  doi:10.1103/PhysRevE.105.064901, doi:10.1631/jzus.A2500203); ILL them before claiming
  novelty in print.
* **It has a pre-registrable prediction with a sharp number.** If H1 is right, the
  required size selectivity `S = (q − ψ)/(q − Ψ(q))` must come out in **0.75–1.03 at every
  CV ≥ 0.3, across three families and four pressure rungs** (item B §3). An S > 1
  refutes H1 outright; an S that varies strongly with family or rung refutes it as a
  clean picture. For the bidisperse family the prediction is sharper still: at CV = 0.5,
  p̂ = 3e−5, **92.4 % of the small grains and 7.3 % of the large grains must be rattlers.**

**Runner-up, and nearly as cheap:** the mobilised contact fraction
`x_m = #{|f_t| ≥ μ f_n}/N_c`, one more scalar in the same measure step. It converts item
D's unexplained z_iso drift from "either x_m falls or marginal rigidity fails" into a
measurement (predicted: x_m falls 0.355 → 0.265 over CV 0 → 0.5 at μ_p = 0.3, lognormal),
and it is the same physical picture as item B's selectivity — small grains being *held*
rather than slipping.

**Third, and it needs new cells:** a three-atom PSD ladder at fixed number-CV
(item C1 §4). K = 3 gives a continuous skewness axis at matched mean, CV and support,
where the binary pair gives only two isolated points; K = 4 is the first design that can
separate skewness from kurtosis. `psd.make_psd` already accepts the
`"discrete:d1/2/4:y.../..."` spec, so this is a composition grid, not new code. This is
what would turn item C2's "skewness is the leading variable" from a two-family
corroboration into a measured coefficient — and it is the test that would let this project
say something D&W could not, since their |S| > 1.17 range is carried by binary packings
alone.

---

## 4. Things found along the way that bear on existing numbers

1. **⚠ `paper/main.tex:222` says a test has not been done, and it now has been.** The
   manuscript writes that "the only test that could exclude [the skewness reading] is the
   two zero-skewness families against each other (+0.0242 against +0.0139), and that pair
   has not been tested at the calibrated threshold." Seed-paired at the design's own
   2 dof it separates at **+0.01031 ± 0.00171, t = +6.02**, and at two other CVs besides.
   The manuscript's "calibrated threshold" caveat is about a max-vs-min pair *selected*
   from the data; this pair is named in advance by the hypothesis, so a plain t-test is
   the right instrument. **Recommend updating the sentence.** It strengthens the paper:
   it is the observation that refutes a two-parameter (δ, S) closure.

2. **⚠ The μ_p = 1.0 and 3.0 arms fail the project's own ladder-health gate** (median r²
   0.386 and 0.142, one negative ladder slope each). The brief's fact 4 — that the CV
   shift is friction-independent at μ_p = 0 / 0.3 / **1.0**, and that μ_p = 1.0 is
   converged against μ_p = 3.0 — rests at its third level on ladders the project rejects.
   The two-level statement (μ_p = 0 vs 0.3) stands and is what items A and C2 use.
   Nothing in `paper/main.tex` uses μ_p = 1.0.

3. **⚠ `out/z_iso_measured.txt` and `analyse_phase8.py:165-169` give a CV drift 36 %
   larger than a matched-ladder per-seed refit** (−0.232 vs −0.1489 for lognormal; CV = 0
   intercept 4.609 vs 4.537 ± 0.009). The obvious cause — one cell carrying an unmatched
   fifth rung — is refuted by its own sign. Cause not traced; the candidate worth checking
   first is the 2026-09-17 `_it{cap}` duplication, which `sweep_descriptive` fixed and
   `measure_z_iso` may not have. **Not in `paper/main.tex`**; it affects the strain-clock
   line of work.

4. **`z_nr_k4` is non-finite in 216 of 217 μ_p = 0.3 records** — the k = 4 rattler
   threshold empties a backbone whose coordination is only ~4.5. Any analysis that sweeps
   k as a robustness check at this friction is really a two-point sweep. Worth a line in
   whichever document claims a three-threshold check.

5. **`MECHANISM_2026-09-08.md` §M3's "rattlers 0.040 → 0.234"** is the CV 0 → **0.3** pair
   at p̂ = 1e−3, not CV 0 → 0.5 (which is 0.039 → 0.448). `paper/main.tex:162` has the
   right endpoints, so this is a stale line in an older document.

6. **arXiv:1303.4627v2 misprints Desmond & Weeks' `c₂` as 0.857**; the published PRE value
   is **0.0857**. Two independent literature passes in this session returned *different*
   values, so it was settled by arithmetic rather than by lookup
   (`geom_skewness_closure_2026-09-19.md` §0). The decisive check needs no source at all:
   at D&W's own largest posted (δ, S) = (0.413, 2.50), c₂ = 0.857 returns
   **φ_RCP = 1.0266, above the FCC limit π/(3√2) = 0.7405** — the closure would be
   inadmissible inside its own fitted range. Anyone citing the preprint for that
   coefficient gets the skewness term ten times too large.

7. **The monodisperse friction ladder agrees with Silbert (2010) to a few times 1e−3**
   (0.63182 vs 0.639 at μ = 0; 0.56148 vs 0.556 at μ = 1.0, on the fit-free p̂ = 3e−5
   channel). That is the first external anchor on this project's *absolute* φ that I am
   aware of in the derivations directory, and it is a pass.

---

## 5. What I would not claim

* That the separation in (CV, μ) is exactly additive. It is not; the residual is +6.5 %
  of the shift at t = 27.7.
* That no monotone transform separates it. The transform ratio *looks* strongly
  CV-dependent (0.70 → 1.14) but the largest pairwise contrast is t = −3.44 against
  t_crit = 4.303. **UNDETERMINED at n = 3**, and more seeds at low CV would settle it.
* That backbone invariance is true. It is merely not impossible.
* That skewness is *the* descriptor. It is the leading one, and §1 step 5 shows it is not
  the only one.
* That anything here explains the 0.6× shape reduction. Two accounts are dead, one has
  the right sign with an asserted premise.
* Any citation not carrying a DOI/arXiv ID. See §6.

---

## 6. References used, all API-resolved

Every DOI below was resolved by a real lookup in this session; nothing is attributed from
memory. Claim-support notes are narrower than the titles in several cases.

**Load-bearing (a number or a statement in §1–4 rests on it):**

| work | DOI / arXiv | what it is used for | status |
|---|---|---|---|
| Desmond & Weeks, PRE **90**, 022204 (2014) | 10.1103/PhysRevE.90.022204, arXiv:1303.4627 | Eq. (4) and its three coefficients (item C2) | SUPPORTS; ⚠ preprint misprints c₂ 10× |
| van Hecke, JPCM **22**, 033101 (2010) | 10.1088/0953-8984/22/3/033101, arXiv:0911.1384 | `d+1 ≤ z_c ≤ 2d` (§IV.1.1, verbatim); `z_c = 2d` "independent of … polydispersity" (§III.3, verbatim); rattler convention | SUPPORTS; ⚠ μ-monotonicity is protocol-dependent, not a law |
| Donev *et al.*, Science **303**, 990 (2004) | 10.1126/science.1093010 | contact number singular at the sphere point; Z 6 → ~10 then saturating (item E1) | SUPPORTS; ⚠ φ(α) is non-monotone and the authors hedge ("putative") |
| Silbert, Soft Matter **6**, 2918 (2010) | 10.1039/c001973a | monodisperse φ_c(μ) table, the external anchor on this project's absolute φ (item A §5) | SUPPORTS |
| O'Hern, Silbert, Liu & Nagel, PRE **68**, 011306 (2003) | 10.1103/PhysRevE.68.011306 | the measurement behind van Hecke's polydispersity remark | SUPPORTS; ⚠ two distributions, not a width sweep |
| Ozawa, Berthier & Coslovich, SciPost Phys. **3**, 027 (2017) | 10.21468/SciPostPhys.3.4.027 | Z = 6 over the whole J-line for continuously polydisperse frictionless spheres | SUPPORTS |
| Monti *et al.*, PRE **106**, 034901 (2022) | 10.1103/PhysRevE.106.034901 | rattler number ≫ rattler volume, stated outright — external support for item B's high selectivity | SUPPORTS |
| Biazzo, Caltagirone, Parisi & Zamponi, JCP **132**, 176101 (2010) | 10.1063/1.3428387 | rattler number AND volume in one figure — the prior art for item B | SUPPORTS; binary only, axes not a CV |
| Cantor, Azéma, Sornay & Radjaï, PRE **98**, 052910 (2018) | 10.1103/PhysRevE.98.052910 | floater proportion per size class vs span — small grains preferentially rattle | SUPPORTS |
| Furnas, *Ind. Eng. Chem.* **23**, 1052 (1931) | 10.1021/ie50261a017 | the ideal binary limit `β(2−β)` (item E2) | ⚠ DOI resolves, **primary text not read** (ACS 403). The form is algebraically exact and is published explicitly in Farr & Groot 2009 Eq. (3) and Brouwers 2025 Eqs. (16)–(17); note those two cite *different* Furnas publications (1931 vs the 1928 Bureau of Mines RI 2894). |

**Cited as context, not load-bearing:** Farr & Groot, JCP **131**, 244104 (2009),
10.1063/1.3276799; Anzivino *et al.*, JCP **158**, 044901 (2023), 10.1063/5.0137111 (a
cleaner citation than Ogarko & Luding for "φ_RCP rises monotonically with relative s.d.
and saturates"); Ogarko & Luding, JCP **136**, 124508 (2012), 10.1063/1.3694030 (⚠ its
*abstract's* "Z independent of the size distribution" is the compressibility factor, not
the contact number — do not cite the abstract for a coordination claim); Donev *et al.*,
PRE **75**, 051304 (2007), 10.1103/PhysRevE.75.051304; Danisch, Jin & Makse, PRE **81**,
051303 (2010), 10.1103/PhysRevE.81.051303 (the citable instance of z_iso being *assumed*
for arbitrary P(R) and friction — use it for "the assumption is standard", never for "it
is established"); Maher & Torquato, PRE **110**, 064605 (2024),
10.1103/PhysRevE.110.064605; Shaebani, Madadi, Luding & Wolf, PRE **85**, 011301 (2012),
10.1103/PhysRevE.85.011301 (⚠ never uses the word "isostatic"; z is an input, not a
prediction); Petit, Kumar, Luding & Sperl, PRE **106**, 054903 (2022),
10.1103/PhysRevE.106.054903; Hoy, PRE **112**, 045417 (2025), 10.1103/3m8z-qxst; Zhang *et
al.*, PRE **91**, 032302 (2015), 10.1103/PhysRevE.91.032302; Xu, Blawzdziewicz & O'Hern,
PRE **71**, 061306 (2005) — the algorithm behind D&W's frictionlessness.

**[UNVERIFIED] — do not cite as-is:**

* **Vaibhav, Anzivino & Zaccone**, arXiv:2605.30185 (2026), "Theory of distribution
  skewness effect on polydisperse random close packing" — a preprint with no DOI and no
  journal reference that directly extends Desmond & Weeks with a closed form in the
  moments and a skewness master curve. **This is the single most relevant paper to item
  C2 and it has not been read.** Read it before writing up the skewness result.
* Furnas 1928, US Bureau of Mines RI 2894, and Bulletin **307** (1929) — exist, no DOI,
  text not read.
* Oquendo-Patiño & Estrada, 10.1007/s10035-020-01043-9 and 10.1103/PhysRevE.105.064901;
  Wang, Shi & Ling, 10.1631/jzus.A2500203 — DOIs resolve, full texts unreachable, rattler
  number-vs-volume convention **unresolved**. These are the live risks to any item-B
  novelty claim.

⚠ **Integrity note carried forward from the citation pass.** A `WebFetch` summariser
**fabricated a quotation** attributed to Hoy (2025) — *"the volume fraction of rattlers
increases with Δ"* — which would have answered item B's literature question outright. A
string search of that PDF returns zero hits for "volume fraction of rattlers". Every quote
marked *verbatim* above was string-matched in a source, not taken from a summary. Run
`citation-checker` on the final `.bib` before any paper build.

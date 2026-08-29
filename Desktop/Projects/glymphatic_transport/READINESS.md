# Readiness — glymphatic_transport

**Assessed:** 2026-07-22  ·  **Canonical artifact:** code — `glymphatic_peclet.py`
**One-line status:** A correct, well-made scoping figure that reproduces the field's existing consensus — which is exactly why it is not near any publication bar.

## Scores (0–5)

| Stage | Score | Primary blocker to next point |
|---|---|---|
| Preprint / arXiv | 1 | It reproduces the accepted synthesis (Kelley & Thomas ARFM 2023) rather than challenging or extending it. There is no new claim to post. |
| Experimental validation | 2 | No falsifiable *new* prediction exists — a Péclet map built from accepted literature values cannot be confirmed or killed by a measurement. |
| Journal submission | 0 | N/A — no manuscript, no novel result, no target. |

**Recommended next stage:** Experimental validation (highest score, 2)  ·  **Why:** rubric-wise it wins on points, but read the honest verdict below — this is not a publication unit and should not be pushed toward one.

**Honest verdict:** this is a scoping artifact and an entry point to a research area, and it succeeds completely at being that. Its README says so ("This is not a simulation"). Grading it against publication bars is grading it against a job it was not built for. The right decision here is not "how do I publish this" but "do I want to enter this field" — and if yes, the README already names the two builds that would matter.

## What exists (evidence)

- **Re-verified tonight: `python glymphatic_peclet.py` → REGIME MAP VALIDATION PASS.** Crossover lengths at Pe = 1: perivascular small molecule L* = 33.3 µm; perivascular amyloid-β 6.7 µm; parenchyma small molecule 390.6 µm; parenchyma amyloid-β 78.1 µm. Pe (perivascular, 250 µm) = 7.5 → advective; Pe (parenchyma, 20 µm) = 0.051 → diffusive. Matches the field synthesis: **True**.
- Figure `glymphatic_peclet.png` current. `requirements.txt` present.
- Inputs are individually literature-cited (Mestre et al. Nat. Commun. 2018 for PVS speed; Holter et al. PNAS 2017 for parenchymal Darcy velocity; Nicholson tortuosity λ ≈ 1.6), and the README explicitly flags them "verify before formal use."
- Dimensional check: Pe = UL/D is dimensionless as used; the quoted crossover lengths are consistent with the stated U and D (e.g. perivascular small molecule, U ≈ 15 µm/s, D ≈ 5e-10 m²/s → L* = D/U ≈ 33 µm ✓).

## Blockers to the recommended next stage (ranked)

1. **[HIGH — structural]** **The output is a restatement, not a result.** The project's own success criterion is "matches field synthesis." A calculation whose validation is agreement with the consensus cannot, by construction, produce a testable new claim.
2. **[MEDIUM]** **The genuinely contested physics is explicitly excluded.** The README notes that oscillatory Taylor–Aris dispersion can raise effective D by "tens-of-% to several-fold" (Asgari 2016) and that this is "an uncertainty band, not modeled precisely here." That dispersion term is precisely where the advection-vs-diffusion controversy actually lives. The map brackets the argument rather than entering it.
3. **[MEDIUM]** Inputs are literature values carrying their own large uncertainties (perivascular speed 10–20 µm/s, parenchymal 0.1–1 µm/s — a full decade), and no uncertainty propagation is done. Pe values are quoted as point numbers.
4. **[LOW]** No tests, no manuscript.

## Carried-forward flags from prior audits (status)

- No prior formal audit exists (built 2026-07-16). This is the first assessment.
- `OVERNIGHT_BUILDS_2026-07-16.md` claim "matching the field synthesis" — **verified tonight, still true.**
- README's "literature-cited — verify before formal use" caveat on all four input groups — **still open** (none independently verified in this pass).
- README's two named next builds (1-D pulsatile-annulus Taylor–Aris calculator for eccentric PVS cross-sections; a PTV/velocimetry + uncertainty pipeline for perivascular imaging) — **both still open.**

## Concrete next actions

1. **Decide whether to enter this field at all.** This is a strategy call, not a readiness gap. The folder's centre of mass is soft matter and suspension rheology; glymphatics is a different community with different reviewers. Nothing else here depends on this project.
2. **If yes — build the Taylor–Aris dispersion calculator** (README's item 1). That is where a contribution exists: the field genuinely disagrees about whether oscillatory dispersion in eccentric perivascular annuli is large enough to close the gap between measured flow speeds and observed clearance rates. A careful eccentric-annulus dispersion calculation *is* a publishable claim, and this Pe map becomes its framing figure rather than the deliverable.
3. **If yes — the differentiating asset is the velocimetry pipeline** (README's item 2), which plays to existing PTV/uncertainty skill rather than competing on modelling.
4. Add uncertainty bands to the current figure regardless — quoting Pe = 7.5 from inputs spanning a decade overstates precision, and it is a five-line fix.
5. **If no — keep as-is.** It is a good one-page orientation to a field and costs nothing to leave alone. Do not invest further.

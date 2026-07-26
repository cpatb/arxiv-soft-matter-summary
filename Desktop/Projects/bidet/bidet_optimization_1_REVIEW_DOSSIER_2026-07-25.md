# Review Dossier — `bidet_optimization_1.JOURNAL.tex` (CFD-validated)

**Date:** 2026-07-25 · **Pipeline:** validation → postdoc rigor → full-professor significance
**Input:** `bidet_optimization_1.JOURNAL.tex` (CFD-validated upgrade of `.REVIEWED`)
**Method:** three independent review subagents (read-only); consolidated fixes applied here.
**Proof:** compiles (MiKTeX pdflatex ×2), **12 pages, 0 undefined refs, 0 errors**, no
placeholders / `a^{?}` / FLAG comments.

## Verdict
**POST-READY as a preprint (arXiv `physics.flu-dyn`).** Not ready for a primary fluids
journal without the multiphase CFD + GCI study. Recommended journal home:
**European Journal of Physics / American Journal of Physics** (pedagogical/applied fluids
case study). The three verifiable errors the review found have been fixed; the remaining
gaps are scoped honestly in-text and listed below as author decisions.

## Author-decision queue (what only you can resolve)
1. **Venue.** Post to arXiv `physics.flu-dyn` now. For a journal, EJP/AJP is the honest
   home; PoF/PRF would require the VOF/free-surface computation + a 3-grid GCI. Decide
   whether to do that work or submit to EJP as-is.
2. **The larger-domain / multiphase standoff check** (the postdoc's "single most valuable
   addition"): the weak-standoff-lever result is firm for H/d≤20 but extrapolated at the
   toilet-scale H/d≈60 (confined single-phase domain). Running H/d=6,20 on a 2–3× larger
   domain + a momentum-flux budget would upgrade it from "extrapolated" to "resolved."
3. **A direct τ_c measurement.** Under the operative closure Q* ∝ τ_c^{1.43} (strong lever),
   so τ_c is now the dominant absolute-flow uncertainty — the Tier-D microfluidic/hot-film
   measurement is well motivated.

## Fixes applied this pass (verifiable errors caught by review)
- **[FIXED · was a contradiction] τ_c sensitivity used the laminar exponent.** The paper
  argued robustness via Q* ∝ τ_c^{1/3} (2.15× per decade) — the *laminar* branch — while
  every number uses the CFD closure, under which Q* ∝ τ_c^{1.43} (~27× per decade). τ_c is
  a *strong* lever; corrected in the abstract, the τ_c-sensitivity subsection, limitations,
  conclusions, and roadmap.
- **[FIXED · load-bearing number] τ_w(r_max)=5.9 Pa was quoted at the ΔP=0.07 MPa nominal.**
  5.9 Pa is at Re_d=10⁴ (≈40 kPa); at 0.07 MPa it is ≈8.6 Pa (τ_w∝V₀^{3/2}). "Just meets /
  near-optimal" replaced with the honest τ_c-band picture: the commercial nozzle cleans to
  r≈25/21/16/12 mm at τ_c=1/3/5/10 Pa — within ~2× of r_max across the band, over-provisioned
  at the low (best-supported) τ_c, "optimal" is τ_c-contingent.
- **[FIXED · algebra] Oblique a*(θ) exponent** had a spurious factor of 2
  (−3/[2(2n_t−1)] → −3/(2n_t−1)); the Q*(θ) exponent and the 95% number were already correct.
- **[FIXED · scope] Standoff weak-lever headline** scoped to the validated H/d≤20; H/d≈60
  flagged as domain-confined (extrapolated). The ~20% monotonic τ_w rise acknowledged.
- **[FIXED · reproducibility] CFD methods** augmented: radial domain (R=25 mm), BCs, inlet
  profile + turbulence intensity, schemes (linearUpwind/SIMPLEC), residual levels, cell count;
  plus an explicit k-ω SST stagnation-anomaly caveat (closure fit over the developed wall jet).
- **[FIXED · venue framing]** de-inflated from "defensible for Physics of Fluids" to the
  honest EJP/AJP + arXiv recommendation, stating what PoF would additionally require.

## Tier findings (condensed)
**Validation** (partial — the dedicated Tier-1 agent hit a session limit; the postdoc tier
covered the verifiable-error checks): the three algebra/consistency errors above; all fixed.
Braces/environments balanced; compiles clean.

**Postdoc rigor** — 5 findings, ranked by referee risk: (1) standoff-independence rests on
confined/unverified large-H cases → scoped + future-work; (2) τ_c exponent contradiction →
FIXED; (3) near-optimal fragile on ΔP mismatch + τ_c choice → FIXED; (4) CFD verification
gaps (single grid, no GCI, SST stagnation, wall-shear magnitude unvalidated) → methods detail
+ caveats added, GCI/multiphase scoped as future work; (5) oblique a* exponent → FIXED.
Well-handled already: the τ_c *framing*, the coherence caveat, the Region-I/II matching honesty.

**Professor significance** — "so what" is modest: application + CFD-corrected design reading,
not new jet physics; the standoff "correction" fixes the paper's own earlier draft, not a
literature claim. Novelty vs the Wilson/Bhagat cleaning-jet literature is thin but honestly
argued. Venue: PoF unrealistic; EJP/AJP + arXiv is the home. Top objections: τ_c fragility
(addressed), closure validated against free-jet decay not wall-shear magnitude (caveated),
single-phase excludes the regime where standoff could fail (caveated). **Go to preprint,
refine before journal.**

## Bottom line
The scaling framework is sound, the CFD removes the calibration blocker, and the three
verifiable errors are fixed. The paper is honestly scoped, internally consistent, and compiles.
Post to arXiv `physics.flu-dyn`; target EJP/AJP for the journal version; the multiphase CFD and
a direct τ_c measurement are the two highest-value follow-ups.

# Three questions, 2026-08-03

Answers to: (1) can the TPM recipe be doubled, (2) research questions for Newtonian drop
impact onto Carbopol, (3) next steps on the complex-substrate simulations.

---

## 1. Doubling the TPM synthesis

Short answer: doubling is safe, but don't do it tomorrow, and it won't fix the yield.

### Why the scaling is safe

Every knob that sets particle size in this route is intensive:

| knob | what it is | preserved by an isotropic 2x? |
|---|---|---|
| pH / [NH3] | primary size control, sets nucleus number | yes, a concentration |
| phi_TPM (monomer vol. fraction) | secondary size control | yes, a concentration |
| temperature, time | ripening plus polymerization | yes |
| the flow field | sets how fast the TPM bolus disperses | **no, and it is the only one** |

van der Wel (2017) runs this chemistry at 15 mL in a 100 mL beaker. The collaborator
protocol runs it at 100 mL in a 1 L beaker and gets the same size range, which is already a
6.7x scale-up on the same recipe logic. That is the strongest evidence you have that
absolute volume is not a size variable. The docx anticipates scaling explicitly: "adjust to
larger beaker, keep ratio of stir bar to beaker diameter the same."

The caveat comes from Crothers (2023), the same paper the SOP already cites. Their finding
is that reproducible size control turns on *how the reagents are mixed*, which is why they
went to a zero-flow regime to get it. Mixing is the variable that scale changes, so it is
the one to pin down.

### What to hold fixed when you do double

Same vessel diameter, same stir bar, same 400 rpm. Going 100 to 200 mL in the same 1 L
beaker takes the liquid from about 1.3 cm to about 2.5 cm deep at fixed diameter, leaving
the bar's swept area and the bar-to-diameter ratio unchanged. That is a smaller perturbation
than moving to a 2 L vessel.

Inject the TPM bolus below the surface, into the stir bar's shear zone, at the same point,
in one go. 4 mL of neat oil floats, and this is the step most likely to disperse more slowly
at 2x. It is exactly the mixing variable Crothers is warning about.

Check that the emulsion is uniformly milky top to bottom at about 15 min. A clear upper
layer means the blend time went up, and the fix is more rpm rather than accepting it.

### What genuinely changes

**Heat-up.** At fixed diameter, doubling the depth takes V/A from about 0.84 to about
1.26 cm, so the bath time constant rises roughly 1.5x. Start the 3 h clock when the internal
probe reads 75 C, not when the vessel goes in. AIBN's half-life at 75 C is about 30 min, so
3 h is six half-lives and an extra 5-10 min of ramp is inside the noise, provided the clock
starts at temperature.

**Washing.** 200 mL means 4-5 Falcon tubes per wash, three washes. Budget 1 to 1.5 h rather
than 30-45 min, and check the bucket count.

**Cleaning.** The solvent soak has to cover the higher fill line.

**Not the oven, at 2x.** Manual periodic agitation of a deeper, denser column is where a
doubled batch will actually go wrong. Use the stirred bath or an oil bath.

### The yield arithmetic, which is the part that matters

Full condensation of TPM (MW 248.35) releases 3 MeOH per monomer, so the solid keeps about
176/248, or 71%, of monomer mass. Per mL of TPM monomer:

    1.045 g x 0.71 = 0.74 g solid  ->  about 0.58 mL of particles  (rho = 1.23-1.31 g/cm3)

At the Master Recipe's 2.0 mL TPM per 100 mL:

| | particles | suspension at phi = 0.44 |
|---|---|---|
| 100 mL batch | ~1.2 mL | **~2.6 mL** |
| 200 mL batch | ~2.4 mL | **~5.2 mL** |

A full 6-7 h day currently buys about 2.6 mL of usable suspension. Your instinct that the
yield is low is right, and it is worse than "a little."

A phi sweep of four volume fractions at ~3 mL each needs about 4.5 mL of particles, which is
about 7.8 mL of TPM, which is about 390 mL of reaction volume. That is roughly **two doubled
batches per sweep**. So 200 mL is not a luxury, it is close to the natural production unit.
But it only pays if both batches give the same particle size, and that needs a calibration
you don't have yet.

### So: run tomorrow's batch at 100 mL

1. Batch 1's job is to produce two numbers that don't exist yet: the [NH3] to size point for
   *your* setup, and the actual yield from the Step 14 mass-loss assay. The SOP says outright
   that batch 1 is what establishes the expected value. Neither needs 2x material, and both
   are needed before a production batch means anything.
2. It is the batch most likely to fail, and the batch where the F108 question gets decided.
   Doubling doubles the cost of a failure that is genuinely likely.
3. Calibrating at 200 mL locks you into always running 200 mL. That is fine as a choice and
   bad as a side effect of wanting more material on day one.

Note the ordering this implies. The other yield lever, raising phi_TPM from 2% to 4-5% (which
van der Wel supports), roughly halves the batch count but *increases the particle diameter*,
by about 1.36x going 2 to 5% at fixed nucleus number, since d scales as phi^(1/3). Diameter
sets gamma_dot_c, which is the entire measurement. You would have to compensate with pH, and
that is a 2-D calibration you don't have. **Volume is the only scale knob that doesn't move
the physics.** Calibrate at 100 mL, then decide about phi_TPM, then scale volume and only
volume.

### Three things that will bite regardless of scale

If you run the original docx cover to cover:

1. **No TPM volume is specified for plain particles.** "If making plain particles (regular
   TPM): Add TPM to solution." The 3 mL figure earlier in the doc is the *hydrolysis* prep
   for microrollers (3 mL into 20 mL MilliQ), not the synthesis charge. Without a TPM volume
   you have neither a size target nor a yield. Take 2.0 mL per 100 mL from the revised SOP;
   it sits mid-range of van der Wel's productive 0.5-3 mL.
2. **AIBN at 1 mg/100 mL is 20 to 67x below both published recipes** (van der Wel
   67 mg/100 mL, Crothers 20 mg/100 mL), about 0.05 wt% of TPM. The failure it produces, a
   clear or tacky film instead of powder, reads identically to dead AIBN, so you would end up
   replacing a good bottle. Use 1 wt% of TPM mass, so 21 mg at 2.0 mL TPM (42 mg doubled).
   That is the conservative end of the published range, not an extrapolation.
3. **60 uL ammonia per 100 mL is the large-particle end** (pH about 10.6). For 0.6-0.9 um you
   want 0.40 mL. Either is defensible, but read the pH back against the table for whatever
   volume you actually used, and log the actual volume rather than the nominal one.

F108 and scaling don't interact. It goes in after the 1.5 h pre-stir, so after nucleation and
ripening have already set the size, which makes it the least scale-sensitive step in the
protocol. 1 mL of 10 wt% stock per 100 mL becomes 2 mL, and that's all. The real F108 issue
is the one already on the card: the brush sits in the frictional contact path and displaces
the shear plane, so zeta is confounded by the layer rather than reporting it. Wash and
measure zeta against wash number.

---

## 2. Newtonian drop impact onto Carbopol: the research questions

The strongest asset here is unusual. The simulations already made quantitative, falsifiable
predictions before any data exists. Frame the questions as tests of those and the paper has a
spine that most experimental papers have to invent afterwards.

### Q1. Which observable reports on the substrate? (the primary question)

The simulation predicts crater depth moves about 80x across the tau_y range while lambda (the
splat radius) moves 33% of its own mean, non-monotonically, and the water control's lambda
still rises 91% over U0 = 1 to 3.

> Across a 10x range in tau_y at fixed impact speed, does the maximum splat radius change by
> less than it changes across the same speed range on a Newtonian bath?

Clean, sign-carrying, numerically pre-registered. And the answer is a methods result with
teeth: the observable most impact papers report is the one carrying the least information
about the substrate.

### Q2. Does the substrate response collapse on J = tau_y/(rho U0^2)?

Your 4 x 4 grid spans J over roughly two decades. The test is that crater depth over D0
against J should collapse 16 cells onto one curve. If it doesn't, the deviations are the
result.

⚠ The single most important design instruction in this document, and the simulation program
already made and caught this mistake: **never compute the collapse on survivors.** If one
concentration fails at one speed, drop the whole row or column from the statistic rather than
averaging over what is left. On the uneven grid the same statistic read 1.46-2.63 against
0.67-0.89 on a balanced one, a sign reversal that was pure grid geometry. Either land all 16
cells or compute the collapse only on a balanced sub-grid, and say which one you did.

### Q3. Is J alone enough?

The simulations found a real secondary correction (J_sigma, about 9%), and the confound
separating it from J is still open. The card notes no clean bench experiment exists for it
*in the solver's framing*. But on the bench you can move the bath's surface tension nearly
independently of tau_y, using a trace of surfactant at fixed polymer. That is a case where
the experiment reaches something the simulation can't. Stretch goal, not the main grid.

### Q4. Does the water control reproduce known Newtonian behaviour?

The four speeds into plain water, using the identical seeded fluid. Prediction: lambda rises
91% over U0 = 1 to 3. This is what turns "lambda doesn't respond to tau_y" from a statement
about noise into a statement about physics, because you show your lambda measurement does
respond to something, just not to the substrate.

### Q5. Is there a side-view proxy for the fluidized fraction?

The card is right that the fluidized fraction is a volume integral with no side-view
equivalent, so don't chase it. But asking whether *any* side-view observable correlates with
the simulated fluidized fraction is worth one figure. A hit is a calibrated bridge between the
two halves of the program, and a miss is a clean statement that the simulation carries
information the experiment cannot.

### The shear-thickening extension, sharpened

Your framing (bridging the power laws between Newtonian fluids and solids) has a specific
shape in the literature. Solid-sphere impact into yield-stress fluids gives penetration
scaling laws. Newtonian drop impact gives different ones, because the drop deforms instead of
penetrating. The empty cell is a projectile that is fluid but becomes rigid during the impact.

> Plot crater depth over D0 against J for (a) water drops, (b) shear-thickening drops, and
> (c) solid spheres from the de Bruyn literature. Prediction: the STF drop lies between (a)
> and (c), and its position is set by the threshold ratio Lambda = sigma\*_drop / tau_y,bath,
> not by U0 independently. At high Lambda it converges on the solid-sphere law, at low Lambda
> on the Newtonian-drop law.

Both materials see the same contact pressure at the same instant, so the physical question is
which threshold is crossed first. Lambda is U0-independent, which makes it a real second axis
rather than a restatement of J.

This is a genuine paper, and it should not start before the 20-run grid closes, because that
grid *is* leg (a). You cannot place the STF drop between two curves you have not measured.

Use cornstarch, not TPM, for the first pass. It is available now, its jamming is strong so the
contrast against water is large, and it is already opaque, which means that arm doesn't
inherit the silica-seeding confound at all. Its weaknesses (settling, batch drift, no clean
phi axis) don't matter for a "does this land between the two curves" question. And spending
early TPM batches here competes directly with the calibrated-phi-axis work that is the whole
point of the synthesis pivot.

### Four operational points that decide whether the data is usable

1. **tau_y per batch is non-negotiable**, and it is worse than the card states. Carbopol's
   yield stress depends strongly on neutralization (NaOH), not just wt%. Two batches at the
   same wt% with different neutralization differ substantially. Flow curve per batch, or
   retain the batch.
2. **Bath depth at least 5 drop diameters**, so at least 15 mm for a 3 mm drop, and 25 mm to
   be safe. This transfers directly from the HBATH = 5 D0 finding: at 2.5 D0 the crater felt
   the rigid floor and floor-limited three of four Newtonian controls. Cheap insurance.
3. **Fix and record the segmentation threshold before acquisition.** This is the exact error
   class that cost Paper 1 five of six contrasts, a per-condition processing parameter that
   varied and was never written down. Decide it now, hold it across all 20 runs, put it in
   the config.
4. **Agitate the seeded reservoir between drops**, and treat opacity drift across a session as
   a concentration-drift suspect first.

---

## 3. Complex-substrate simulations: next steps

One correction first. The grant card's resume says "Everything from 03 AUG is UNCOMMITTED."
That is stale, and `14cadfc` has it. The only untracked file in the repo is
`basilisk-drop-impact/DOMAIN_SIZE_2026-08-03.md`. Worth fixing on the card, since an
inaccurate "your work is at risk" line trains you to ignore the ones that are true.

### The ordering I would use, and where I depart from the card

**1. Fix `paper_stf`.** Unchanged, still first. A cleanly-building manuscript that typesets a
withdrawn claim is the most dangerous artifact in the repo, because it looks finished. The
sentence-by-sentence delta already exists.

**2. Port the crater fix into `hb_drop_on_bath.c`, and move this up.** The card has it at #2
behind the m ladder's follow-up, but the Carbopol acquisition starts tomorrow, and
`paper_yieldstress` runs on the HB solver, which carries the same dropout pathology with no
flag at all. The paper that is about to acquire matching experimental data is the one running
the unfixed detector. Comparing a fresh dataset against a known-broken diagnostic is how you
manufacture a discrepancy and then spend a week chasing it.

**3. Settle the regularization question, and run the cheap diagnostic first.** This is the
real blocker. None of the five candidate papers can quote an absolute depth until L10 is
reachable at a defensible setting, and per-cell m is not one. "A per-cell tuned parameter
cannot be defended in print" is the card's own verdict and it is right.

The benchmark is now well posed: 10 of 12 at L10 at one fixed setting. Before committing to
the expensive route, run Bercovier-Engelman, mu = mu_p + tau_y/(|gamma_dot| + eps), roughly a
one-line change from the clamp, purely as a diagnostic. If BE and Papanastasiou fail on the
*same* cells, the problem is the grid or the physics rather than the regularization, and no
amount of regularization work will help. That costs about an hour and it tells you whether
the expensive route is worth taking.

If it is, augmented Lagrangian (Uzawa) is the principled choice, because it converges to the
true Bingham solution instead of approximating it. That leaves no regularization parameter to
defend, which is exactly the objection that killed per-cell m. More implementation work,
slower per step.

**4. Texture BC, run concurrently on other cores.** The brainstorm ranks this highest by
leverage per line and the reasoning holds: spec written, axisymmetric, no rheometry, and it
unblocks both the STF-on-texture paper and the third leg of the Theta-collapse test. The
reason it can run now rather than queuing behind item 3 is that it is a wall BC on the
rigid-solid side, so it shares neither the rheometry blocker nor the bath solver's `u.x`
viscous problem. It is the one track blocked by nothing that is currently blocking.

**5. Write the observable-selection note (§3.1) this week, specifically because of the bench
schedule.** Every number already exists, and the Carbopol experiment is about to test its
central claim. Writing it before the data lands pre-registers the prediction instead of
rationalising the result afterwards.

### Where I would push back

The Theta-collapse test (§1 of the brainstorm) is the most attractive item in that document
and the weakest near-term bet. It needs three solvers all trustworthy at one common
resolution, which is precisely what is not established. It is the right destination, not the
right next step, and it gets cheaper the moment item 3 lands.

Separately, the card is 870+ lines and says so itself. Push §2's history down into the session
docs it already cites and keep §2 to the state.

---

## Sources

- van der Wel et al., *Preparation of Colloidal Organosilica Spheres through Spontaneous
  Emulsification*, Langmuir 33, 8174 (2017).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5569361/
- Crothers, Orr, van der Meer, Dullens, Yanagishima, *Characterization and Optimization of
  Fluorescent Organosilica Colloids ... Prepared Under "Zero-Flow" Conditions*, Langmuir 39,
  5306 (2023). https://arxiv.org/abs/2304.07049
- `Desktop\Driscoll\Protocols\TPMSynthesisProtocol.docx` (original collaborator route)
- `Desktop\LaTeX Docs\tpm_synthesis_SOP.tex` and `..._BENCH_READINESS_03AUG26.md` (revised)
- `Status\Drop Impact - Paper 1.md` §2b, §2c; `Status\Drop Impact - Substrate Grant.md` §3
- `basilisk-drop-impact\NEXT_PAPERS_COMPLEX_SUBSTRATES_2026-08-03.md`

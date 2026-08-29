# One permeability, two clocks: bridging silica and cornstarch drop impact through pore size

**2026-08-10 overnight session.** Companion to `PORE_PHYSICS.pdf` (pore-relaxation
analysis, 2026-08-10) and `PORE_UNJAMMING_ADJUDICATION_2026-08-08.md`. New inputs
tonight: the published ID paper (Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar,
Driscoll & Cheng, PRL 136, 148201 (2026), plus its SI), Brown & Jaeger (J. Rheol.
56, 875 (2012)), Royer, Blair & Hudson (PRL 116, 188301 (2016)), and the 18-cell
`trsweep_2026-08-10` DNS campaign.

## 1. The claim in one paragraph

The Kozeny-Carman permeability kappa(a, phi) = (1-phi)^3 a^2 / (180 phi^2) controls
both ends of the jamming cycle, with opposite powers of particle size. Jamming under
impact is set by the dilatant timescale t_d ~ rho*kappa/(A mu_w dphi), which goes as
a^2 (Mobaseri et al., Eq. S23): tighter pores resist Darcy backflow harder, the
dilatant particle pressure is larger, and momentum is absorbed faster. Un-jamming by
drainage is set by the poroelastic time tau_p ~ mu_w L^2/(kappa K_p), which goes as
a^-2: the same tight pores that jam a silica drop instantly also strangle the
backflow needed to relax it. The ratio of un-jam to jam times therefore scales as
tau_p/t_d ~ a^-4, and that single scaling puts the two materials at opposite corners
of the phase space:

| | a | ID = t_d/t_i | jams | tau_p (contact scale) | un-jams |
|---|---|---|---|---|---|
| silica colloid | ~0.5-1 um | ~1e-5 | instantly (always dilatant regime) | 790-1430 s | pore-limited creep, seconds and up |
| cornstarch | ~14 um | ~0.2-20 | slowly enough to go inertial | ~ms (De = 7) | fast, ms scale |

This is the quantitative form of the observation that silica and cornstarch sit at
two ends of the shear-thickening drop-impact spectrum. The ends are the two limits
of one permeability, not two different kinds of physics.

## 2. What the DNS said tonight (trsweep_2026-08-10, 18 cells, L9)

The un-jamming clock t_r of the measured relaxation law (WC_CLOCK=1,
dfhat/dt = -(gd/gamma0 + 1/t_r)(fhat - fss)) was swept over four decades, 1 ms to
10 s, across both calibrated presets. Gates are in the `analyze_trsweep.py` output.

1. The impact phase does not feel the clock at all. Peak jam_frac is identical to
   six digits across t_r = 0.1-10 s (si arm: 0.997951 in every cell), and the
   beta_max spread is exactly zero. Pre-registered P1 passes; master-plan item 4
   ("sweep t_r over the measured bracket") is closed. During impact the strain term
   gd/gamma0 is a channel roughly a thousand times faster, and everything else is
   invisible next to it.
2. After impact, relaxation is still strain-dominated unless t_r is faster than
   about 30 ms. The post-peak decay rate of the mean fhat saturates at a floor
   k = 0.75/tau_ic for every t_r of 30 ms or slower; only the 3 ms and 10 ms cells
   decay faster (2.3x and 1.9x the floor). The crossover is bracketed in
   t_r* = 10-30 ms, right at the cornstarch pore-drainage scale (their measured
   front time t_d = 3.8 ms) and three to five decades away from silica's tau_p.
3. The silica arm's decay is not a single exponential (all nine fits fail the rms
   gate at 0.50-0.65, partly from a recoil re-jamming bump at t/tau near 1). That
   is an honest UNRESOLVED, and a hint that even the flowing relaxation of the
   stiff arm has structure worth its own campaign.
4. Cornstarch retains more state than silica at every matched clock (8 of 9
   pairs): the sigma*/phi_M asymmetry again. The lower onset stress wins,
   independent of t_r.

The pre-registered falsifier fired, and the conclusion is sharp: in flowing or
deforming material, un-jamming is strain-driven, and a pore clock can only express
itself in arrested material. That reconciles two things that looked contradictory:

- The 2026-08-08 adjudication found that flow-cessation release times do not scale
  with particle size, which ruled out an a^-2 clock. But those experiments watch
  material that is still relaxing by flowing, where the strain channel dominates
  for any t_r above 30 ms. An a^-2 dependence in t_r would have been invisible to
  them.
- The bench phenomenology: silica solidification fronts hold for seconds, and the
  common-support front-speed work reads their slow collapse directly. Arrested
  fronts have gd near zero, which is exactly where 1/t_r is the only channel left.

Consequence: setting t_r = tau_p(a, phi), the pore-limited-creep hypothesis, is
safe and testable. It changes nothing for any material with tau_p above 30 ms in
any flowing configuration (which is why the literature cannot see it), it
reproduces cornstarch's fast relaxation, and it makes a hard prediction only in
arrested-state holds: release time proportional to a^-2 at fixed chemistry and
phi/phi_M. The adjudication's proposed experiment (particle-size ladder,
stress-drop release time) remains the decisive test, but it must be done on an
arrested configuration, not a flow-cessation one. A null in flow cessation no
longer rules the hypothesis out.

## 3. Where the bridge lives in the solver

A general STF-impact solver spanning silica to cornstarch needs these parameters,
now each with a physical anchor:

| knob | physics | a-scaling | status |
|---|---|---|---|
| sigma*, delta, phi_M, phi_0, eta_s | WC flow curve | sigma* ~ F_rep/a^2 (Royer 2016) but chemistry-dominated | calibrated (si_fit/cs_fit); Royer's silica set (phi_m = 0.592, phi_0 = 0.711, f = f_max exp(-sigma*/sigma)) is an external cross-check |
| gamma0 = 0.2 | strain scale of the state clock | none (universal; Richards, Han) | measured |
| t_r | rest clock | hypothesis: t_r = tau_p, goes as a^-2 | swept tonight; invisible above 30 ms in flow; needs the arrested-state test |
| DIL_COEFF | dilatant particle pressure | ~ mu_w/kappa, goes as a^-2 | implemented 2026-08-09; DIL_crit measured, Kozeny-Carman phi-ratio confirmed to 1.5-2.8% |
| DIL_CEIL | saturation of the dilatant stress | bulk: K_p dphi/phi (a-free); boundary: 0.1 gamma/a (Brown & Jaeger 2012) | implemented tonight, 5 gates pass; the smaller scale binds |
| FORCE_MODE | substrate force F(t) | - | implemented tonight; the L8 smoke test already reads Fmax = 0.851 at t = 0.155 against the universal 0.85 +/- 0.05 at 0.17 +/- 0.02 |

Two identifications worth writing down:

- DIL_COEFF is essentially 1/ID. DIL_COEFF = (mu_w/kappa) eps_dil L_d^2/(rho U0 D0)
  and 1/ID = A mu_w dphi D0/(rho kappa U0) differ only by the order-one geometric
  factor (L_d/D0)^2 eps_dil/(A dphi). Sweeping the amplitude, as the dilatancy
  campaigns do, is sweeping the inertial-dilatant number, and the measured jamming
  threshold DIL_crit(phi = 0.44) = 0.178 is this programme's DNS estimate of their
  ID = 1 regime boundary, up to that geometric factor.
- The ceiling is also a-dependent, in the opposite direction to the drive. The
  dilatant demand grows as a^-2 while the free-surface confining ceiling
  (0.1 gamma/a) shrinks only as a^-1, so smaller particles both demand more
  dilatant stress and can sustain more of it (silica: about 7 kPa interface
  ceiling vs cornstarch's 0.5 kPa), unless the bulk packing stiffness K_p dphi/phi
  (a-free, about 27 Pa at face value for our silica) is the compliant element.
  Which ceiling binds is exactly the K_p measurement question on the Paper 1 card.

## 3b. Late addition (same night): the ceiling experiment ran, and it repairs yesterday's "50x" claim

`dilceilsweep_2026-08-10` (7 cells, phi = 0.44, DIL_COEFF = 1.0 fixed, ceiling
swept in half-decades, jam_frac read at t/tau = 2):

| DIL_CEIL (rho U0^2) | 0 (none) | 0.0035 | 0.011 | 0.035 | 0.11 | 0.354 | 1.12 |
|---|---|---|---|---|---|---|---|
| jam_frac | 1.000 | 0.008 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

- C1 and C2 pass as pre-registered: no ceiling reproduces saturated jamming, and
  the face-value bulk ceiling (0.0035 = 26.8 Pa) abolishes it.
- C3's quantitative expectation failed in the informative direction. CEIL_crit is
  bracketed in [0.0035, 0.011] (geometric midpoint 0.0062 rho U0^2, about 47 Pa),
  a factor of roughly 29 below DIL_crit = 0.178, and sigma*_code = 0.0092 sits
  inside the bracket. The dilatant pressure that must actually be delivered to jam
  the drop is just the onset stress sigma*, which is obvious in hindsight (fss
  needs sigma above sigma*) and repairs a real error: PORE_PHYSICS compared the
  swept coefficient (0.178, which multiplies a strain rate that is often much less
  than 1) against the packing's pressure ceiling, and concluded a 50x
  impossibility. Comparing pressures with pressures, the requirement is
  K_p at least CEIL_crit * rho U0^2 * phi/dphi, which is 1.3-1.8x the face-value
  K_p estimate. Dilatant jamming under a physical ceiling is marginal, not
  excluded, and the K_p bench measurement decides it. (Ladder resolution:
  half-decades, so the bracket is a factor 3.2. Quote the bracket, not the
  midpoint.)

## 4. What to do next, in order

1. Arrested-state hold in the DNS (cheap, decisive for the model): take a jammed
   deposit, let the flow die, and measure the fhat release time as a function of
   t_r alone. This is the configuration where the clock is the only channel, the
   DNS analogue of the silica front holds. If the release time equals t_r there,
   the solver is ready to carry t_r = tau_p(a, phi). [Ran the free-running version
   the same night: H2 negative, deposits do not arrest cleanly within 20 tau, so
   frozen-dump restarts are the required instrument.]
2. F(t) validation campaign: Newtonian control at L9/L10 against the universal
   curve (quantitative gate: 0.85 +/- 0.05 at 0.17 +/- 0.02, early (3 sqrt6/2)
   sqrt(t) with the Re^-1/2 correction), then cs_fit cells at Mobaseri's (phi, U0)
   grid against their Figs. S10-S15. This is the first external, published,
   quantitative validation target the STF arm has.
3. K_p measurement (bench; two cheap routes already on the Paper 1 card). It
   decides the binding DIL_CEIL and converts the dilatant term from swept to
   predicted.
4. Particle-size ladder on arrested holds (bench): chemistry and phi/phi_M fixed,
   a varied, release time of an arrested front or stress-drop measured. An a^-2
   result confirms pore-limited creep; a null is publishable against it.

## 5. Caveats

- Everything DNS-side above is L9, one impact speed (U0 = 2), one phi per preset.
  The t_r crossover bracket of 10-30 ms will move with U0 and TEND because the
  strain floor gd_res is set by how fast the deposit is still creeping.
- The si-arm decay shape is unresolved. Do not quote a silica decay time from this
  campaign.
- The tau_p numbers carry the unmeasured K_p through them (the De verdict is
  robust over three decades of assumed modulus, but absolute times are
  order-of-magnitude).
- sigma* comparisons across materials fold in chemistry (F_rep), so the a^-2 in
  Royer's sigma* scaling must not be used to predict si vs cs onset ratios.

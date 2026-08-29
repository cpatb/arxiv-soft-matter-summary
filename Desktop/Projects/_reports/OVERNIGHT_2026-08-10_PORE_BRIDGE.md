# Overnight 2026-08-10/11: the pore-size bridge, the ceiling, and the force curve

Session directive: work the silica-to-cornstarch bridge from the
inertial-dilatant-number side, fold in the new papers in Downloads, and push the
carbopol-bath, complex-substrate, and STF-on-rigid solver arms. Commits `946f013`,
`781723d`, `1418687` on `overnight/2026-06-29-autorun`. Six campaigns run and
analyzed; one theory note written (`impact_stress\PORE_BRIDGE_2026-08-10.md`, which
holds the physics synthesis).

## The night in six results

**1. The un-jamming clock is not the bridge (trsweep_2026-08-10, 18/18 cells).**
t_r was swept from 1 ms to 10 s on both calibrated presets. Impact phase: peak
jam_frac is identical to six digits across t_r = 0.1-10 s and the beta_max spread
is exactly zero, so the measured-bracket threat (master-plan item 4) is closed.
Post-impact: the mean-fhat decay rate saturates at a strain-set floor for every
t_r of 30 ms or slower; the clock is visible only below 10-30 ms, which is the
cornstarch pore-drainage scale and more than three decades from silica's. The
pre-registered falsifier fired. In flowing material the strain channel owns
relaxation, and a pore clock can only express itself in arrested material. This
reconciles the 08-08 adjudication (no a^-2 in flow-cessation release times, which
could never have seen one) with the bench phenomenology (silica fronts hold for
seconds, and those are arrested). Pore-limited creep (t_r = tau_p, going as a^-2)
survives as an arrested-state hypothesis and is now safe to wire in. Silica's
decay is not a single exponential (recoil re-jamming bump), so its decay time is
deliberately not quoted.

**2. The dilatant ceiling repairs yesterday's "50x" (dilceilsweep, 7/7 cells).**
`DIL_CEIL` is implemented and five-gated (the no-ceiling path is byte-identical to
HEAD). Sweeping the ceiling at saturated DIL_COEFF: jamming is abolished at the
face-value bulk ceiling (0.0035 rho U0^2) and fully restored by 0.011, so
CEIL_crit brackets sigma*_code = 0.0092. The pressure the dilatant term must
deliver to jam is just the onset stress, not the swept coefficient scale;
PORE_PHYSICS's 50x impossibility compared a coefficient to a pressure. The
requirement now reads K_p at 1.3-1.8x the face-value estimate. That is marginal,
and the K_p bench measurement decides it.

**3. The solver measures impact force now, and it lands on the universal curve
(FORCE_MODE + forcevalid).** F(t) = 2 pi int p y dy on the wall, with its own
0.005 D0/U0 cadence, gated byte-identical when off. The Newtonian control at L9
reads Fmax = 0.845 at t = 0.18 against the field's 0.85 +/- 0.05 at 0.17 +/- 0.02
(Cheng, Sun & Gordillo, ARFM 2022). This is the first external force anchor the
solver has hit. Amplitude is level-converged (delta Fmax = 0.015 L9 to L10); the
peak time rides a broad plateau (0.15-0.25), so quote the amplitude and flag the
time. Under the field's own t0-free protocol the early exponent is 0.495, so the
sqrt-t law is confirmed; the prefactor anchor needs the finer cadence cell now
running. Cornstarch cells: Fmax falls with U0 (0.93, 0.84, 0.80) and the
late-force fraction falls (1.03, 0.94, 0.79 of Newtonian), meaning higher U0
moves cornstarch toward the inertial side, the ID-proportional-to-U0 direction of
the regime map. U0 = 1 shows a mild late-force excess over Newtonian (+3%), the
dilatant-regime signature.
Two solver landmines were found the hard way, both now encoded in comments: an
event's schedule breaks byte-identity even with an inert body, and a non-zero
event return terminates a Basilisk run. The second killed four cells at t=0 with
rc=0 and only the instant-death check caught it; the same `return 1` idiom is the
solver's own blowup STOP a few hundred lines up.

**4. The texture solver's flat wall is Wagner-true; its grooves are not converged
(texwagner, 5/6 cells).** Brinkman wall at L10: exponent 0.485, pinned prefactor
2.445 vs sqrt6 = 2.449, agreement to 0.18%, monotone over L8-L10. The two
flat-wall debts from the 08-09 card are paid. The groove pair fails: u_solid_max
reaches 1.6% then 4.0% of U0 from L9 to L10 (the bound is 0.2%), so the leak
grows under refinement at ETA_PEN = 1e-8 under impact loading. No impalement
number is quotable. The next instrument is an ETA_PEN ladder under impact, not
more resolution.

**5. Theta does not unify yield stress and onset stress (theta_verdict on the
corrected hoop campaign).** T0 is formally VOID by the pre-registered letter (the
U0=1 anchors differ by exactly 2 cells against a declared "much less than 2";
U0=2 differs by 1). But T1 fails 0/5 with 40-91-cell separations, which is 20-45x
the anchor defect, and T2 passes 5/5: the yield-stress bath craters far less at
matched Theta, while sigma* moves the WC crater at most 5 cells across its whole
ladder. Same moral as the viscoelastic arm: the substrate classes are different
physics, not two settings of one knob. Owed: an L10 anchor pair to clear the VOID
(running now).

**6. The jam_frac factor-0.44 mystery is solved with zero compute.** The 29 JUL
theta campaign ran `WC_DRIVE=inertial` (`night_2026-07-29.sh:86`); the 08-08
rerun ran stress drive. u^2 = 4 code units against sigma*_code = 0.009 saturates
fss wherever material moves, so the two jam_frac sets were never comparable, and
all three suspects on the card (TEND, tau, SIGMA_STAR) were wrong. beta_max
reproducing while the state variable moved 2.3x is exactly what a drive-only
difference does. The cross-campaign rule is sharpened: diff the printed banners
first. Also fixed tonight: TAU_REST now enters the tag under WC_CLOCK=1, after an
18-cell sweep produced 18 identical filenames.

## Literature folded in (Downloads triage: four new plus three known)

- Mobaseri et al., PRL 136, 148201 (2026) and its SI (the ID paper, published):
  ID = rho kappa U0/(A mu_w dphi D0) = t_d/t_i with Kozeny-Carman kappa, so ID
  goes as a^2; unified force F = (1-w_F) F_i + w_F F_d; front law
  zdot_F/z_F = A mu_w dphi/(rho kappa); the force plateau peaks at ID = 1. Their
  F(t) curves (SI Figs. S10-S15) are now a direct validation target via
  FORCE_MODE. Identification: DIL_COEFF is 1/ID up to (L_d/D0)^2 eps_dil/(A dphi).
- Brown & Jaeger, J. Rheol. 56, 875 (2012): the DST stress is a boundary-confining
  stress with tau_max = 0.1 gamma/a. That is the second DIL_CEIL scale (going as
  a^-1) and the reason the free surface is a load-bearing part of STF impact
  physics.
- Royer, Blair & Hudson, PRL 116, 188301 (2016): a silica WC closure with
  phi_m = 0.592, phi_0 = 0.711, f = f_max exp(-sigma*/sigma), an external
  cross-check for si_fit.
- Cheng, Sun & Gordillo, ARFM 54 (2022) and Sun et al., Nat. Commun. 13, 1703
  (2022): the protocol checklist (universal force curve, Philippi self-similar p
  and tau, Wagner turning line 1.22 sqrt(t), benchmark case Re = 5000, We = 250).
  Several anchors are now wired into gates. Wood's law (Ultrasonics 2026): 0.5%
  entrained air cuts the bath sound speed tenfold, so any wave-speed argument must
  use the Wood c.
- None of the above contains un-jamming kinetics. The 1/t_r term remains this
  programme's genuine differentiator, now with a measured statement of when it can
  matter.

## Late-night completions

**Arrest-hold (4 cells, TEND=20): H2, an admissible pre-registered negative.**
Only one cell per preset truly arrests by t/tau = 20 (a different clock each
time), and where arrest happens the mean fhat has already decayed below the
admissibility floor. Free-running holds cannot isolate the clock; frozen-dump
restarts (the Stage 5b machinery) are the required instrument. One suggestive but
inadmissible signal: cs with t_r = 10 s holds 6x more state at tau = 20 than with
0.24 s (fhat 0.022 vs 0.0035), the right direction, to be tested properly by
restart.

**Force validation, full verdicts (`analyze_forcevalid.py`):** covered in result 3
above; the gates JSON lives with the campaign.

## night2 early returns (three of seven cells)

- **The ETA_PEN ladder REFUTES P1 at L9: the leak RISES as the penalization
  stiffens under impact loading** (u_solid_max at t/tau = 2: 0.016 at 1e-8,
  0.049 at 1e-9, 0.063 at 1e-10). The linear-in-ETA_PEN behavior of the 08-09
  static gate inverts under impact. Brinkman penalization does not converge the
  grooved-impact problem at any knob setting; the fix is a scheme change
  (embedded boundaries, e.g. Basilisk's embed.h, or solid-fraction masking with
  groove-wall refinement), not a parameter. The L10 pair will complete the
  picture but cannot rescue the method.
- **The early-force prefactor is unresolved at L9, for stated reasons**: fine
  cadence (FORCE_DT=0.001) with t0 free gives C = 2.91 vs the predicted 3.79
  (77%), and the free exponent at the earliest times is 0.32. At t = 0.004 D0/U0
  the contact patch is about 5 L9 cells and air cushioning screens the
  singularity, so the t -> 0 anchor needs L10/L11. The [0.01, 0.08] exponent
  (0.495) stands.

## night2 final returns (campaign COMPLETE 00:48)

- **Theta L10 anchors: the T0 defect HALVES under refinement.** hb 1.0117 vs wc
  1.0273 at t/tau = 1.4, a gap of 0.0156 D0 against 0.0313 at L9. The two
  solvers' Newtonian limits converge toward each other in physical units, which
  reads as a discretization difference between the code paths rather than
  physics. In L10's own quanta the gap is again exactly 2.0, so the formal
  "much less than 2 cells" bar stays marginal at every level; against the
  40-91-cell T1 separations, the T1 (no collapse) and T2 (yield-stress bath
  craters far less) verdicts stand in substance. If a formally clean T0 is ever
  needed, the residual 2-quanta offset itself wants a one-cell diagnosis (likely
  the two arms' different viscosity-branch code at fhat = 0).
- **The L10 groove pair confirms the ETA_PEN refutation**: u_solid_max = 0.043
  (1e-9) and 0.073 (1e-10) of U0, rising with stiffer penalization at L10 just
  as at L9. The Brinkman route for grooved impact is closed at both levels;
  embedded boundaries are the next scheme.

Owed beyond compute: the K_p bench measurement (decides the dilatant ceiling), an
arrested-state particle-size ladder on the bench (decides pore-limited creep), and
a frozen-dump restart campaign (decides whether the clock ever gets the stage).

Owed beyond compute: the K_p bench measurement (decides the dilatant ceiling), an
arrested-state particle-size ladder on the bench (decides pore-limited creep), and
a frozen-dump restart campaign (decides whether the clock ever gets the stage).

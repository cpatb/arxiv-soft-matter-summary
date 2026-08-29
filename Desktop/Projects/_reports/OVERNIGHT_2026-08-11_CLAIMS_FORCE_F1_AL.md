# Overnight 2026-08-11/12 — claims rewritten, force deficit diagnosed, F1 built, AL spike passed

Session branch: `worktree-sim-dev-session` in the Drop Impact repo (pushed;
based on the overnight branch tip `8c8926a`; merge with
`git merge worktree-sim-dev-session` from `overnight/2026-06-29-autorun`).
Commits `da19b5c` → (evening/overnight sequence). Status card §2 has the
condensed version; this is the full account.

## 1. The three claim rewrites (asked for, done, humanized)

`CLAIMS_REWRITE_2026-08-11.md` is now the one place claim wording lives:

- **Novelty**: opens on the empty intersection (WC field with its own
  dynamics + sharp capillary interface + free-surface impact) and the 1/t_r
  un-jamming term with measured constants; the three refuted "firsts" are
  listed with their counterexamples so they cannot creep back.
- **Stage 4**: rewritten as an onset-stress (sigma*/phi_M) ordering. Under
  the measured clock cornstarch retains more in 8/8; the old two-clocks
  wording restated its own inputs. Paper-sentence included.
- **Stage 5a**: the b-ratio gate is retired in prose AND in code —
  `analyze_stage5a_clock.py` now prints the retirement and descriptive
  statistics instead of PASS/FAIL (re-run against the 08-09 campaign;
  numbers reproduce exactly).

Living docs corrected so nothing still tells the old story:
`WC_SOLVER_METHODS.md` §2.3 (measured law now documented as the default,
legacy law demoted to reproduction-only; both old "physics wins"
re-statused — the two-clocks one withdrawn, the arrest-and-hold one marked
legacy-clock pending re-derivation), `PAPER_STF_CLAIMS_2026-08-02.md`,
master plan pointer added.

## 2. Force prefactor: the 23% deficit is three known things, not a defect

New L10 + L11 fine-cadence Newtonian cells (`forceF1_2026-08-11` ARM F) plus
night2's L9 cell, analyzed by `analyze_force_prefactor.py` (+ figure +
`wc_analysis/force_prefactor.md`):

1. **Resolution dominates**: local ratio R(t) = F/(C_pred sqrt t) rises with
   level in 7/8 bins; at t = 0.009 it reads 62% / 80% / 93% (L9/L10/L11).
2. **The fit was misspecified**: F starts at a finite ~0.2 rho U0^2 D0^2
   (air-film pressure before wetted contact) instead of falling to zero, so
   free-t0 sqrt fits absorb the offset and bias C low.
3. **The old window read the rollover**: 3.785 sqrt(t) crosses the actual
   peak value 0.85 at t = 0.050, so [0.01, 0.08] was mostly bend.

Delta→0 extrapolation (linear in Delta, stated assumption): R = 0.876,
range [0.860, 0.919] over t in [0.006, 0.03]. The residual ~12% is
consistent with the physical air cushion. The 10% gate is not passed at any
single level; what is established is the trend and the decomposition. Bonus:
at fixed t = 0.17 the despiked peak converges 0.805 / 0.807 / 0.842 upward
toward the 0.85 ± 0.05 anchor. Two data-hygiene findings: the raw probe
carries single-sample regrid spikes (despike before quoting anything), and
the L11 "mid-time structure" is correlated spike clusters (4+ consecutive
high samples around t = 0.30), not physics.

**The L12 rung COMPLETED (rc = 0, 24 min) and the linear-in-Delta check
passes with scatter**: per bin, L12 lands within 0.1-4.0% of the
three-level linear prediction, both signs, no systematic overshoot toward
C_pred (the pre-registered signature of faster-than-linear convergence).
Four-level extrapolation: R = 0.884 [0.858, 0.937] — the fourth level moves
the three-level 0.876 by less than 1%. **The ~12% residual survives
Delta -> 0.**

**And then the air-cushion attribution DIED under its own test (~00:00).**
The thin-air control was launched rather than noted (`aircush_2026-08-12`:
new ETA_AIR_FACT / RHO_AIR_FACT knobs, byte-identity gated at factor 1.0,
A1/A2 pre-registered): mu_air/10 leaves the early-force band ratio UNCHANGED
(0.822 -> 0.822) and the startup offset does not fall. Root cause found in
the source: **gap = 0 — the drop is initialised tangent to the wall** (the
GAP comment claimed "default 0.05"; the variable says 0.0; beta_wet = 0.21
at t = 0 is the discrete tangency patch). No approach phase exists in any
force record, so there was never a film to cushion. The residual and the
~0.2 startup offset are the tangent-sphere impulsive start. The GAP comment
is fixed, and the correct experiment — GAP = 0.1 with a real approach, with
and without the mu_air cut (P1/P2 pre-registered) — ran as
`gapprobe_2026-08-12` and closed the story by ~00:15:

- With a real approach the pre-contact force is 0.013 against gap=0's 0.2:
  the impulsive-start attribution is CONFIRMED, and contact lands at
  t = 0.097-0.100 as geometry predicts.
- The post-contact prefactor does NOT recover (C = 3.12, 82%, on the
  earliest clean window at L11 — same deficit as gap = 0), so the deficit is
  the under-resolved contact-line singularity itself.
- mu_air/10 still changes nothing (1-3%, within fit spread): the air film's
  physical scale is microns, far below any grid here. **Air cushioning is
  outside this solver's resolved physics at every tested configuration** —
  a limitation to state in the paper, not an explanation to reach for.

What to quote for the V2 anchor: exponent 0.495 (pass); peak 0.842 at L11,
converging upward into 0.85 +/- 0.05; prefactor as the level trend plus the
Delta->0 value 0.884 [0.858, 0.937] with the mechanism named. Not a pass.
Unexplained and parked: the rho_air/10 arm moves R DOWN (0.777) and the
startup UP (0.32).

The L12 cell only exists because of a guard bug found on its first attempt:
death at t=0 with rc=2 in 7 s — the dt-collapse guard's fixed 3e-5 floor
sits ABOVE the legitimate L12 timestep (1.26e-5). Same defect class as
detector-vs-its-own-quantum. Fixed in `wc_drop_impact.c` (floor now scales
as 3e-5 · 512/2^LEVEL), gated by byte-identity of the L10 force cell under
the new binary before the rerun was allowed to start.

## 3. Stage 7 F1 exists (`wc_analysis/fig_F1_jamfront.png` + `.md`)

si and cs at L10, measured clock, six times through impact, fields exported
from dumps by `dump2blender` (exact axisymmetric reconstruction). The figure
shows the corrected Stage 4 story in field form: both presets jam through
spreading; silica un-jams by t/tau = 0.5, transiently re-jams in the axis
collision at 1.0, and rebounds fluid; cornstarch holds jammed state through
1.0. L9 companions ran overnight and close the resolution caveat for the
impact phase (jam_frac within 0.02, beta within 0.023 at t/tau <= 0.3);
the one non-converged entry is cs at t/tau = 1.5 (0.087 vs 0.277), which
joins the theta question as a second reason to end the published F1 at
t/tau <= 1.0.

Decision left for a human: keep the retraction-phase panels (needs theta
justified from experiment) or cut at t/tau <= 0.5/1.0 and be submit-safe now.

## 4. The AL question ("best way to go about it") — answered with a plan and a passed spike

`AL_PLAN_2026-08-11.md`: scope split (AL is for the H-B family; the WC
eps_min cap is F5's problem and explicitly out of scope), ALG2 with the
plastic-term-only split so the velocity solve stays constant-coefficient,
staged gates G1–G8, a 12–16 h go/no-go spike before the 40–60 h commitment,
and a Papanastasiou continuation fallback.

The spike ran tonight and **passed both gates**
(`al_spike/al1_poiseuille.py`, `al1_gates.json`, `fig_al1.png`):

- G1: velocity to 4e-6 of the analytic Bingham Poiseuille solution across
  tau_y = 0.5/2/8; the plug is EXACT (q ≡ 0.0, not small) and its width
  lands within 0.8 cells of tau_y/G.
- G2: warm-started iterations to 1e-6: median and max 19 across a driving
  ramp, at r = mu.
- The augmentation parameter optimum is SHARP at r = mu (24 iters vs 153 /
  1135 a decade below/above): in AL2, r must track the local plastic
  viscosity per phase.
- Formulation landmine, now encoded in the script header: with the viscous
  term kept in the velocity solve, the shrinkage denominator is r alone;
  mu + r in both places double-counts mu and solves -2 mu u'' = G. The G1
  gate caught it as an exact 50% error on the first run.
- AL0 scan: the local Basilisk pin contains no prior AL/Bingham code.

**Decision: GO for AL2** (the Basilisk port), as its own focused session per
the plan's sequencing note — not started tonight by design.

## 5. Session-tooling notes (for the next agent)

- The bg-isolation worktree was cut from the repo's DEFAULT branch, months
  behind `overnight/2026-06-29-autorun`; first act had to be
  `git reset --hard overnight/2026-06-29-autorun`. If this recurs, that
  reset is step zero.
- wsl.exe mangles `$`-expansions even inside quoted payloads (bit twice:
  a variable and a `$(seq ...)`); the reliable pattern remains literal paths
  in tiny script files, per repo CLAUDE.md.
- A WSL background launch dies if the instance shuts down before nohup
  detaches — keep the launching command alive a few seconds (sleep + a
  liveness check) when no other WSL process is running.

# Overnight plan, 2026-09-27 → 06:30 2026-09-28: the σ\* axis of particle tuning, measured in a lubricated suspension

**Written 21:05 CDT, before any cell was launched.** The author asked for "the most high-impact overnight plan
that gets us closer to the heart of the particle tuning in dense suspension drop impact", devised and run
autonomously. This file is the design record (the brainstorming skill's spec, kept in the project's own
convention beside `OVERNIGHT_PLAN_2026-09-23.md`).

## Why this and not something else

The programme reduces the particle question to three numbers — x = φ/φ₀, χ = φ_m/φ₀ and σ\* — plus the
activation exponent δ (`PROGRAMME_2026-09-22.md`; exact WC covariance). Where each axis stands tonight:

| axis | state | decisive? |
|---|---|---|
| χ (window) | measured on CV, shape, doublets, chains, facets, rolling; the DNS switch (`cvmap`, `kmermap`) is ENTIRELY through χ | saturated: the critique (`CRITIQUE_AND_FUTURE_WORK_THREE_ARMS_2026-09-24.md` §2) says stop |
| x | set by weighing | not a research question |
| σ\* | dry DEM: "σ\* is contact geometry", chain predicts CV 0.2→0.5 lowers σ\* ×0.768 (11-claim ledger, 6 referee passes). **The one lubricated test (`lfsig`, 28 cells, one configuration per arm) came back INCOMPLETE as registered and, post hoc, 0.96–1.17 — never near 0.768.** | **the live contradiction** at the heart of the σ\* axis |
| δ | dry DEM: CV broadens the activation (ratio 0.869); doublets soften it under shear | never measured lubricated |
| jam-vs-thicken observable | the DNS has none; the WC jammed state is the viscosity cap by construction | a capability question, not an overnight one |

The drop is stress-controlled (`jamwin`, P2 PASS), so σ\* is where particle tuning at fixed φ acts on impact
besides χ. A dry-DEM chain that transfers to a lubricated suspension would make σ\*(PSD) predictable from
geometry; one that does not says the σ\* axis is (nearly) inert to polydispersity at fixed number-mean size —
either is a result the ledger has to carry. The `lfsig` post hoc cannot carry it: one configuration per arm
(configuration scatter unbounded), a half-decade stress ladder with the crossing always in the unsampled
3 → 10 gap, and a normalisation-dependent readout.

Ruled out for tonight, with reasons: a jam-vs-thicken solver capability (elastic jammed phase) — too large
and too risky unattended; the DEM material-axis calibration to the lab silica — no trusted rheometry target
exists (`si_fit` withdrawn); the frictionless faceted edge — the level-set build cannot converge it (card §3);
`dilsize` — 29 Basilisk cells whose every solver handle pushes the wrong way (card 2026-09-25), and the author's
own two questions on it are unanswered.

## What runs

`lfsig2` — `PREREG_LFSIG2_2026-09-27.md`. LF-DEM (Seto/Mari/Morris), the same construction as `lfsig`
(repulsive force ∝ reduced radius, λ = 0.02 a, μ 0.3, N 500, stress-controlled shear, strain 4, window [2, 4]):

- CV {0.2, 0.35, 0.5} at x_m 0.90 (Yade μ_p 0.3 φ_m convention, as `lfsig`; CV 0.35's φ_m interpolated and
  flagged), seeds {1, 2, 3}, stress ladder {1, 2, 3, 4.5, 7, 10, 15, 30, 100} σ₀ — 81 cells, ~30–60 min each,
  P = 10 on an idle 22-thread VM (load 0.00 at 20:47).
- Queue ordered so the crossing-region cells (4.5, 7, 3, 10) of every arm run first.
- Ten seed-1 cells coincide with `lfsig` cells → a bit-identity wiring gate (P4) for free.

Readouts that are DEFINED for a suspension that thins before it thickens (the `lfsig` failure): the WC
activation form fitted to the frictional-contact fraction, its midpoint, an absolute-Z crossing, the contact
share of the stress at the midpoint, and the thickening ratio. Both candidate answers (chain 0.768, inert
1.000) are registered with a resolving-power gate, and NEITHER is an admissible verdict.

## While it runs (analysis track, in order)

1. Write and dry-run `lfdem/grade_lfsig2.py` on the `lfsig` cells (reproduce their known numbers first).
2. Mechanism on the existing `lfsig` interaction dumps: the Love–Weber contact-stress split in LF-DEM and
   what the chain's geometric factor (N_c/V)⟨aᵢaⱼ⟩_c does with CV under lubrication.
3. Grade as cells land; figure L2; RESULT block; referee; Status card; commit; handoff by 06:30.

## What I will not do

Loosen a bar, add a normalisation after seeing the numbers, or quote a magnitude the seeds cannot resolve.

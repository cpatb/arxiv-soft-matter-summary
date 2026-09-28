# Paper opportunities from what is already on disk — 2026-09-27 (overnight session)

Written from two read-only surveys (vault Status cards; project folders) plus tonight's
work. "No WSL" = runs in Windows Python on data already on disk. Polydispersity jamming is
excluded (another session owns it tonight).

## Chosen for tonight

**1. bidet — the one journal-ready manuscript, blocked on one numerical defect.**
`C_P` diverges under mesh refinement (`GRID_STUDY_2026-09-03.md`), so the JFE
numerical-accuracy policy is unmet. Tonight: `PREREG_so_triplet_2026-09-27.md`, which tests
whether the first-order k/omega advection is the cause, on an r = 1.5 triplet with a
first-order control on the same grids. Status in that file and in the bidet card.

**2. dst_fit — "the Wyart–Cates strain clock depends on the protocol".** Draft in
`dst_fit/paper_strain_clock/` (J. Rheol. short communication). Built entirely on the
open Lin et al. (2015) deposit and gates that already existed. Prior-art check
(`dst_fit/PRIOR_ART_STRAIN_CLOCK_2026-09-27.md`): no one has fitted a WC strain clock to
experimental reversal data. It also found that the "0.2 from Richards' corn starch"
baseline was **misattributed** (Richards used O(0.1); 0.197 is Han et al. 2018's
start-up scale) — corrected in the gate and three docs.

## Ranked backlog (not started tonight)

| # | Project | Core claim | Blocker | Est. |
|---|---|---|---|---|
| 3 | tier0 + dst_fit | What one steady flow curve can identify in WC: `sigma*` yes (validated), `phi_0` is gauge, `phi_m` a band 6-24x the `phi` tolerance, clocks absent | none; could fold into #2 as a companion | 12-20 h, no WSL |
| 4 | Paper 1 (`front_collapse_2026-07-29`) | Impact as a rheometer: solidification onset 632 -> 479 -> 195 1/s vs salt, 3.1-3.3x invariant across constructions | **author decision #1** (goes in the PRL?) and the planned 1.5 um retake | 6 h as a section |
| 5 | onset_clocks | Stress-activation clocks cannot explain drop arrest: f moves 1-4 % while footprint moves 6.35x, 224/224 corners | uses unpublished Paper 1 data (**author decision**); negative result | 12 h, no WSL |
| 6 | CARBTRI (`Github/Drop Impact/paper_carbopol_projectiles`) | Gel crater freezes by t/tau 0.5-0.8; 3-constant drag model 2.5 % rms over 12 DNS cells | cells still to grade; sphere arm never cratered; novelty narrowed (Peng 2021, Ioannou 2023) | 8 h + WSL |
| 7 | biorheo | Newtonian R^-4 needle rule fluid-specific; published blood Casson tau_y is a fit-objective convention | **founder decision** (IP) | 10 h, no WSL |
| 8 | thrombosis_jamming | Sharp falsifiable critical shear rate, testable on published data | prior-art + data scout first; may die | 1-2 d |
| — | Slime Mold, Frisbee | manuscripts build clean | author decisions only (repository; three confirmations) | < 1 h each |

Not recommended: bloodstain_bpa (prior art fatal), the four new product folders
(design-only, no result to write), `paper_yieldstress`/`paper_stf` (HB bath numbers
carry the unsigned HOOP bias).

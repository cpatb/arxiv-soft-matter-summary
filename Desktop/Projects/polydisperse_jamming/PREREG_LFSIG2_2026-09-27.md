# PREREG `lfsig2` — does polydispersity move σ\* and δ in a LUBRICATED suspension? (registered 2026-09-27, ~20:50 CDT, before any cell was launched; queue driver started 20:56:18 CDT)

Companion: `PREREG_LFSIG_2026-09-24.md` (+ RESULT), `CLAIMS_SIGMA_GEOMETRY_2026-09-24.md`,
`derivations/sigma_star_geometry_2026-09-24.md`, `out/predict_cvmap.txt`, `out/grade_actshape.txt`.
Grader: `lfdem/grade_lfsig2.py` → `out/grade_lfsig2.txt`. Figure L2: `lfdem/make_lfsig2_figures.py`.

## Why

`lfsig` (28 cells) asked whether the dry-DEM chain's prediction — CV 0.2 → 0.5 lowers σ\* by ×0.768 at fixed
number-mean radius, through Zφ — holds with lubrication. Its registered readout (σ\*_η, the midpoint of ln η_r)
was UNDEFINED on 3 of 4 arms because a repulsive suspension shear-THINS before it thickens, and its post-hoc
contact-number readout gave 0.96–1.17 depending on the normalisation, with one configuration per arm and the
crossing always inside the unsampled σ₀ 3 → 10 gap. So the σ\* axis carries a dry-DEM prediction and a
lubricated non-answer. This campaign is designed so that every readout is defined for a thinning-then-thickening
curve, the configuration scatter is measured (3 seeds), and the crossing region is sampled.

Two things the dry DEM also predicts and nobody has measured lubricated: the activation exponent δ falls with
CV (`predict_cvmap.txt` δ ratio 0.938 at CV 0.2, 0.815 at CV 0.5 → δ(0.5)/δ(0.2) = **0.869**; from
`actshape`, Derjaguin weight, μ 0.3), and — a mechanism the `lfsig` single-seed rows suggest — the share of
the stress carried by CONTACTS at activation rises with CV (0.49 vs 0.69 at σ₀ 3), i.e. lubrication carries less
of the load in the polydisperse pack.

## Design (identical construction to `lfsig` except where stated)

- LF_DEM commit 00fa7cc0, binary md5 `5dba0579…` (~/lfdem/LF_DEM/LF_DEM/LF_DEM), `lfdem/params_base.txt`
  unchanged: repulsive force F_R e^{−h/λ} scaled by the reduced radius, λ = 0.02 a, tangential lubrication,
  Coulomb friction μ 0.3, kn = max(1e4, 1000 σ) F_R/a (the `lfsig` rule), stress-controlled shear.
- Stress unit σ₀ = F_R/a² with a = 1 = the number-mean radius (`lfdem/make_conf.py`; d10 convention).
- N = 500, lognormal PSD (project sampler, truncation 0.2–5 d10), strain 4, window [2, 4], 4 blocks.
- **Arms:** CV {0.2, 0.35, 0.5}, all at x_m = 0.90 in the Yade μ_p 0.3, p̂ 3e-5 φ_m convention:
  φ = 0.52880 (CV 0.2), 0.57300 (CV 0.5) as in `lfsig`; **CV 0.35: φ_m interpolated** by the quadratic through
  the three measured φ_m (0.57535, 0.58755, 0.63667 at CV 0, 0.2, 0.5) → φ_m 0.60749, φ = 0.54674 (linear
  interpolation would give 0.6121; the 0.0046 difference is 0.007 in x_m, ≈ 1 % in σ\* by the `lfsig`
  x-sensitivity). The CV 0.35 arm is therefore graded on SIGN/ordering only (P5), never on the chain's number.
- **Seeds {1, 2, 3}** per arm. The radius draw is seeded by the seed alone (`packing.py` convention: the same
  seed is the same normal draw rescaled), so seeds are PAIRED across CV and every CV contrast is a paired
  per-seed difference (lesson 2026-09-19/22).
- **Stress ladder** {1, 2, 3, 4.5, 7, 10, 15, 30, 100} σ₀ (9 points; the crossing region 3–10 now has 4.5 and 7).
  The thinning region below σ₀ 1 is not sampled: no registered readout uses it.
- 81 cells; queue order per seed: 4.5, 7, 3, 10, 100, 2, 15, 30, 1 (decisive cells first; the plateau cell,
  the slowest, before the outer rungs). P = 10 workers, TMAX 4 h per cell.
- Interaction dumps ON (`out_data_interaction = true`, ~1.2 MB per cell) for the post-hoc contact-stress
  decomposition; they do not affect the trajectory.

## Health rule (per cell, as `lfsig` Amendment 3)

rc ∈ {0, 124 (wall cap after the window)}, strain reaches ≥ 3.95, and max overlap −min_gap < 0.05.
An unhealthy cell EXCLUDES its arm's readout for that σ (the arm is INCOMPLETE if any ladder point is missing);
it never counts as a failure of any prediction.

## Observables (per cell, block-averaged over the window, 4 blocks, s.e. = sd/2)

η_r (col 4), η_contact/η_r, η_hydro/η_r, η_rep/η_r (cols 5, 7, 8 over col 4), Z_c = frictional contact number
(col 20; equals the contact number here because every h < 0 contact is frictional), −min_gap (col 16).

## Readouts (per arm = (CV, seed), needing all 9 ladder points)

- **R1 (primary) WC activation form.** f(σ) = Z_c(σ)/Z_c(100). Fit f = exp(−(σ\*/σ)^δ) by unweighted least
  squares in f over all 9 points, p₀ = (5, 1.5). Returns **σ\*_WC** and **δ_WC**. (On the `lfsig` seed-1 rows this
  gives 3.31/1.45 at CV 0.2 and 3.02/1.09 at CV 0.5, max residual 0.022 — the form fits.) The fit must have
  max |residual| < 0.05 or the arm's R1 is UNDEFINED.
- **R2 midpoint σ½.** The stress at which f crosses 0.5, log-linear interpolation between the bracketing rungs
  (requires f monotone non-decreasing, else UNDEFINED). This is the `lfsig` post-hoc "self" readout.
- **R3 absolute crossing σ_Z1.** The stress at which Z_c crosses 1.0 (same interpolation).
- **R4 contact share at activation.** η_contact/η_r log-linearly interpolated in σ at σ½.
- **R5 thickening ratio.** η_r(100)/η_r(1).
- **Convention factor.** R_d32 = R_d10 × (a32(0.5)/a32(0.2))² from each config's `conf.json` `a32` (the
  fixed-Sauter-diameter reading of the same data; `predict_cvmap` gave 1.145/0.735 = 1.558 for this factor).
  Printed beside R1, not graded.
- **Arm statistics.** For each readout and CV pair (0.5 vs 0.2; 0.35 vs 0.2; 0.5 vs 0.35), the paired per-seed
  ln-ratio; mean, sd, s.e. = sd/√3, paired t = mean/s.e. (df 2, t_crit 4.303 at 95 %). "16–84 % band" = mean ± s.e.

## Predictions and gates (all decided on x_m 0.90, CV 0.5 vs 0.2, R1 unless stated)

Two registered candidates for R = σ\*_WC(0.5)/σ\*_WC(0.2):
**H_geo** (the dry chain, `CLAIMS_SIGMA_GEOMETRY` C9–C11): R = 0.768. **H_inert**: R = 1.000.
Separation ln(1/0.768) = 0.264.

- **G0 resolving power.** P1 is INFORMATIVE only if s.e.(ln R) < 0.066 (a quarter of the separation, so both
  candidates cannot be consistent at 2 s.e.). Otherwise P1 prints NOT INFORMATIVE with the measured s.e., and
  no candidate is graded.
- **P1.** A candidate is CONSISTENT iff |mean ln R − ln R_cand| ≤ 2 s.e. Verdicts: GEO, INERT, **NEITHER**
  (the `lfsig` seed-1 point estimate 0.912 sits between them — a partial effect is an admissible outcome and is
  reported as the measured R with its band), or BOTH (only possible if G0 fails). I do not adjust either
  candidate after seeing the data.
- **P2 (δ).** H_δ: δ_WC falls with CV — paired ln δ(0.5)/δ(0.2) < 0 on 3/3 seeds AND mean − 2 s.e. < 0. Beside it,
  descriptive: the dry-DEM reference 0.869 is CONSISTENT iff within 2 s.e. of the mean.
- **P3 (mechanism).** The contact share at activation rises with CV: R4(0.5) > R4(0.2) on 3/3 seeds (paired).
- **P4 (wiring, bit-identity).** The ten seed-1 cells that `lfsig` already ran (CV 0.2 and 0.5 at σ₀ 1, 3, 10,
  30, 100) must reproduce `out/lfsig/<label>/data_*.dat` row-for-row, bit-identically, for every row with
  strain ≤ 4 (LF_DEM is deterministic without Brownian motion; `lfsig` ran some of them to strain 6). PASS iff
  10/10. A miss is a wiring fault and stops the grade until explained.
- **P5 (ordering, CV 0.35 arm; descriptive).** σ\*_WC(0.35) lies between σ\*_WC(0.2) and σ\*_WC(0.5) on ≥ 2/3 seeds.
  (H_geo, with ln σ\* ∝ CV² between the two measured chain points, would put R(0.35/0.2) at 0.901; printed, not graded.)
- **P6 (the old readout, for the record).** R2's ratio and its band, and R3's, are printed beside R1 so the
  `lfsig` post-hoc table (0.956 self / 1.107 absolute) can be checked against 3-seed data.

## Stated confound (registered, not corrected for)

The two arms sit at the same Yade-relative x_m but at different LF-DEM-relative distances from their own
jamming point (`lfsig` P2 inversion: φ_m,LF ≈ 0.626 / 0.66 → x_LF ≈ 0.845 / 0.87). The `lfsig` x-sensitivity is
+5 % (CV 0.2) / +9 % (CV 0.5) in σ½ per +0.05 in x_m, so a 0.025 mismatch can move R by ≤ ~4.5 % (ln 0.044),
less than the candidate separation but not necessarily less than 2 s.e. The grader prints R and, beside it,
R with this sensitivity applied as a bracket; the verdict is on the raw R. If the bracket straddles a candidate
boundary, the verdict says so.

## What would change the plan

Nothing after launch except health exclusions. If the queue driver is reaped, cells are relaunched from
scratch with the same seeds (`relaunch_lfsig2.sh`, skip-finished rule as `lfsig`), which cannot change a
deterministic cell's answer.

## Erratum 1 (2026-09-27 ~21:22 CDT, before any cell finished; found by the grader's dry run on `lfsig`)

The convention factor was registered as R_d32 = R_d10 x (a32(0.5)/a32(0.2))^2. That is wrong for this force law:
LF_DEM's repulsion is F_R0 x a0 a1/(a0+a1) (`RepulsiveForce.cpp:32`), so at fixed F_R0 (a force per unit length, the
Derjaguin form) the stress unit is sigma_0 = F_R0/a and the fixed-Sauter reading is R_d32 = R_d10 x a32(0.5)/a32(0.2),
FIRST power. The sampled, truncated configurations give a32 = 1.0589 / 1.4681 (factor 1.386), where
`predict_cvmap.txt` used the analytic untruncated (1+CV^2)^2 for CV 0.5 vs 0 (1.558) -- a provenance difference to
carry into the RESULT, not a grade. Printed only; no verdict depends on it.

## Amendment 1 (2026-09-27 ~21:35 CDT; the queue is running, no lfsig2 cell has finished, none has been read) -- a MECHANISM prediction from the `lfsig` dumps

`lfdem/mech_lfsig.py` on the eleven `lfsig` seed-1 x_m 0.90 cells with interaction dumps (`out/mech_lfsig.txt`) gives,
for CV 0.5 / CV 0.2 at matched sigma: the chain's geometric factor G = (N_c/V)<a_i a_j>_c is **1.13 at sigma_0 3,
0.84 at 10, 0.77 at 30, 0.85 at 100**, with <a_i a_j>_c 1.74x at sigma_0 3 falling to 1.50x at the plateau. So the dry
chain's ~0.77 is reproduced on the FULLY ACTIVATED network (sigma >= 10) and not at the onset, where the polydisperse
pack's first-activating contacts are large-large pairs. In a lubricated suspension the contact network is BUILT by the
stress (N_c/V rises with sigma from 0), where the dry chain evaluates its factor on a network that exists at every stress.

Registered for the lfsig2 dumps (3 seeds, 9 rungs; graded in `mech_lfsig.py lfsig2`, descriptive, not a P-gate):
- **M1 (size-selective activation).** <a_i a_j>_c(sigma)/<a_i a_j>_c(100) at sigma_0 3 exceeds 1 on 3/3 seeds at CV 0.5,
  and the excess at CV 0.5 is larger than at CV 0.2 on 3/3 seeds.
- **M2 (the chain's factor is a function of sigma, not of the pack).** The CV 0.5/0.2 ratio of G is > 1.0 at sigma_0 3 and
  < 0.90 at sigma_0 >= 10, on 3/3 seeds each.
- **M3 (what the chain would have needed).** At each arm's own sigma_half (R2), G(0.5)/G(0.2) is within +-0.15 of 1.0
  (i.e. the geometric factor does not carry the 0.768 at the onset). Reported with its seed band.
If M1-M3 hold, the reading is: the dry chain's Z*phi mechanism is real but describes the plateau network, and the onset
stress is set by a size-selected sub-network whose geometric factor is CV-neutral -- which is why sigma* barely moves
(P1) while Z_c(100) falls x0.8 with CV. If they fail, the discrepancy is unexplained and is reported as such.

## Amendment 2 (2026-09-27 ~21:45 CDT; no lfsig2 cell has finished or been read) -- the stated confound gets a MEASURED arm

The "stated confound" section argues the x_LF mismatch (CV 0.5 arm ~0.025 closer to its own LF-DEM jamming point)
moves R by at most ~4.5 % from the lfsig single-seed x-sensitivity. An argued bound is weaker than a measured one, so a
fourth arm is added beside the running queue (second driver, P = 5; `lfdem/make_queue_lfsig2x.py`, `run_cell2x.sh`,
`launch_lfsig2x.sh`, run set `lfsig2x`, 15 cells):
- **CV 0.5 at the CV 0.2 arm's LF-relative x**: x_LF(0.2) = 0.52880/0.6258 = 0.84500 (lfsig P2 inversion, x_m 0.90
  arms), phi_x = 0.84500 x 0.6585 = **0.55643** (registered arm 0.57300; Yade x_m 0.874). Seeds 1-3, rungs
  {3, 4.5, 7, 10, 100} only (the crossing region plus the normalising plateau).
- **Readouts on this arm**: R2 (sigma_half, 5-rung log-linear) and R1 refitted on the SAME five rungs for the 0.2, 0.5
  and 0.5x arms (like for like; the 9-rung R1 stays the primary). Paired per seed as before.
- **X1 (registered).** ln[R2(0.5x)/R2(0.5)] -- the confound's own size -- is reported with its seed band. Prediction from
  the lfsig sensitivity: **-0.047 (ratio 0.954)**. The verdict of P1 stays on the raw 9-rung arms; X1 says whether the
  raw R could be moved across a candidate boundary by the mismatch: if |X1| < the P1 margin to the nearest candidate
  boundary, the confound is CLOSED; otherwise the P1 verdict is reported as CONDITIONAL on it.
- **X2 (descriptive).** R on the x-matched pair, 5-rung R1 and R2, printed beside the registered R.
The phi_m,LF values come from the lfsig inversion that failed its own two-phi consistency check for CV 0.2 (P2 NOT USABLE,
0.6945/0.6810 for phi_0 -- the phi_m halves agreed, 0.6258/0.6254); the x-match is therefore approximate to the ~0.5 %
those halves differ by, which is far below the 2.9 % shift being tested.

## Amendment 3 (2026-09-27 ~21:50 CDT; physics-referee pass 1 on the registration, acted on before any cell finished or was read)

- **Coverage (referee MATERIAL 1).** The "2 s.e." bands in G0/P1/P2 cover ~82 % at df 2 (t.cdf(2, 2)), not the 95 % the
  "t_crit 4.303" bullet implied, and T_CRIT was never used. The REGISTERED verdicts are now at 95 %: G0 requires
  se(ln R) < ln(1/0.768)/(2 x 4.303) = **0.0307** (so the two candidates' 95 % intervals cannot both contain the mean), and
  P1 CONSISTENT means within 4.303 s.e. The 2-s.e. (82 %) reading is printed beside it as descriptive, never as the verdict.
  X1's margin uses the same 95 % boundaries. This is STRICTER than what was registered, not looser.
- **P2's side.** "mean - 2 s.e. < 0" was the wrong side for a predicted FALL (a fall is resolved when the upper end is below 0).
  Registered now: 3/3 seeds negative AND mean + 4.303 s.e. < 0. The as-written test and the 82 % one are printed beside it.
- **Size-blind kn (referee MATERIAL 2).** kn is one spring per cell while a pair's activation force scales with R*
  (x~25 across the PSD), so large-large contacts overlap more at the same load per threshold. M1 (size-selective activation)
  is graded only beside an overlap check CONDITIONED ON R*: `mech_lfsig.py` prints, per cell, the mean and max -h/R* of
  touching contacts in the top and bottom R* quartile. If the top quartile's max -h/R* exceeds 0.05 at CV 0.5 on any rung M1
  uses, M1 is reported CONDITIONAL on the contact stiffness rather than as evidence for size selection.
- **Minor, fixed:** the loader tolerates a partial last row of a live data file; M3's interpolation refuses to extrapolate
  (np.interp clamps); the identity-check range is quoted as measured (0.76-1.02; 0.88-1.02 at sigma >= 3).
- **Unresolved (referee minor 1), carried into the RESULT:** the Design section and params_base.txt call the stress unit
  F_R/a^2 while Erratum 1 argues F_R0/a. LF_DEM's pair force is F_R0 x R*, so at a = 1 the two coincide numerically; only the
  printed-not-graded R_d32 factor depends on the exponent.

## Amendment 4 (2026-09-27 ~22:15 CDT; registered before any of its cells exist; lfsig2 cells so far read only for health and P4)

A **monodisperse arm (CV 0)** gives the onset question its longest lever arm: the dry chain predicts
sigma*(0.5)/sigma*(0) = **0.735** and sigma*(0.2)/sigma*(0) = 0.957 (`out/predict_cvmap.txt`, d10). Run set `lfsig2m`
(`lfdem/run_cell2m.sh`, `launch_lfsig2m.sh`, `queue_lfsig2m.txt`), launched when the x-matched arm's five workers free up.
- CV 0 at x_m 0.90: phi = 0.9 x 0.57535 = **0.51782**; seeds 1-3 (positions differ, radii identical); rungs {3, 4.5, 7, 10, 100};
  otherwise identical to lfsig2. 15 cells.
- **Readouts (like for like, five rungs):** 5-rung R1 (sigma*_WC, delta_WC) and R2 on the CV 0 arm, against the same 5-rung
  refits of the CV 0.2 and 0.5 arms already printed by the Amendment 2 block.
- **Q1 (registered).** R0 = sigma*_WC(0.5)/sigma*_WC(0), paired by seed index, 95 % (t_crit 4.303). Candidates: GEO 0.735,
  INERT 1.000 (separation ln 0.308). Informative iff se(ln R0) < ln(1/0.735)/(2 x 4.303) = 0.0358. Verdict GEO / INERT /
  NEITHER / NOT INFORMATIVE, as P1.
- **Ordering confound (stated, and checked).** Monodisperse spheres can crystallise under shear near phi 0.5-0.6; an ordered
  arm would not be the random reference the chain assumes. Check: the global bond-orientational order Q6 (Steinhardt;
  neighbours within 1.2 x contact distance) on each cell's last configuration snapshot, against the same number on the CV 0.2
  arm. If any CV 0 cell has Q6 > 0.3 (random packings sit near 0.03-0.05 for N = 500; fcc is 0.575), that cell is excluded
  from Q1 and reported as ORDERED.

## Amendment 5 (2026-09-27 ~22:40 CDT) -- per-grain activation threshold, found POST HOC on lfsig, registered for lfsig2

**Status of the data at registration:** the lfsig2 dumps synced so far have been read ONLY by the health/P4 grader (data
rows, not interaction dumps). `lfdem/sizeact.py` has been run on the OLD lfsig seed-1 cells only (`out/sizeact_lfsig.txt`).

**The post-hoc finding (lfsig, one seed per arm, descriptive).** Activation is size-selective at onset: at CV 0.5, sigma_0 1,
the largest radius quartile is 12.8x more likely than the smallest to carry a touching contact, falling to 1.9x at sigma_0 100;
the per-quartile mean-degree exponent p falls from ~2.5 at onset to ~1.4 on the plateau. Normalising each quartile's
activated fraction by its own plateau value (small grains rattle even on the plateau), the activation collapses best onto
sigma * a^m with **m = 1.05 (CV 0.2) and 1.50 (CV 0.5)**, against a size-blind m = 0 whose residual is 1.2-2.0x larger.
The Derjaguin argument predicts **m = 1**: the load a stress puts on a contact scales as sigma a^2 and its repulsive threshold
as R* ~ a, so a grain activates when sigma a exceeds a constant.

**Registered on the lfsig2 dumps (3 seeds, 9 rungs; `python lfdem/sizeact.py lfsig2`, collapse block; the grid and the
plateau normalisation exactly as in the file at this commit):**
- **S1 (size selectivity).** The pooled best m lies in [0.5, 2.0], and the residual at m = 0 exceeds the best residual by
  >= 1.3x, on the pooled set AND on each of CV 0.2 and CV 0.5 separately.
- **S2 (Derjaguin).** m = 1 is within 1.15x of the best residual on each CV (i.e. the one-parameter Derjaguin law is nearly
  as good as the free exponent). If S1 passes and S2 fails, the activation is size-selective but NOT by the simple
  Derjaguin load/threshold ratio, and that is reported.
- **S3 (does it explain sigma*?, descriptive).** From each arm's measured contact-weighted radius a_k = sum(k a)/sum(k) at
  sigma_half, the per-grain law predicts sigma*(0.5)/sigma*(0.2) = (a_k(0.2)/a_k(0.5))^m_best. Printed beside the measured
  R (P1); agreement within the 95 % band would make the per-grain threshold the mechanism of whatever P1 finds.

## Event (2026-09-27 ~23:10 CDT, execution only)

Throughput check: 17 of 81 main cells done at 22:58, cells 26-75 min (mean ~45), the VM shared with another session's five
simpleFoam jobs. To finish all three seeds before morning, `run_cell2.sh` gained an ATOMIC claim (mkdir without -p; a cell
whose directory exists is skipped), deployed by rename so running cells kept their copy (tested on a finished and a running
cell: both skipped, nothing touched), and a reverse-order helper (`lfdem/helper_lfsig2.sh`, P 5) is chained to start when the
CV 0 arm finishes. Cells are deterministic, so which driver runs a cell cannot change its answer; P4 still checks that.

## Amendment 6 (2026-09-27 ~23:20 CDT; lfsig2 interaction dumps still unread) -- prior art found; the exponent is what is testable

A literature check (`LIT_SIGMA_CV_2026-09-27.md`; verified by me against arXiv:1901.02066v2) found that **large-first
frictional activation in bidisperse suspensions is published**: Guy, Ness, Hermes, Sawiak, Sun & Poon, Soft Matter 16, 229
(2020), doi:10.1039/c9sm00041k -- "frictional contacts first form amongst the large species", small particles "largely
redundant for stress transmission". Amendment 5's size selectivity (S1) is therefore a CONFIRMATION in a continuous PSD,
not a finding. Their contact model takes the critical load F* "independent of the size of the contacting particles", which
puts a grain's activation stress at F*/a^2: **m = 2** in Amendment 5's notation (their sigma*_1 ~ 0.06 sigma_0 vs
sigma*_2 ~ sigma_0 at size ratio 1/0.25 = 4 is a factor ~16 = 4^2). LF_DEM's repulsion here is Derjaguin-scaled (threshold
~ R*), which predicts **m = 1**. The exponent is the discriminating, and new, quantity.

- **S4 (registered).** On each of CV 0.2 and CV 0.5 separately, rms(m = 1) < rms(m = 2), and on the pooled set
  rms(m = 2)/rms(m = 1) >= 1.15 (Derjaguin beats constant-force). Verdicts: DERJAGUIN (all three), CONSTANT-FORCE (m = 2
  better on both CVs and pooled ratio <= 1/1.15), else MIXED.
- **Prior information, stated honestly:** on the old lfsig seed-1 cells the same comparison gives rms(m 1)/rms(m 2) =
  0.557/0.629 at CV 0.2 (m = 1 better) and 0.583/0.561 at CV 0.5 (m = 2 slightly better), best m 1.05 and 1.50. So a MIXED
  outcome is live, and it is reported as such if it happens.
- Claim 10 of `CLAIMS_SIGMA_GEOMETRY` is rescoped the same hour (Guy 2020 is prior art for the mechanism; the surviving
  novelty is a continuous PSD, a Derjaguin threshold, and the exponent).

## Interim note (2026-09-28 ~00:10 CDT; seed 1 complete, seeds 2-3 running; NOT a grade)

Seed 1 alone (no verdicts are drawn from n = 1): R(sigma*_WC, 0.5/0.2) = 0.951; delta ratio 0.793; contact share at
activation x1.154; x-matched arm moves sigma_half by ln -0.055 (predicted -0.047). Grain-level best m climbs with CV
(1.10 / 1.35 / 1.70), S4 reads MIXED so far; S3's per-grain prediction (0.672) over-predicts the measured onset shift.
POST HOC, labelled (`lfdem/forcesize.py` -> `out/forcesize_lfsig2.txt`): the pair LOAD scales mean-field, F ~ R*^beta with
beta 2.05 / 1.76 / 1.95 (not-yet-touching pairs, sigma <= 10) at CV 0.2 / 0.35 / 0.5, so load routing does NOT explain the
rising grain-level m; and a PER-PAIR activation collapse (no neighbour multiplicity) cannot discriminate m = 1 from m = 2
here -- rms 0.418 vs 0.404 pooled, a 3 % difference against the 15 % the registered S4 uses. The exponent question is
therefore likely to end MIXED or NOT INFORMATIVE on this design; a wider size ratio (bidisperse 1:4) is what would separate it.

## Amendment 7 (2026-09-28 ~00:30 CDT; registered before any of its cells exist) -- a binary arm that can decide the exponent

The continuous-PSD exponent test (Am. 5/6) is likely to end MIXED because the pairs' reduced radii span too little (interim
note). A BINARY suspension at size ratio 2 has three discrete pair types, R* = 0.5 (SS), 2/3 (SL), 1 (LL), and the two
threshold laws predict activation stresses in the ratio sigma_half(SS)/sigma_half(LL) = **2** (Derjaguin: threshold ~ R*,
load ~ sigma R*^2, the load exponent measured 1.8-2.1 in `out/forcesize_lfsig2.txt`) against **4** (size-independent
critical load, Guy et al. 2020).
- **Design.** Radii 1 and 2 (deterministic, 444 small + 56 large; small species 49.8 % of the solid volume), phi 0.55, mu 0.3,
  N 500, stress ladder {1, 1.5, 2, 3, 4.5, 7, 10, 15, 100} sigma_0 (sigma_0 = F_R0/a_s^2, a_s = 1), strain 4, window [2, 4],
  seeds 1-3 (positions only): 27 cells, run set `lfsig2bi` (`make_conf_bi.py`, `run_cell_bi.sh`, `launch_lfsig2bi.sh`),
  launched when the main lfsig2 queue completes. phi 0.55 is a judgement (no measured phi_m for this binary; the readout
  is a ratio between pair types inside one suspension, so it does not need a matched x).
- **Readout and gate** exactly as coded in `lfdem/pairtype.py` at this commit: per pair type, the touching fraction of listed
  pairs, normalised by its sigma-100 plateau; sigma_half by log-linear interpolation. **B1 (registered, 95 %):** ln ratio
  SS/LL, INFORMATIVE iff se < (ln 4 - ln 2)/(2 x 4.303) = 0.0805; DERJAGUIN / CONSTANT-FORCE / BOTH / NEITHER as coded.
  SL/LL (1.5 vs 2.25) printed as a consistency check.
- **Why it matters.** Guy et al.'s mechanism (large pairs activate first) is prior art; its SIZE LAW is not tested for a
  Derjaguin (electrostatic) repulsion, which is what the lab's silica has. A DERJAGUIN verdict halves the log-selectivity
  relative to the constant-force picture and is the part of Paper C that would be new.

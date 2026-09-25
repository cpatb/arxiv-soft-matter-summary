# PREREG — `dilsize`: can the drop solver carry the size half of the ball → hybrid transition?

Written 2026-09-25, before any cell (commit time is the record). Scaling estimate:
`derivations/ball_hybrid_scaling_2026-09-25.md` (numbers in `out/ball_hybrid_scaling.txt`).
**Not launched.** It runs in `basilisk-drop-impact`'s campaign tree (another session's workspace) and needs the
author's OK after tonight's harness memory reap (PREREG_LFSIG Event).

## The question, and why it is framed this way
The bench says: at similar chemistry, larger particles and wider size distributions move an impacting drop from a
solid, rebounding ball to an inertia-dilatant hybrid. The project has shown the CV half in the DNS through the
jamming window (`cvmap`). The size half has only been argued, through permeability (ID ∝ κ ∝ a²).

The scaling estimate finds two things that shape this registration.
- **(S1)** In `wc_drop_impact.c` the dilatant pressure enters only as `sigma += p_dil`, the stress that RECRUITS
  frictional contacts (lines 1119–1133). It does not change the jammed state's viscosity or strength.
- **(S2)** Along a size ladder, every size handle the solver has — σ\* ∝ a⁻¹ or a⁻², and the dilatant activation
  p_dil — makes contacts easier to activate as a grows. So the solver is expected to move TOWARD the ball with
  size, i.e. the wrong way.

This registration therefore tests whether the solver can produce the size trend at all, and maps what its dilatancy
term does. It is designed so that a failure is informative: if P3 holds, the size trend needs physics the closure
lacks, and §6 names the minimal candidate.

## Base cell and fixed settings
`cvmap`'s in-window monodisperse cell `c0.0_x0.932`: PHI 0.57784, PHI0 0.62, PHIM 0.56459, DELTA 1.3169,
SIGMA_STAR 300 Pa, cs_fit preset, WC_CLOCK 1, THETA0 5, L9, TEND 6, CFL 0.10, SHAPE_MODE 1, **EPS_MIN 0.03** (the
cap at which `cvmapE03` had no deaths), JAM_FHAT = 0.7609 (this cell's f_c, as in `cvmapE03`; unchanged by σ\* and
DIL), U₀ = 4 m/s (σ\*_code = σ\*/ρU₀²; banner-printed, never assumed).
Binary: the G0-gated `jamwin` build (md5 db413e87dd21), which already carries DIL_COEFF/DIL_CEIL.

## Observables (from `spread_*.dat`; all DIL settings verified from the banner line `# DIL_COEFF=`)
- J_max: maximum jammed volume fraction (col 7 at JAM_FHAT = f_c); Rg_peak: maximum of col 11; t_95: first
  t/τ at which Rg reaches 0.95 Rg_peak; Rg_end at TEND.
- **Class** (fixed now): FLUID if J_max ≤ 0.01 (the `kmermap`/`cvmap` J_ON); **BALL** if J_max > 0.01 and
  t_95 ≤ 0.5; **HYBRID** if J_max > 0.01 and t_95 > 0.5 (jams, yet keeps spreading well past impact).
  ⚠ The project has already found that no DNS observable separates a jammed drop from a near-edge thickened one
  (card, 2026-09-24); the class is a solver-internal reading, not a bench signature.
- Health: rc 0, t/τ reached TEND, no energy injection (late umax ≤ 2 U₀ and CoM below release, the admissibility
  gate added after `cvmap`'s 27th cell). Unhealthy cells are excluded, never graded as a class.

## Arms (21 cells)
**Arm A — does dilatancy do anything where shear already activates? (4 cells)** Base σ\* 300 Pa,
DIL_COEFF ∈ {0, 0.1, 1, 10}, no ceiling. For scale: the face-value microscopic silica p_dil (~3×10⁵ Pa at γ̇ ~ U₀/D₀,
the solver's own DIL_COEFF comment) is DIL_COEFF ≈ 15 at ρU₀² ≈ 1.9×10⁴ Pa, so the range brackets it from below.

**Arm B — can dilatancy activate where shear cannot? (9 cells)** σ\* raised ×30 to 9000 Pa (so the shear route is
starved), DIL_COEFF ∈ {0, 0.01, 0.1, 1, 10}, no ceiling; plus DIL_COEFF ∈ {1, 10} with DIL_CEIL = 0.0035 (the
silica S at φ 0.44) and with DIL_CEIL = 4 (the cornstarch S).

**Arm C — the size ladder through the solver's own handles (8 cells; the a₀ cells carry a ceiling and so are new
cells, not Arm A's).** a/a₀ ∈ {0.5, 1, 2, 4} at fixed (PHI, PHIM, DELTA), U₀ 4:
- [D] Derjaguin: σ\* = 300·(a/a₀)⁻¹ Pa, DIL_COEFF = 1·(a/a₀)⁻², DIL_CEIL = 0.0035 (S is size-independent in [D]).
- [F] fixed repulsive force: σ\* = 300·(a/a₀)⁻² Pa, DIL_COEFF = 1·(a/a₀)⁻², DIL_CEIL = 0.0035·(a/a₀)⁻¹.
These are the two conventions `ball_hybrid_scaling` prints; neither is claimed to be the material's.

## Predictions and gates (graded by a grader committed before the first cell is read)
- **G0 (wiring).** Arm A's DIL_COEFF 0 cell reproduces `cvmapE03`'s `c0.0_x0.932_U4.0` `spread.dat` byte for byte
  (same binary, same settings; the reference exists: `out/cvmapE03/c0.0_x0.932_U4.0`).
- **P1 (inertness where saturated).** Arm A: every healthy cell has the same class as DIL_COEFF 0, and
  |J_max − J_max(0)| ≤ 0.1·J_max(0). Reason: σ_shear/σ\* ≫ 1 already; exp(−(σ\*/σ)^δ) is at 1 either way.
- **P2 (activation threshold).** Arm B, no ceiling: the DIL_COEFF 0 cell is FLUID or has J_max ≤ ½ Arm A's; the class
  switches to BALL or HYBRID at the first DIL_COEFF with Π_d ≡ DIL_COEFF·F(φ)/σ\*_code ≥ 1 (γ̇_ref = 1 in code
  units, F(φ) from the banner), and not below Π_d = 0.3. **P2b:** with DIL_CEIL 0.0035, Π_d ≤ 0.0035/σ\*_code < 1,
  so no switch at any DIL_COEFF; with DIL_CEIL 4, the switch happens as without a ceiling.
- **P3 (the registered NEGATIVE: the solver's size handles go the wrong way).** Arm C, both conventions: J_max is
  non-decreasing and Rg_peak non-increasing in a/a₀ (tolerance 2 % of the a₀ value), and no class moves from BALL
  toward HYBRID or FLUID as a grows. **PASS means the WC closure, dilatancy included, cannot produce the bench's
  size trend.** FAIL (a larger-a cell is more HYBRID/FLUID) would mean the solver has a size mechanism the scaling
  estimate missed, and that cell is examined before anything else is said.
- **P4 (descriptive).** Arm C's Rg_peak and t_95 are printed against Π_d(a) and σ\*_code(a) to show which handle
  moved them.

## What would change the plan
- If more than ⅓ of cells die at EPS_MIN 0.03, the arm is re-registered at L8 before any class is read.
- If G0 fails, nothing is graded until the difference is explained.

## 6. If P3 passes: the minimal missing term (not part of this registration)
The scaling estimate places cornstarch where drainage, dilatant strength and inertia are all O(1) (De ≈ 7, S ≈ 4)
and silica where the pack is undrained with a weak network (De ≈ 7×10⁵, S ≈ 0.004). A closure that can carry that
difference needs the jammed state's stress capacity to depend on drainage, e.g. cap the jammed-region stress at
μ_eff·min(p_Darcy(De), K_pΔφ/φ) instead of the fixed viscosity cap. It would be a new solver term with its own
reduction gate and pre-registration. The bench measurements that would pin its inputs are K_p for the lab silica
(the λ_D series) and whether the observed size trend is at fixed chemistry.

## Amendment 1 (2026-09-25, pre-data referee acted on; no cell has run)
- **F1 — class rule replaced.** A first-crossing t_95 called a drop that jams, UNJAMS and oscillates a BALL
  (the referee's on-disk counterexample, `out/jamwin/N_s300_d1.3169`: J 0.318 → 0 by t/τ 0.3, Rg retracting to 46 %
  of its peak and oscillating). The rule is now `dilsize_classify.py` (committed with this amendment):
  FLUID J_max ≤ 0.01; **TRANSIENT** J_max > 0.01 but J_end ≤ 0.01 (mean over the last 0.5 τ); BALL J_end > 0.01 and
  Rg_end ≥ 0.9 Rg_peak and t_95 ≤ 0.5; HYBRID J_end > 0.01 otherwise. Checked before any dilsize cell exists, on
  on-disk cells: the counterexample → TRANSIENT; `cvmapE03` c0.0_x0.932_U4.0 and c0.0_x0.950_U2.0 → BALL;
  c0.5_x0.932_U4.0 and c0.5_x0.920_U2.0 → FLUID; c0.2_x0.932_U4.0 → TRANSIENT (jams, then unjams). ⚠ BALL cells'
  J_end is small (0.031, 0.061: the rest clock relaxes f̂ after arrest), so J_end > 0.01 has a factor ~3–6 of
  margin; every Arm C cell's full J and Rg series is plotted and looked at before its class is quoted.
  P1–P3 are restated on these classes; P3's "moves from BALL toward HYBRID or FLUID" now reads "toward HYBRID,
  TRANSIENT or FLUID".
- **M2 — health gate corrected to the project's own** (`grade_capvol.py`): after t/τ 0.5, umax ≤ U₀, zc ≤ 0.5 D₀
  and h ≤ D₀ (the registration had typed "≤ 2 U₀" and dropped h). In `dilsize_classify.py`.
- **M1 — a DIL-off control ladder (Arm C0, +8 cells, 29 total).** Arm C varies σ\* and DIL_COEFF together, so a
  P3 PASS could be σ\*'s size scaling alone. Arm C0 repeats Arm C's σ\*(a) schedule with DIL_COEFF = 0 and no
  ceiling. P3 is graded on Arm C as registered AND reported per arm; the dilatancy contribution is Arm C minus
  Arm C0, cell by cell. **The registration's claim is restated as "the WC closure's size handles (σ\* and the
  dilatant activation term, separately and together)", not "dilatancy".**
- m4: Arm C holds PHI, PHI0, PHIM fixed, so the loss of (φ, φ₀, φ_m) covariance with DIL_COEFF > 0 (PHI_REF is
  fixed at 0.49) never enters; extending the ladder in φ would need it handled.

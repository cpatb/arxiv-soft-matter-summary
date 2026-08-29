# Session 3, 2026-08-09 — dilatancy implemented, and a bug in five solvers

Follows `OVERNIGHT_2026-08-09b_STAGE_RERUN.md`. Read §1 and §5.

---

## 1. The short version

Three things, in descending order of consequence.

1. **Reynolds dilatancy is implemented and it makes a measurement, not an assertion.** The
   permeability factor genuinely controls *where* dilatancy bites: the amplitude at which
   the drop tips into jamming falls with volume fraction, in the direction and roughly the
   magnitude Kozeny–Carman requires. **This is the pore-size physics you pointed at,
   finally doing testable work.**
2. **The axisymmetric hoop term was wrong on every radial face, in five solvers**, three of
   them unconditionally. It computed the *axial* velocity over the radius instead of the
   radial one, and the error **grows** under refinement. Fixed and gated.
3. **Stage 4's claim and Stage 5a's gate are rewritten in the master plan**, and the
   novelty paragraph is replaced — three claims it made are refuted by prior work.

---

## 2. Reynolds dilatancy

**The mechanism.** To un-jam, the frictional network must dilate. The drop is
incompressible, so the new pore volume must be fed by liquid dragged through the packing.
That Darcy flow is resisted by the permeability, and the resulting particle pressure **adds
to the stress that recruits frictional contacts**:

```
p_dil ~ (eta_f / k) * eps_dil * gammadot * L_d^2,   k = (1-phi)^3 a^2 / (180 phi^2)
```

This is referee objection #1 for the whole DNS argument — Athani, Metzger, Forterre & Mari
(PRF 10, 043301, 2025) show the stress in a locally shear-jammed region is set by this and
not by the flow rule, and their closing sentence names impact as the target. It is *also*
exactly where the 08 AUG pore-size adjudication said permeability belongs: **setting a
stress, not the un-jamming clock.** Two independent lines of reasoning, same missing term.

**One swept amplitude, not four microscopic parameters**, and the reason is arithmetic.
Evaluated at face value with silica numbers the expression gives p_dil ~ 3×10⁵ Pa against
σ\* = 44 Pa — a factor of 10⁴, which would jam the whole drop at t = 0. That is the same
order-of-magnitude problem already visible in Mobaseri et al.'s published dilatancy number.
Choosing four uncertain quantities that multiply out to whatever is wanted would not be a
model, so the solver exposes `DIL_COEFF` and `F(phi)` and the sweep does the work.

**Four gates pass** (`dilatancy_gate.sh`): reduction, effect, sign (monotone — dilatant
pressure must *raise* jamming), and the permeability factor checked against arithmetic done
outside the solver.

> ⚠ **Gate B initially failed and the gate was wrong, not the term.** It asserted on
> `jam_frac`, which is a *thresholded* volume fraction and saturates at exactly 1 wherever
> the drop jams through — so an inert term and a saturated observable are indistinguishable
> on it. Re-asserted on the continuous `fhat_max`, at a φ where the system is marginal. The
> gate now also *prints* when `jam_frac` is saturated, so a future reader cannot repeat it.

### 2.1 The result: the permeability factor sets the threshold

`DIL_crit` is defined **before the runs**, in the launcher, as the geometric midpoint of the
rungs bracketing `jam_frac = 0.5` — interpolation-free, so it cannot be tuned afterwards.

| φ | F(φ) | DIL_crit | bracketed by |
|---|---|---|---|
| 0.40 | 0.4092 | **1.212** | 1.0 – 1.47 |
| 0.44 | 0.6091 | **0.562** | 0.464 – 0.681 |
| 0.54 | 1.6551 | *running (stiff)* | — |

**Measured ratio 2.157 against Kozeny–Carman's predicted 1.488.**

The **direction is unambiguous**: a denser packing jams at a *lower* dilatant amplitude,
because the pores are tighter and the backflow is harder. That is the pore-size physics
producing a falsifiable prediction and surviving it — and it needed **no knowledge of the
absolute amplitude**, which is the one quantity we cannot currently justify.

The magnitude is 45% steeper than the bare permeability ratio, right at the edge of what an
11-rung ladder at 6 rungs/decade resolves (1.47×). So it passes, **marginally**, and should
not be over-read as a confirmed exponent.

⚠ **The deviation has a physical explanation the gate does not model, and its sign is the
expected one.** The competing viscous stress `eta_wc(fhat)*gammadot` carries its own
φ-dependence in the *opposite* direction, so the threshold ratio *should* exceed the bare
permeability ratio. It does. That is a consistency check passing, not a fit.

⚠ **φ = 0.54 is the stiff end and its cells log `u.x` convergence failures.** Raising the
dilatant stress drives `fhat` → 1, the viscosity to its ceiling, and the implicit solve into
the regime that is already this suite's binding numerical constraint. Worth recording: the
new physics and the known numerical weakness meet in the same corner.

---

## 3. The hoop term was wrong on every radial face

**Confirmed from `qcc -source`, not inferred.** `D_tt = u_r/r` was written as
`0.5*(u.y[] + u.y[-1,0])/y` inside a **generic** `foreach_face()`. That loop iterates both
face families and **rotates field components but not coordinates**:

```
x-faces:  ur = 0.5*(val(u.y,0,0,0) + val(u.y,-1,0,0))    correct
y-faces:  ur = 0.5*(val(u.x,0,0,0) + val(u.x,0,-1,0))    the AXIAL velocity
```

The other three terms are immune because `2(D11²+D22²)+4D12²` is symmetric under the 2-D
component swap. `D33` is the only asymmetric one, which is why it survived.

**Why it is worse than a typo.** For pure axial shear `u_z = c·y`, `u_r = 0`, the true hoop
term is zero and the buggy one returns `c` — inflating γ̇ from `c` to `√3 c`, a 73%
over-estimate. And the spurious term is `u_z/r ~ U0/Δ` near the axis, so it **grows under
refinement**: a defect that gets *worse* on better grids, in the opposite direction to the
correction `HOOP` exists to make. The existing `hoop_gate.sh` cannot catch it — its sign
test passes for the wrong reason.

**The fix** is one idiom once the problem is named: build the hoop rate as a cell-centred
**scalar** `d33c = u.y[]/y`, where `u.y` unambiguously means the radial component, then
average onto faces. A scalar's *stencil* rotates with the loop — which is what we want, it
picks the two cells straddling whichever face we're on — while its *value* does not.

**Five solvers affected.** `hb_drop_on_bath.c` and `wc_drop_on_bath.c` had it gated behind
`HOOP` (default off, so only campaigns that set `HOOP=1` are affected — **`theta.sh`'s
campaign did**). `rheo_substrate.c`, `sphere_impact_hb.c` and `sphere_settle_hb.c` carried
it **unconditionally with no knob**, so they have been computing it wrongly for their whole
existence. All five now default to the corrected form with `HOOPFACE=0` restoring the old
one for reproducing earlier runs. All 13 solvers build clean.

Gates at L8 on `hb_drop_on_bath.c`: **A, B, C all PASS** (HOOP=0 byte-identical;
HOOPFACE=0 reproduces the old behaviour exactly; the correction changes the answer). L9 and
Gate D — *does the discrepancy grow with refinement, as a 1/Δ term must* — were still
running at the time of writing.

---

## 4. Master plan rewritten

- **Stage 4's claim** now leads with the fact that the ordering is a σ\*/φ_M asymmetry, not
  a clock asymmetry, and points the *other* way. The legacy numbers stay as the control.
  It also tells the reader not to cite `fig_stage4_sics_twoclocks.png`, which draws the
  withdrawn attribution **in pixels**, where no macro guard can see it.
- **Stage 2's "KEY PHYSICS WIN"** is struck through, and it is the more embarrassing of the
  two: *"cs melts on its 20 ms clock, si holds flat on its 2 s clock"* is not a finding, it
  is a restatement of two numbers that were typed in.
- **Stage 5a's b-ratio gate is retired**, not re-tuned. It fails by 5.09σ and was already
  weak three ways.
- **The novelty paragraph is replaced.** Three claims are refuted (Baumgarten–Kamrin,
  Angerman, Athani each got somewhere first). What survives is the *intersection*, plus the
  `1/t_r` term — both prior WC implementations use a pure strain clock and so **cannot
  un-jam material that has stopped moving.**

---

## 5. What to do next

**Compute, in order:**
1. Finish φ = 0.54 for the permeability test — it is the long lever (F = 1.66 against 0.41,
   a 4× span) and it is the stiff end, so it may need `MUMAX_WC` or a finer CFL rather than
   just more time.
2. **Gate D on the hoop fix.** If the discrepancy grows L8→L9 the diagnosis is confirmed
   and every `HOOP=1` campaign needs re-running; if it does not, the term is still wrong
   (proved from source) but that observable does not expose it.
3. **Re-run `theta.sh`'s campaign** — it used `HOOP=1` and is therefore affected.
4. A `t_r` sensitivity sweep over the measured 0.1–2 s bracket, which would retire the
   clock as a threat to any conclusion if impact-phase results are flat across it.

**Still open and unresolved:** the `jam_frac` factor-0.44 mismatch against the 29 JUL
campaign, where the flow reproduces and the state variable does not.

**New landmine worth knowing:** a state variable of order **1e-25** — 25 orders below
anything physical — is enough to change this solver's mesh and move `umax` by 13%, because
`fj` is one of the four fields `adapt_wavelet` refines on. Nothing here depends on it, but
do not attribute a small difference between two campaigns to physics without checking this
first.

# PRE-REGISTRATION — constant-PRESSURE shear: the shakeout

**Written 2026-09-12, before any constant-pressure cell was launched.** This registers
the SHAKEOUT only (6 cells, N = 600). The production ladder is a separate
pre-registration that may only be written after the bars below are read.

Discharges `NEXT_STEPS_2026-09-11.md` §2 item 5. Code: `shear.shear_at_constant_pressure`,
`shear.servo_edot`, `shear.dilatancy_rate`. Grader: `grade_constp.py` →
`out/constp_shakeout.txt`.

---

## 1. Why this protocol is worth the risk

At **constant pressure** the steady-state solid fraction **is** `φ_c(p)`. No ladder, no
power-law root, no fitting window, no isostatic contact target. Every `φ_c` estimator
this project has used needs at least one of those, and at the production friction all of
them end up reading the isotropic `φ_J` the shear ladders are centred on
(`PREREG_ROOT_PINNING_2026-09-09.md`; restated in `PREREG_ZISO_2026-09-12.md` §8, where
correcting the target moved `Δφ_c` by 0.0003 and therefore could *not* have been the
cause of the E4/E4m split). **This is the first estimator here that is structurally
independent of `φ_J`.** It also returns `α`, the dilatancy coefficient, which
`derivations/dilatancy_sign_mechanism.md` currently **imports rather than derives**.

And it is the one protocol whose answer can be corrupted by its own controller, which is
why every bar below is about the controller and not about the physics.

## 2. The mechanism, and where `α` comes from

The grains are incompressible, so `V φ` is constant and `ε̇_v = d ln V/dt = −d ln φ/dt`.
With the Reynolds law `ε̇_v = α γ̇ (φ − φ_c)`,

```
    -(1/φ) dφ/dt = α γ̇ (φ − φ_c)        =>    dφ/dγ = −α φ (φ − φ_c)
```

which, linearised about `φ_c`, is an exponential in the **strain** with rate `α φ_c`. So

```
    α = −[ d ln|φ − φ_c| / dγ ] / φ_c
```

and the sign structure is a check, not an assumption: a pack denser than `φ_c` must
dilate (`φ` falls) and a looser one must contract (`φ` rises), so the fitted rate must be
**negative in `|φ − φ_c|` from either side**, and the two sides must return the same `α`.

## 3. The servo

`ε̇_v = K γ̇ (p/p* − 1)`, clamped at `|ε̇_v| ≤ 0.5 γ̇`, applied as the diagonal of
`velGrad` while the off-diagonal stays `γ̇`. It is **strain-clocked**: the gain `K` is
dimensionless and `γ = γ̇ t` remains the only clock. Sign: Yade's `velGrad` diagonal is
positive in expansion and `_stress_hat` returns `p` positive in compression, so an
over-pressured cell must expand. (Unit-tested in `test_constp.py` without Yade.)

## 4. The shakeout grid — 6 cells

`cv = 0`, `μ_p = 0.3`, lognormal, `N = 600`, seed 1, `γ_max = 3`, `shear_steps = 200000`
(so `γ̇` = 6e-05 and `I ≈ 3.5e-03`, inside the quasi-static bar), `p* = 3e-04` — **a rung
of the isotropic ladder**, so the answer can be set beside an independent measurement at
the same pressure.

- `φ_start ∈ {0.55, 0.64}` — deliberately one below and one above the expected `φ_c`.
- `K ∈ {0.25, 0.5, 1.0}`.

⚠ `N = 600` is the smoke size, not the production size (1500). Nothing here is a
measurement of `φ_c`; the shakeout decides whether the protocol is admissible.

## 5. Bars, fixed here

**S1 — the servo holds.** `|⟨p⟩_ss/p* − 1| ≤ 0.05` on **≥ 5 of 6** cells. *Predicted:
6/6.* Fewer than 5 means the gain range is wrong and the grid is re-chosen, not the bar.

**S2 — 🔑 the answer does not depend on the gain.** At fixed `φ_start`, the spread of
`φ_ss` across `K ∈ {0.25, 0.5, 1.0}` must be **≤ 0.005** (the seed scatter in `φ_J` is
0.0011–0.0048, so this is that scale). *Predicted: ≤ 0.005.* **If it fails, `φ_ss` is
tracking the controller rather than the material and the protocol is not admissible** —
no production run, and the failure is the result.

**S3 — 🔑 the critical state forgets its preparation.** At matched `K`,
`|φ_ss(0.55) − φ_ss(0.64)| ≤ 0.005`. *Predicted: ≤ 0.005.* This is the defining property
of a critical state and it is being tested rather than assumed. A failure means `γ_max`
= 3 is not enough strain (diagnosable: the `φ` series will not have flattened), which is
a different fix from S2's.

**S4 — `α` is admissible and consistent.** The fitted rate must be negative on **6/6**
cells, `|φ − φ_ss|` must move **≥ 0.3 decades** for a rate to be quoted at all, and the
two starting sides must agree within a **factor of 2** at matched `K`. *Predicted:
negative 6/6; agreement within 2× is genuinely uncertain — the dilating and contracting
branches are not required by anything to have the same coefficient, and a clean
disagreement is a physical result rather than a failure.*

**S5 — the cell is a valid measurement at all.** `I ≤ 1e-2` and both steadiness channels
(`p` halves and `φ` halves within 15 %) on ≥ 5 of 6. Reported per cell.

**S6 — the admissibility sanity check, which cannot be passed by a broken servo.**
`φ_ss` must lie **strictly between** 0.55 and 0.64. A `φ_ss` outside the bracket means
the servo drove the volume somewhere neither starting state nor the material asked for.
*Predicted: inside.*

**S7 — the external-ish comparison, reported not gated.** At `p* = 3e-04` the project's
own **isotropic** ladder gives `φ = 0.583314 ± 0.000725` (n = 9, `cv` = 0, `μ_p` = 0.3,
N = 1500). Steady shear should reach a **denser** state than isotropic compression at the
same pressure, so *predicted `φ_ss` > 0.5833*. And `E4m` on the common ladder reads
`φ_c(cv = 0)` = **0.5902**; if `φ_ss` lands within **0.015** of that, an estimator
needing no isostatic target has corroborated one that does — which is the strongest thing
this shakeout could produce. ⚠ Reported with the `N` = 600 vs 1500 caveat attached, and
it is not a gate: finite size alone could move it.

## 6. What makes this uninformative

If S2 fails, nothing else in the run means anything and the honest outcome is "the servo
is the measurement". If S1 fails while S2 passes, the gain range needs widening and the
protocol is still alive. Those are different diagnoses and the grader must say which.

## 7. OUTCOME — 2026-09-12, `out/constp_shakeout.txt`, 10 cells, 0 failures

**The protocol is admissible on both deciding gates, and the gate that failed came with
its own diagnosis.** 6 cells as registered, plus the 4 low-gain cells §5's S1 prescribes
in advance ("fewer than 5 means the gain range is wrong and the grid is re-chosen, not
the bar").

| gate | bar | measured | verdict |
|---|---|---|---|
| S1 servo holds | ≤ 5 % on ≥ 5/6 | **4 of 10**, worst 23.6 % | **FAIL** (diagnosed below) |
| S2 gain-independent | spread ≤ 0.005 | 0.0027 (all 10 cells) | PASS |
| **S2′** gain-independent, S1-passing cells only | ≤ 0.005 | **0.0023** | **PASS** |
| S3 forgets its preparation | gap ≤ 0.005 | 0.0032 | PASS |
| **S3′** same, S1-passing cells only | ≤ 0.005 | **0.0032** | **PASS** (n = 1 gain) |
| S4 `α` consistent | sign 6/6, sides ≤ 2× | sign **10/10**, ratio up to **4.41×** | NOT PASSED |
| S5 valid measurement | ≥ 5/6 | 7 of 10 | PASS |
| S6 inside the bracket | strictly inside | 10 of 10 | PASS |
| S7 comparison | reported | see below | — |

### 🔑 What the deciding gates say

`φ_ss` = **0.5935 ± 0.0012** over the seven valid cells, and it moves by **0.0023 across
a 5× range of servo gain** and **0.0032 across preparations 0.09 apart in `φ`**. Those
are the two ways this protocol could have been meaningless, and it is neither. At
constant pressure the steady solid fraction is set by the material.

### 🔴 S1 failed, and the recorded series say exactly why — twice over

The failure is **not** controller droop. Splitting `p_hat` and `z_series` over the
averaging window:

```
   K = 0.05 / 0.10 / 0.25 :  z_win_min 3.8-4.8 ,  p_win_min 4e-05 - 1e-04   STAYS JAMMED
   K = 0.50 / 1.00        :  z_win_min 0.02-0.05, p_win_min 1e-08 - 4e-08   UNJAMS
```

**Above K ≈ 0.25 the servo drives the cell through unjamming** — `p_hat` visits 1e-08
against a 3e-04 set point, six decades low, and the contact network empties (`z_nr` to
0.02). The mean pressure then reads **+20 % HIGH** rather than low, because the
excursions are one-sided: `p` is bounded below by zero and unbounded above. And below the
boundary the opposite limit bites: at K = 0.05 the servo is too slow to track and the
offset returns (+8.7 %), with the `p` halves disagreeing by 28–33 %. So a *proportional*
servo has a usable window of roughly **K ∈ [0.1, 0.25]** here, and its residual offset is
a systematic that must be priced rather than assumed away.

**Priced, on the project's own isotropic ladder:** `dφ/d ln p` = **0.00906** (0.583314 at
p = 3e-04 to 0.594219 at 1e-03), so a 20 % pressure offset is worth **0.0017 in `φ`** —
a third of the S2 bar and a sixth of the effect in S7. A production run should carry an
integral term or a tighter clamp, and must **report `z_win_min` per cell**, because a cell
can hold `⟨p⟩` on average while visiting an unjammed state.

### 🔴 And S2 as written had passed on a different question

S2 and S3 average over **every** cell, including the two gains that unjam. A cell that
visits `p` ~ 1e-08 is not a measurement of `φ_c` at 3e-04, so the pass was compatible with
the very thing the shakeout exists to rule out. **S2′/S3′ restrict to the cells that
passed S1's own pre-registered per-cell criterion** — no new bar — and only then is the
verdict earned. On the first six cells alone S2′ would have been **UNREADABLE**: one
stable gain existed. That is what the four low-gain cells bought.
⚠ **S3′ still rests on a single shared held gain (K = 0.25).** It is readable, not thick.

### 🔑 S7 — and this is the result worth carrying

`φ_ss`(p = 3e-04) = **0.5935**, against:

- the **isotropic** ladder at the same pressure, **0.583314 ± 0.000725** (N = 1500, n =
  9): steady shear reaches a state **denser by +0.0102**, i.e. 14 seed scatters, in the
  direction predicted. Constant-pressure shear densifies past isotropic compression.
- **`E4m`'s `φ_c`(CV 0) on the common ladder, 0.5902**: `|Δ|` = **0.0033**.

So an estimator that needs **no ladder, no root, no fitting window and no isostatic
contact target** lands 0.0033 from the one that needs all four. That is the first
corroboration of `φ_c` in this project from a route that is structurally independent of
the isotropic `φ_J`, and it answers the standing complaint in
`PREREG_ZISO_2026-09-12.md` §8 and `PREREG_ROOT_PINNING_2026-09-09.md` §7.
⚠ **N = 600 here against 1500 there**, one seed, one pressure, one CV. Finite size alone
could account for a move of this size. It is a corroboration, not a measurement, and the
production ladder is what settles it.

### `α`, the dilatancy coefficient — first numbers, and a split

The fitted rate is negative on **10/10** (the pack always approaches its own steady state)
and `|φ − φ_ss|` moves 0.66–4.08 decades, so the fits are admissible. `α` lands in
**1.5–7.3**, and the **dilating and contracting branches disagree by up to 4.41×** at
matched gain. §5 pre-registered that a clean disagreement here is a physical result rather
than a protocol failure — a pack that must dilate to flow and one that must contract are
not required to do it at the same rate — but with one seed, one `N` and a servo whose
offset varies cell to cell, **no `α` is quotable yet**. What is established is that the
protocol *returns* an `α` with the right sign structure, which is what
`derivations/dilatancy_sign_mechanism.md` needs to stop importing the coefficient.

### Verified this session, separately

`shear._prepare` was extracted verbatim from `shear_at_constant_volume` so both drivers
share stages A–C. `check_prepare_refactor.py` runs the **pre-refactor** module (out of
git, via a sidecar) and the current one on the same cell: **bit-identical on all five
continuous channels** (`p_ss_hat` = 0.0030137350194985426 both). So the constant-pressure
driver's preparation is the same code the whole shear campaign was measured with.

### What the production pre-registration must fix

1. `K` in [0.1, 0.25], or an **integral term** so the offset is not a systematic.
2. `z_win_min` and `p_win_min` reported per cell; a cell that unjams is refused.
3. N = 1500, ≥ 3 seeds, and the CV ladder — the whole point is `Δφ_c`(CV) from a route
   that does not inherit `φ_J`.
4. `p*` on the isotropic rungs {3e-05, 1e-04, 3e-04, 1e-03}, which also gives
   `φ_c(p)` — the shear-minus-compression offset at four pressures with no estimator on
   either side.
5. More strain, or a demonstration that `γ_max` = 3 suffices: the two K = 0.05 cells
   failed the `p`-steadiness channel at 28–33 %.

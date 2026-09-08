# Pre-registration — the containment CURVE, not its two endpoints

Written **2026-09-08, before the runs.** Not rewritten; errors go in an addendum.

## Why

`RA_LADDER_2026-09-02.md` closed the Rayleigh axis and left one line open on the Status
card:

> ⬜ Not done: `η` is sampled at **two points, not swept**, and `flux = 0` throughout, so
> the provisioning axis is not crossed with this one.

`run_ra_ladder.py` says so in its own docstring: *"Only the two containment endpoints that
the ROI actually needs — `eta` = 0 and `eta` = 0.95 — because the question is whether the
GAP between them survives at high Ra, not the shape of the curve in between."*

**The gap is the wrong deliverable for a buyer.** Nobody buys containment as an on/off
switch: they buy blanking panels, aisle doors and ceiling curtains incrementally, and what
decides the spend is **the shape of `φ(η)` between the endpoints** — whether the return is
proportional, or whether it saturates at some `η*` beyond which more containment buys
nothing. Two points cannot distinguish a straight line from a curve with a knee.

## What is already known, and what it predicts

From the existing ladder, at each Ra, a **purely proportional** leak model
`φ = φ(0)·(1 − η)` predicts `φ(0.95) = 0.05·φ(0)`:

| Ra | `φ(0)` | proportional prediction at η = 0.95 | measured `φ(0.95)` | excess |
|---|---|---|---|---|
| 4e3 | 0.2411 | 0.0121 | 0.0729 | **6.0×** |
| 1e6 | 0.1824 | 0.0091 | 0.0249 | **2.7×** |

So `φ` is **not** proportional to the leak, there is an **excess above proportional**, and
that excess **shrinks as Ra rises** — which is exactly what the 2026-08-31 Lewis-number
retraction says it should, because the excess is a *diffusive floor* and raising Ra
strengthens advection against diffusion.

## The model under test, fixed now

    phi(eta)  =  A * (1 - eta)  +  F        with A, F >= 0

`A` is the advective leak term, proportional to the opening; `F` is the floor. Both are
fitted per Ra rung from the sweep. This is a two-parameter form fitted to nine points, so
it can fail.

## Gates

### GD1 — the positive control, and it comes first

The η = 0 and η = 0.95 points in this sweep must reproduce the existing ladder's values at
the same (Ra, N). The solver is deterministic, so this is close to an identity check.

**Bar: `|Δφ| ≤ 2 × PHI_TOL` = 1.6e-3**, the solver's own convergence tolerance. If it
fails, the two scripts are not running the same cell and nothing below means anything.

### GD2 — is `φ(η)` proportional to the opening?

Fit `φ = A(1 − η) + F` on each Ra rung.

**Bar: the one-parameter proportional model (`F ≡ 0`) is REJECTED on every rung**, in the
sense that its residual RMS exceeds the two-parameter form's by more than 3×.

**Registered prediction: it is rejected on every rung**, because the excess above
proportional is already 2.7–6.0× at the one interior-free point available.

### GD3 — does the floor fall with Ra, and does the leak term not?

**Registered predictions**, both directional and both falsifiable:

1. `F` **falls monotonically** with Ra across the four rungs.
2. `A` **does not fall by more than 25 %** over the same range — i.e. the Ra dependence of
   `φ(0.95)` is carried by the floor, not by the leak.

If (2) fails, the "the floor is diffusive" reading is not the whole story and the
2026-08-31 retraction's replacement account needs revisiting.

### GD4 — the tracer channel, which is the fit-free version of the same question

`phi_tracer` is a passive scalar with no buoyancy coupling. If the floor is diffusive, the
tracer's own floor should be far smaller and should fall faster with Ra.

**Reported, not barred**: `F_tracer` per rung beside `F`. No bar because the tracer's
diffusivity is a solver parameter rather than a physical one here, and a bar on it would be
a bar on a knob.

### GD5 — the deliverable: where is the knee?

Define `η*` as the containment efficiency at which `φ` has fallen to within 10 % of its
`η → 1` limit, i.e. `A(1 − η*) + F = 1.1 F`, so `η* = 1 − 0.1 F/A`.

**No bar. This is a number to report, not a hypothesis** — its value is the answer a buyer
needs, and inventing a bar for it would be inventing a target.

## What this may NOT claim

- **Nothing about a real hall.** `Ra` here reaches 1e6 against 1e9–1e12 in a hall, and
  above `Ra` ~ 1e8 a 2-D steady laminar solver is the wrong model. Those limits are
  `run_ra_ladder.py`'s and are unchanged.
- **Nothing about provisioning.** `flux = 0` throughout, as in the Ra ladder. Crossing η
  with flux remains open and is explicitly not attempted here.
- **`η` is a model parameter, not a measured containment efficiency.** No vendor's
  containment product has been mapped onto it.

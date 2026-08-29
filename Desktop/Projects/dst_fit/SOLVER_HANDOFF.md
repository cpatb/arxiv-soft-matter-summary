# What the Lin deposit can do for the Basilisk WC solver

Written 2026-08-28, answering "can yesterday's data better constrain the STF solver?"
Companion: `DATA_INVENTORY.md` §2 (what the deposit is), and the vault card
`Status\Drop Impact — Simulation.md` (what the solver needs).

**Short answer: yes, in two places, and neither of them is `phi_m`.** The solver's own
card says the two things below are exactly what it is missing. Both are cheap.

Reproduce with `python scripts/probe_lin_fig2b_recovery.py` from the package root
(needs the gitignored DataShare deposit).

---

## 1. The strain clock `gamma_0` — the solver's biggest unpinned knob

The solver runs Richards et al. Eq. 4 under `WC_CLOCK=1`:

```
dfhat/dt = -(gammadot/gamma_0 + 1/t_r) (fhat - fhat_ss)
```

with **`gamma_0` ~ 0.2 borrowed from Richards' corn starch** (and Han et al.'s 0.197 on
another system), and `t_r = 0.24 +/- 0.05 s` bracketed from four flow-cessation studies.
The sim card states plainly that these "cannot be pinned without step-down rheometry" and
that `tau_relax`/`tau_up` remain PLACEHOLDER.

**Lin `Fig2b` is that measurement, on silica.** Instantaneous relative viscosity against
accumulated strain after flow reversal, 2 um silica, four rates 0.5 to 5 1/s, ~10k to 25k
points each. Shear reversal destroys the contact network and rebuilds it over accumulated
strain: that rebuild *is* `gamma_0`.

### What the data already says

Landmarking the 63% point of the min-to-steady recovery, minimum restricted to
`gamma < 2` (the global argmin lands on a one-point dip at the very end of the 2.0 and
5.0 records — that bit me on the first pass):

| gammadot [1/s] | gamma at min | eta_min | eta_ss | gamma_63 | t_63 [s] |
|---|---|---|---|---|---|
| 0.5 | 0.051 | 4.85 | 8.09 | 1.595 | 3.191 |
| 1.0 | 0.348 | 4.56 | 9.81 | 1.320 | 1.320 |
| 2.0 | 0.264 | 4.25 | 12.46 | 1.084 | 0.542 |
| 5.0 | 0.244 | 4.56 | 14.70 | 0.981 | 0.196 |

**Across a 10x span in rate the recovery collapses on strain, not on time:**

- strain: `gamma_63 = 1.245 +/- 0.273`, spread **1.63x**
- time: `t_63 = 1.31 +/- 1.34 s`, spread **16.3x**

A pure time clock would have given a 10x strain spread and a flat time spread. It is the
other way round. **The strain-clock term is the one that governs in this window** — which
is an independent confirmation of the law the solver already runs, on the solver's own
material rather than on Richards' corn starch.

### Two things not to overstate

1. **`gamma_63 ~ 1.2` is not `gamma_0 = 0.2`.** It is the e-folding strain of the
   *viscosity*, and `eta` is strongly nonlinear in `fhat` near jamming. Converting one to
   the other means fitting through the Wyart-Cates relation, which is what this package
   does and what the landmark does not. The factor ~5 gap is the reason to do the fit
   properly, not a result yet.
2. **The residual 1.63x strain spread has the wrong sign for the `1/t_r` term.** Richards'
   form gives an apparent strain `gamma_0 / (1 + gamma_0/(gammadot * t_r))`, which
   *increases* toward `gamma_0` with rate. The measured `gamma_63` *decreases*
   monotonically with rate (1.60, 1.32, 1.08, 0.98). Either the nonlinearity in (1)
   accounts for it — higher rate means higher `fhat_ss` and a compressed apparent strain —
   or the additive-time form is wrong here. Worth knowing which before quoting a `t_r`.

Also: the `gammadot = 0.5` record's minimum sits at `gamma = 0.051` against 0.24 to 0.35
for the other three, and it is shallow. Its `gamma_63 = 1.60` is the least trustworthy of
the four and it is also the one driving the spread.

**Why this is load-bearing rather than tidy-up.** The sim card's L9 row reads peak
`jam_frac` = 0.912 / 0.997 / 0.671 at `gamma_0` = 0.05 / 0.2 / 0.8. The observable is not
flat in `gamma_0`, and a factor of a few in `gamma_0` moves it. A measured `gamma_0` on
silica replaces a corn-starch literature value in the one term the solver claims as its
differentiator.

---

## 2. `sigma*` — and a 179x discrepancy that wants one gate

Gate G10 already fitted `sigma*` from Lin's 2 um silica flow curve and matched it against
the shear-reversal contact/hydro crossover to **0.73x**. That is the validated half of
Wyart-Cates. The solver's `sigma*`, by contrast, comes from third-party rheometry the repo
itself flags as "provisional (not paper-grade)", on a fluid it suspects is not water.

Comparing them has to be done in reduced stress `Sigma = sigma/eta_s = eta_r * gammadot`,
because Lin never states a solvent viscosity — the same trick that made G10 work:

| quantity | value |
|---|---|
| solver fitted silica, `sigma*/eta_s` = 44.2327 / 0.004736 | **9340 1/s** |
| Lin 2 um silica, fitted `Sigma*` (G10) | 8.70 1/s |
| scaled 2.0 um -> 0.78 um via Guy's `d^-1.9` (x5.98) | 52.1 1/s |
| **ratio** | **179x** |

To agree, the solver's `sigma*` would have to be **0.25 Pa, not 44.2 Pa**.

### RETRACTED, 2026-08-28 evening. The 179x was an artefact of the comparison.

Gate `scripts/gate_g12.py` now does this properly and the conclusion above does not
survive. Three things were wrong with it.

**1. Reduced stress does not cancel `eta_s` across two datasets.** It cancels within one,
which is exactly what makes it correct in G10 and wrong here: G10 compares two quantities
measured on the *same sample*, so the unknown appears on both sides. Comparing the solver
against Lin in reduced stress silently assumes the two solvents have the same viscosity.

**2. Lin's solvent viscosities are not unstated — the deposit contains them twice.**
`Fig4a` is absolute viscosity in Pa.s and `Fig4b` is relative viscosity for the same
sample, so their ratio at matched stress is `eta_s`; and `FigS3` measures the bare PMMA
solvent directly. The two independent routes to PDV-2331 agree to **0.29%** (1.770 against
1.765 Pa.s), which is what licenses the same method on the silica:

| system | solvent (named in the Fig4 README) | `eta_s` at 20 C |
|---|---|---|
| 2.0 um silica | glycerol + water | **36.3 mPa.s** |
| 1.4 um PMMA | PDV-2331 | **1.770 Pa.s** |

With `eta_s` measured, **179x becomes 24.6x**, and the factor that vanishes is exactly
`eta_s(Lin)/eta_s(solver) = 36.3/4.736 = 7.66x`.

**3. `eta_s` was the wrong prime suspect anyway.** A 2.1x error in `eta_s` moves the fitted
reduced `sigma*` by about 5%, because `eta_s` enters the fit and both axes, not just the
divisor. Closing 179x that way would need `eta_s ~ 0.85 Pa.s`. There is no such solvent.
Demote it.

### What the data actually says, and it is not a discrepancy

The deposit holds **two** silica anchors, and they disagree with each other by more than
the quantity being tested:

| anchor | `sigma*` at 2 um | scaled to 800 nm (`d^-1.9`) |
|---|---|---|
| `Fig2a`, dense (`eta_r` 6.4 to 19.1) | 0.316 Pa | 1.80 Pa — solver is 24.6x **above** |
| `Fig4b`, dilute (`eta_r` 1.5 to 1.85) | 57.0 Pa | 325 Pa — solver is 7.4x **below** |

Wyart-Cates says `sigma*` is a repulsive force scale and therefore `phi`-independent, and
on the **PMMA** pair the same two routes agree to **1.9x** — so the method is sound and the
silica spread is not an analysis artefact. The likely reason they differ: PMMA is
sterically stabilised in a non-polar oil, while aqueous silica is electrostatically
stabilised and its counterion concentration rises with `phi`, so its repulsion is screened
at high `phi`. **`sigma*` for aqueous silica is not a transferable material constant** —
which is what this lab's own salt series says from the other direction.

The solver's 44.2 Pa sits **inside** that 180x span, at every size exponent from 1.8 to
2.2. So this data cannot impeach the solver's `sigma*`. That is a very different statement
from "it would have to be 0.25 Pa", and the honest verdict is that the literature is too
widely spread to test it.

**What would test it:** a `sigma*` measured on the solver's own silica at its own ionic
strength. Nothing transferred from another preparation will do, because the 180x above is
the size of the preparation dependence.

Two candidates from the original list survive as worth checking, neither of them `eta_s`:
the **phi label convention** (`SILICA_CONC_IS_VOLPCT` is still a live flagged knob in
`calibrate_from_rheometry.py`), and — the one this list missed — **whether the solver's own
stress sweep brackets its transition at all**. G5 characterised exactly that failure mode:
when the sweep sits on one side of the onset, `sigma*` does not run to a bound, it migrates
into the measured window and returns a plausible error bar while being wrong by orders of
magnitude. `transition_coverage()` already exists and answers it in one line. Run it on the
solver's rheometry before trusting `44.2 Pa` for anything.

---

## 3. What the data explicitly CANNOT do: `phi_m`

The solver's most sensitive input is the one this package just spent a session proving is
not measurable from a flow curve.

- Lin is **a single volume fraction per system**, so it cannot test `phi_m` consistency at
  all — that is already recorded in `DATA_INVENTORY.md`.
- More importantly, the project's central result applies directly: `phi` and `phi_m` enter
  only through `phi/phi_J`, so near jamming a flow curve pins the ratio and neither number.
- And the sim card's own lamella result says **"the driver is phi vs phi_M, not sigma*"**,
  with the fit putting operating `phi = 0.51` at the jamming edge and noting
  "since phi ~= phi_M within fit uncertainty the model sits in its most sensitive regime,
  so exact numbers are soft".

**"Soft" is the wrong word for it, and that is the real transfer between the two projects.**
It is not that `phi_m` is imprecise and better data would sharpen it; it is that a steady
flow curve does not contain it, so no amount of better steady rheometry will. The solver's
`PHIM = 0.5176` is not a measurement with an error bar, it is a number the fit had to
choose. Two consequences:

1. The solver should carry `phi_m` as a **swept** parameter with the ordering quoted, not
   as a calibrated one with a value quoted — the same discipline the sim card already
   applies to `jam_frac` ("quote the ordering, not the value").
2. `params_fitted.yaml` sets `PHI0 = 0.7`. This package measured `phi_0 = 0.6076(22)` and
   `0.6123(10)` on two PMMA sizes and showed that **`phi_0` sets the level of `phi_m` and
   leaves the spread bit-for-bit identical**. So `PHI0 = 0.7` is directly biasing
   `PHIM = 0.5176` by a knowable amount, and it is not a fitted quantity.

---

## Suggested order, if this gets picked up

1. **G11** — fit `gamma_0` from Lin `Fig2b` *through* the Wyart-Cates relation rather than
   by landmarking `eta`, using all four rates jointly. Deliver `gamma_0` for 2 um silica
   with an interval, and settle whether the wrong-sign rate trend is the `eta`-vs-`fhat`
   nonlinearity or a defect in the additive-time form. Cross-check against `Fig2d` (1.4 um
   PMMA) for a second material, and against `Fig3` (simulation, contact/hydro evaluated
   directly) as the control on the method itself.
2. **G12** — the reduced-stress `sigma*` comparison in §2, as a gate, so the 179x either
   resolves or becomes a stated finding about the third-party rheometry.
3. Hand `gamma_0` and the `phi_m`/`phi_0` caveat to the solver as a parameter note.

Both gates are new cross-checks on parameters this package does not currently fit, which
is the same thing that made G10 the strongest result in the project.

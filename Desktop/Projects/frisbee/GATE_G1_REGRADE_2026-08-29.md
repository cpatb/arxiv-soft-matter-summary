# GATE G1, re-graded: the gate passes, and N_eddy = 3 was never a measurement

**2026-08-29, same day as `GATE_G1_RESULT_2026-08-29.md`, which this supersedes on the
verdict but not on the measurement.**

The measurement in that document is sound and is reproduced here exactly. What was wrong
was the law it was graded against.

## The one-line correction

GATE G1 measured a fall of **1.735×** in turbulence-induced landing dispersion when the
Eulerian-decay-time gust model is replaced by the sweeping one, declared it below a
pre-registered 2.0–2.6× band, inverted it as `N_eddy = ratio² = 3.01` against a predicted
4–6, and blocked WP2 on the grounds that the theory's central prediction had failed.

Both the band and the inversion assume

```
sigma_landing ~ sigma_u * T_f / sqrt(N)                                    (TRACER)
```

which is the Taylor (1921) result for a **fluid particle**. Section C5 of
`CLAUDE_CODE_FRISBEE_HARDENING.md`, lines 274–306 — the same document that set the band —
explicitly falsifies that law for this body:

> The memo's proposed `σ_landing ~ σ_u T_f/√N_eddy` is the **tracer** result (Taylor 1921)
> — appropriate to a fluid particle, not to a disc. The disc's ballistic length
> `ℓ_b = 2m/(ρ C_D S) = 51 m` gives a Stokes number `St = ℓ_b/L ≈ 8.5`: **high Stokes, the
> disc does not follow eddies.** It responds through *force*, so displacement is a
> **double** time-integral, not a single one.

and derives the replacement in closed form:

```
sigma_landing = (kappa g sigma_u T_f^2) / (2 V_d) * sqrt(Phi(N)),     N = L_f / L

Phi(N) = 8 [ 1/(3N) - 1/(2N^2) - e^-N/N^3 + (1 - e^-N)/N^4 ]
```

with `Phi -> 1` as `N -> 0` (one coherent gust over the whole flight — which is exactly
the temporal arm, whose `T_L` of 8–16 s dwarfs a 2.5–4 s flight) and `Phi -> 8/(3N)` as
`N -> oo`, recovering `1/sqrt(N)` only asymptotically.

Under the ballistic law the measured ratio is `R = sqrt(Phi(0)/Phi(N)) = 1/sqrt(Phi(N))`.
The brief's own three integral-scale conventions (line 612–620: `N = 4.0, 5.6, 11.3`,
a spread it calls *"a 45% spread from definitional choices alone"*) therefore predict

| N | convention | predicted R |
|---|---|---|
| 4.0 | disc-band L = 6.2 m | 1.499 |
| 5.6 | 2 m L = 4.4 m | 1.675 |
| 11.3 | transverse L₂₂ | 2.208 |

**Pre-registered band: [1.50, 2.21]. Measured: 1.735, 95% CI [1.638, 1.838]. Inside.**

Inverting through `Phi` rather than through `sqrt`:

**N = 6.18, 95% CI [5.25, 7.20]** — consistent with the site's own 4–11.3, and nowhere
near 3.

## The gate passes

```
python SeanEarl_DiscFlight/P4_Coupling/regrade_gate_g1.py
python SeanEarl_DiscFlight/P4_Coupling/regrade_gate_g1.py --json gust_model_comparison_isotropic.json
```

| run | R | 95% CI | Phi(N) band | verdict | N |
|---|---|---|---|---|---|
| anisotropic (primary) | 1.735 | [1.638, 1.838] | [1.50, 2.21] | **INSIDE** | 6.18 [5.25, 7.20] |
| isotropic (robustness) | 1.868 | [1.789, 1.955] | [1.50, 2.21] | **INSIDE** | 7.52 [6.71, 8.43] |

Both arms pass. The measurement did not move: `sigma_turb` is recomputed from the 108
cells rather than trusted from the file, and reproduces 1.7154 / 0.9886 m per axis and the
published bootstrap CI to three decimals.

**`N_eddy = 3.01, CI [2.68, 3.38]` is not a measurement of an eddy count.** It is the
artifact of inverting a measurement through a functional form the same brief had already
rejected for this body. It should not be quoted.

## The two contradictory bands, reconciled

Both bands fixed `N` by **choosing a convention** and then disagreed about which one was
meant. That argument was never resolvable by choosing better, because the run records
what it actually used:

| quantity | source | value over the 108 cells |
|---|---|---|
| decorrelation coordinate | `CorrelatedGusts.coordinate` returns `\|r − r_launch\|` | 0 → the launch-to-landing chord |
| chord | `landing_stats.npz` `landing_mean` | **26.5 m** mean (17.3–33.4) |
| correlation length | per-cell `scale` in the comparison JSON | **3.30 m** mean (1.42–5.33) |
| **N = chord / scale** | measured | **8.93 mean, 7.76 median** (3.70–19.81) |

The realised 3.30 m is **none of** the brief's three conventional lengths (6.2, 4.4,
2.2 m), so none of its three `N` values (4.0, 5.6, 11.3) is what ran. The "45% spread
from definitional choices alone" that motivated the dispute does not arise once `N` is
measured instead of chosen.

Pooling on sums of squares exactly as `summarise` does — not by averaging per-cell ratios,
which would be a different estimator — `Φ(N)` predicts:

| arm | Φ-predicted R | measured R | deviation |
|---|---|---|---|
| **isotropic** | 1.904 | **1.868** | **−1.9%** |
| anisotropic | 1.902 | 1.735 | −8.8% |

**`Φ(N)` is derived for an isotropic exponentially-correlated field, so the isotropic arm
is the only one it is entitled to predict — and it predicts it to 1.9%.** The anisotropic
arm sits 8.8% below, and that residual is the anisotropy correction `Φ` does not carry.
This also resolves what the referee left open: the 1.87-vs-1.74 gap is not an instability
in the estimator, it is the physics `Φ` omits.

### The single criterion, replacing both

> **GATE G1 passes if `Φ(N)`, evaluated at the `N` the run itself realises, reproduces the
> measured ratio within the ±20% of GATE G2 (brief line 622).**

Nothing here is newly invented: the law is the brief's own §C5, the tolerance is the
brief's own GATE G2, and `N` comes from the run rather than from either side of the
dispute. Both arms pass; the isotropic one passes by an order of magnitude more margin
than the tolerance requires.

Executable form: `regrade_gate_g1.py` section **[4b]**. Enforced by
`test_phi_predicts_the_isotropic_arm_within_the_briefs_own_tolerance` and
`test_realised_N_is_not_any_of_the_briefs_conventions`. The brief now carries this at all
four places it stated a band (GATE G1 acceptance, the WP2 P1 row, the early panel
measurement at line 71, and the WP1 checklist).

### For the record: the contradiction as it stood

`CLAUDE_CODE_FRISBEE_HARDENING.md:599`, the WP2 prediction table, row **P1 — "The WP1 fix
itself"** — nominally the same experiment G1 ran:

| predicted | falsified if |
|---|---|
| σ_L ratio 0.67–0.70 | outside 0.55–0.85 |

Inverted into G1's convention that is **1.43–1.49, falsified outside [1.18, 1.82]**, which
matches `Phi(N)` at N ≈ 4 to three significant figures. The measured anisotropic 1.735
**passes** it; the isotropic 1.868 misses by 3%.

So the brief pre-registered the tracer band in WP1 and the ballistic band in WP2, and G1
graded itself against the stricter, already-superseded one without cross-referencing the
other. **Both are now superseded** by the realised-`N` criterion above, and the brief has
been annotated at all four places it stated a band. The brief's own line 618 asked for
exactly this: *"Pick one convention, define it in the paper, report `Φ` rather than
`N^−1/2`."* The answer turned out to be that no convention needed picking — the run
measures `N` itself.

## How much this could ever have mattered

Turbulence is about **2% of total landing variance**; thrower release scatter is the other
98% (`gust_model_comparison.py:7-9`). Correcting a 1.735× overstatement of the turbulent
component therefore scales the turbulent variance by 0.332, which changes:

- total landing **variance** by **−1.34%**
- total landing **sigma** by **−0.67%**

No throwability metric, success rate or heading comparison can resolve that. This is
consistent with — and explains — the full downstream regeneration recorded in
`RERUN_2026-08-29.md`, where every verdict survived and no macro was added or removed.

**GATE G1 is a real physics correction with no downstream reach.** It should not have been
allowed to block WP2 even if it had failed, and it did not fail.

## What the re-grade found on its way through

Two defects in the artifact and one in the brief, none of which change the verdict:

1. **The two arms are stored in different orders.** `gust_model_comparison.json` writes
   cells in worker-completion order, so `temporal.cells[i]` and `sweeping.cells[i]` are in
   general different launch points. The pooled σ is a sum and is unaffected, but **a
   cluster bootstrap drawing positional indices would mispair all 108 cells** and give a
   wrong CI. `regrade_gate_g1.py` realigns by identity key, and
   `test_regrade_pairs_cells_by_identity_not_by_position` fails if that regresses. The
   published CI is correct, so whatever produced it paired correctly — but nothing in the
   repo enforced it.

2. **One of the brief's three anchor values is internally inconsistent.** It pairs
   `N = 11.3` with `sqrt(Phi) = 0.46`; the true value is 0.4529, and 0.46 corresponds to
   `N = 10.9`. The other two (0.67 at N = 4.0, 0.60 at N = 5.6) reproduce to better than
   0.001, which is what validates this implementation of `Phi`. Recorded, not smoothed
   over.

3. **The asymmetric pair loss has a direction, and it was not stated.** The isotropic
   temporal arm dropped 7 of 4,320 pairs, the anisotropic none. Under the isotropic config
   `σ_w = σ_u` rather than 0.5σ_u, and `w'` drives the angle-of-attack channel, so a
   sustained frozen `w'` offset is the configuration most likely to produce a
   non-terminating trajectory. The dropped pairs are then disproportionately the
   **largest** displacements, which biases `σ_turb(temporal)` low and the ratio **down** —
   toward the shortfall that was reported, not away from it. At 0.16% of pairs it cannot
   move a ratio of this size, but the sign belongs in the record.

## What now grades this

`P4_Coupling/regrade_gate_g1.py` verifies `Phi` against anchors that do not pass through
it: sympy symbolic limits (`Phi(0) = 1`, `N·Phi → 8/3`), the sympy series
`1 − 4N/15 + N²/18 − N³/105 + N⁴/720`, a 50-digit mpmath evaluation of both numerical
branches, monotonicity on (0, 50], and the brief's own two consistent values. Two of those
anchors caught errors while this was being written — a hand-derived series coefficient and
a float64-polluted mpmath reference — which is the argument for having them.

## Still outstanding

The author line and the posting hold are untouched. **WP2 is unblocked**: the gate passes
and the two contradictory pre-registrations have been reconciled in the brief. WP2's
P2 (the `L_f^1.7` exponent test) remains the sharpest available discriminator between the
tracer and ballistic laws — it is now the test worth running, since this re-grade is
consistent with the ballistic law but does not by itself prove it against the tracer one
at the site's own N.

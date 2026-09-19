# Item A — Does φ_J(CV, μ_p) separate?

**VERDICT: REFUTED for the multiplicative form; REFUTED (but by a factor 2.3 less) for
the additive form. The supportable statement is "additive to ~6 % of the shift".**

**What this anchor can and cannot decide.** It can RANK the two candidate forms and put a
number on the residual interaction, because their predictions are 36 se apart where the
data is strongest. It CANNOT identify a mechanism — it reads only φ_J, and two different
microscopic stories give the same two-way table. It CANNOT speak for μ_p > 0.3, whose
ladders fail the project's own health gate.

Every number below is printed by `geom_separability_2026-09-19.py` →
`geom_separability_2026-09-19.txt`. Run it; nothing here is transcribed from memory.

---

## 1. The statistic, derived

Write the seed-paired CV increment `D(c, m) = φ_J(c, m) − φ_J(0, m)`. SymPy (§1 of the
script) gives

| form | D(c, m) | prediction for `D(c,m1) − D(c,m0)` |
|---|---|---|
| additive `F(c) + G(m)` | `F(c) − F(0)` | `0` |
| multiplicative `F(c)·G(m)` | `(F(c) − F(0))·G(m)` | `(r_m − 1)·D(c, m0)` |

with `r_m = φ_J(0, m1)/φ_J(0, m0)`. Both are parameter-free once the monodisperse cells
are measured, and their **separation** is `|r_m − 1|·D(c, m0)`, which is printed beside
the residual so the test can declare itself uninformative.

**The transform-free version.** If φ_J separates additively after *any* smooth monotone
`h`, then to first order

    D(c, m1)/D(c, m0) = h′(φ_J(0,m0)) / h′(φ_J(0,m1))   for EVERY c.

That right-hand side is a single number fixed by `h`. So **if the measured ratio depends
on `c`, no monotone transform whatsoever separates the data.** SymPy evaluates the ratio
for four candidate transforms:

| h | required ratio |
|---|---|
| identity (additive) | 1 |
| log (multiplicative) | `φ₁/φ₀` |
| `1/φ` (specific volume) | `(φ₁/φ₀)²` |
| `(1−φ)/φ` (void ratio) | `(φ₁/φ₀)²` |

Specific volume and void ratio give the *same* requirement, so a free-volume reading of
the separation is one hypothesis, not two.

---

## 2. Which friction arms are admissible — a correction to fact 4

The project's own ladder-health gate (`sweep_descriptive.ladder_health`: no negative
ladder slope, median r² ≥ 0.85) is implemented at `μ_p = MU_REF = 0.3` only. Evaluated at
every arm on disk (§2):

| family | μ_p | n cells | neg slopes | median r² | verdict |
|---|---|---|---|---|---|
| lognormal | 0.0 | 18 | 0 | 0.9115 | HEALTHY |
| lognormal | 0.3 | 18 | 0 | 0.9468 | HEALTHY |
| **lognormal** | **1.0** | **9** | **1** | **0.3856** | **FAILS** |
| **lognormal** | **3.0** | **3** | **1** | **0.1421** | **FAILS** |
| uniform | 0.0 | 18 | 0 | 0.9282 | HEALTHY |
| uniform | 0.3 | 18 | 0 | 0.9207 | HEALTHY |
| bidisperse | 0.3 | 18 | 0 | 0.9210 | HEALTHY |

**Consequence for the brief's fact 4.** The three-level friction-independence statement
(+0.05803 / +0.06207 / +0.05621 at μ_p = 0 / 0.3 / 1.0) rests at its third level on
ladders the project itself rejects, and the convergence claim ("μ_p = 1.0 is converged,
t = 2.78 against μ_p = 3.0") compares two arms that *both* fail. This is not new to the
repo — `transfer_law.py:79` records "0.85 admits them and refuses the μ_p = 1.0 arm,
whose median is 0.386", and `phase11_tail.py:35` records "the μ_p = 1.0 arm on disk is
0/6 monotone" — but it is not carried by the fact as stated in the brief. Nothing in
`paper/main.tex` uses μ_p = 1.0, so the manuscript is unaffected.

Everything in §3–4 therefore runs on μ_p = 0 and 0.3 only. The μ_p = 1.0 arm is reported
in §5, labelled, never pooled.

---

## 3. The test

Seed-paired, seeds {1, 2, 3}, t_crit at 2 dof = 4.303. `r_m = 0.90854` (lognormal),
`0.90949` (uniform).

**Lognormal:**

| CV | D(μ=0) | D(μ=0.3) | resid ADD | t_ADD | resid MULT | t_MULT | separation | sep/se |
|---|---|---|---|---|---|---|---|---|
| 0.1 | +0.00512 | +0.00358 | −0.00154 | −2.24 | −0.00108 | −1.61 | 0.00047 | 0.68 |
| 0.2 | +0.01302 | +0.01352 | +0.00050 | +1.35 | **+0.00169** | **+4.68** | 0.00119 | 3.24 |
| 0.3 | +0.02595 | +0.02945 | +0.00350 | +3.94 | **+0.00588** | **+6.38** | 0.00237 | 2.67 |
| 0.4 | +0.04025 | +0.04517 | +0.00491 | +3.44 | **+0.00860** | **+5.76** | 0.00368 | 2.58 |
| 0.5 | +0.05803 | +0.06207 | **+0.00404** | **+27.70** | **+0.00935** | **+65.43** | 0.00531 | 36.38 |

**Uniform:** no residual of either form exceeds t_crit at any CV (worst additive
t = −1.94 at CV = 0.4; worst multiplicative t = +2.46 at CV = 0.3), and the candidate
separation is 0.59–3.54 se. That arm is **not informative**: both candidates pass.

**Reading.** The lognormal family at CV = 0.5 is where the design has power, and there
**both forms are refuted**. Additive is wrong by +0.00404 = 6.5 % of the shift;
multiplicative by +0.00935 = 15.1 %, and with the opposite sign to what it needs. So the
answer to "additive or multiplicative" is *additive, by a factor 2.3, and neither is
exact.*

The brief's suggestion that the additive form already looks better at CV 0 → 0.5 is
**verified**, and the size of its failure is measured rather than assumed.

---

## 4. No monotone transform separates it — NOT RESOLVED

The transform-free test of §1 needs `R(c) = D(c, 0.3)/D(c, 0)` to be constant.
Lognormal:

| CV | R | se | t vs 1 (additive) | t vs 0.9085 (multiplicative) | implied k in h = φ^k |
|---|---|---|---|---|---|
| 0.1 | 0.6985 | 0.1223 | −2.47 | −1.72 | 4.741 |
| 0.2 | 1.0381 | 0.0285 | +1.33 | +4.54 | 0.611 |
| 0.3 | 1.1351 | 0.0334 | +4.04 | +6.77 | −0.321 |
| 0.4 | 1.1221 | 0.0334 | +3.65 | +6.39 | −0.201 |
| 0.5 | 1.0696 | 0.0055 | +12.63 | +29.22 | +0.298 |

`R` looks strongly CV-dependent (0.70 → 1.14) and the implied exponent `k` even changes
sign — but the largest pairwise contrast, CV 0.1 vs CV 0.3, is **t = −3.44 against a
t_crit of 4.303**. So:

> **The claim "no monotone transform can separate φ_J(CV, μ)" is UNDETERMINED at n = 3.**

The CV = 0.1 rung, where `R` is most deviant, is also where `D` is smallest (0.00512) and
the ratio is therefore noisiest. This is the honest state and it is a design statement:
resolving it needs more seeds at low CV, not a different analysis. Uniform is worse
(largest |t| = 1.57).

---

## 5. The μ_p = 1.0 arm, and an external anchor for the friction ladder

Reported, gated out, never pooled. On `φ_J_lin` the interaction against μ_p = 0.3 is
+0.00008 ± 0.00139 (t = +0.06) at CV = 0.2 and −0.00586 ± 0.00077 (t = −7.65) at CV = 0.5.

On the **fit-free** channel φ at the matched lowest rung p̂ = 3e−5 (lognormal):

| μ_p | φ(CV=0) | φ(CV=0.5) | D | se | n |
|---|---|---|---|---|---|
| 0.0 | 0.63182 | 0.69428 | +0.06338 | — | **1** |
| 0.3 | 0.57535 | 0.63667 | +0.06131 | 0.00240 | 3 |
| 1.0 | 0.56148 | 0.61486 | +0.05339 | 0.00373 | 3 |
| 3.0 | 0.56054 | — | — | — | 0 |

(The μ_p = 0 row has only one seed at that rung and therefore no error bar — printed, not
hidden.)

**External anchor, API-resolved:** Silbert, *Soft Matter* **6**, 2918 (2010),
doi:10.1039/c001973a, Table 1, monodisperse 3D frictional spheres: φ_c = 0.639 / 0.595 /
0.556 / 0.544 at μ = 0 / 0.2 / 1.0 / 10, saturating near 0.55 for μ ≥ 1.

| μ | this work, φ at p̂ = 3e−5 | Silbert 2010 | difference |
|---|---|---|---|
| 0.0 | 0.63182 | 0.639 | −0.0072 |
| 1.0 | 0.56148 | 0.556 | +0.0055 |

The monodisperse friction ladder agrees with an independently published one to a few
times 1e−3, and both say μ_p = 0.3 is *not* saturated while μ_p = 1.0 essentially is.
That is what makes item A non-trivial: **the CV effect is stable across a friction range
over which φ_J itself moves by 0.06 — the same size as the whole CV effect.**

---

## 6. The mechanism does NOT come for free from depopulation

The rattler decomposition (`geom_rattler_decomposition_2026-09-19.md`) writes

    φ_J(CV, μ) = φ_B(μ) / (1 − ψ(CV, μ)).

**If ψ depended on CV alone, this would be EXACTLY multiplicative** — and §3 refutes
multiplicative at 65 se. So ψ must be friction-dependent. The rattler *number* fraction q
is on disk at both arms and confirms it (lognormal, §6):

| rung | CV | q(μ=0) | q(μ=0.3) | ratio | ΔQ(μ=0) | ΔQ(μ=0.3) | ratio |
|---|---|---|---|---|---|---|---|
| 1e−3 | 0.5 | 0.2771 | 0.4476 | 1.615 | +0.2700 | +0.4082 | **1.512** |
| 1e−3 | 0.3 | 0.0778 | 0.2211 | 2.843 | +0.0707 | +0.1818 | **2.572** |
| 3e−5 | 0.5 | 0.3473 | 0.5218 | 1.502 | +0.3320 | +0.4393 | **1.323** |

**The rattler count's CV response is 1.3–2.6× larger at μ_p = 0.3 than at μ_p = 0, while
φ_J's is 1.07× larger.** So whatever makes φ_J nearly separable, it is *not* that
depopulation is friction-blind. Under the backbone-invariance hypothesis the compensation
must sit in the conditional rattler VOLUME, i.e. in the selectivity S. At the lowest rung
the required S is close to friction-independent (lognormal, CV = 0.5: 0.780 at μ = 0
against 0.806 at μ = 0.3); at the top rung it is not (0.750 against 0.848). Printed in
§6, not concluded from.

---

## 7. What physics would produce additive separation

Nothing in this file derives one, and it is worth saying so plainly rather than offering a
story. What the file does establish is the shape of the constraint any such derivation
must satisfy:

1. It must be additive in φ itself, not in log φ, not in the void ratio `(1−φ)/φ`, and
   not in the specific volume — the last two are refuted hardest (they need R = 0.825
   where R is measured at 1.07 ± 0.006 at CV = 0.5, i.e. 44 se away).
2. It must leave a residual of about +6 % of the shift, of consistent sign at CV ≥ 0.2.
3. It must be compatible with a rattler-count response that is *not* separable.

Point 1 is the informative one. A free-volume account — "friction and polydispersity each
add an independent amount of excess void" — predicts additivity in `(1−φ)/φ`, and that is
the most strongly refuted of the four transforms tested. **The separation is additive in
the DENSITY, not in the void.**

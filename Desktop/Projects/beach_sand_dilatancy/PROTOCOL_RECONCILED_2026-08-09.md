# PROTOCOL.md reconciled against DESIGN_FINDINGS_2026-08-07.md

**2026-08-09.** `DESIGN_FINDINGS_2026-08-07.md` listed five corrections. None had reached
`PROTOCOL.md`, which is the document that goes to the beach with the hardware. All five were
re-derived from `dilatancy_design.py` before anything was edited — the findings doc was not
taken on trust either — and all five hold. They are now applied.

Test suite re-run after editing: **22 passed**. Only docstrings were touched in the module;
no computed value changed.

---

## 1. The five corrections, as re-verified

Every number below was recomputed from `dilatancy_design.py`, not copied from the findings
doc. Values agree to the digits quoted there.

### C1. Porosity and packing fraction were the same symbol

`PROTOCOL.md` §2 wrote `φ ≈ 0.4` for the **porosity** in Kozeny-Carman while every other
section used `φ` for the **packing fraction** (φ₀ = 0.560 vs 0.604, φ_c, φ_min/φ_max).

Verified: `permeability_kozeny_carman(0.2 mm, n=0.40) = 3.9506e-11 m²` against
`(…, n=0.60) = 3.0000e-10 m²`. Ratio **7.5937**, and the closed form is
(0.6/0.4)³·(0.6/0.4)² = 1.5⁵ = 7.59375. The findings doc's "7.6×" is right; the module's own
comment said "7.5" in one line and 7.6 two lines later, so that was fixed too.

**Holds.** Applied: `n` for porosity, `φ` for packing fraction only, with a Notation block
at the top of the document stating n = 1 − φ.

### C2. The rig never reached the footstep regime

Hardening ratio Δp/(ρ_s V²) = ηψD/(k ρ_s V) — at fixed sand the intruder enters only through
D/V.

- Footstep: D/V = 0.10 m ÷ 1.0 m/s = **0.1 s**.
- Old spec, 20 mm sphere at 2 m/s: D/V = **0.01 s**. Mismatch **10.0×**
  (`match_field_regime` returns `mismatch_factor = 10.0`, `V_lab_required = 0.20 m/s`).
- The old spec's *slowest* drop, 0.5 m/s, still only reached D/V = 0.04 s — 2.5× short.

**Holds.** Applied: speed range corrected, derivation below.

One caveat found during verification and now written into §5 and into the
`match_field_regime` docstring: matching D/V matches the **drainage** group only. The
cavitation demand goes as V·D, which is 0.1 m²/s for a footstep against 0.004 m²/s at the
matched lab point — **25× lower**. Matching both groups forces D = 0.1 m and V = 1 m/s, i.e.
full scale. A real footstep on 0.2 mm sand is cavitation-limited; the matched lab drop is
not. This does not weaken C2 (the 10× mismatch is real and worth fixing) but the protocol
should not claim more than it has.

### C3. Q4 was demoted, and it should not have been

Verified the whole chain:

- Cavitation ceiling V_cav = k(p_atm−p_vap)/(ηψD), ∝ d². At ψ = 0.1, D = 20 mm:
  0.49 / 1.10 / 1.96 / 4.40 / 12.2 / 48.9 m/s for d = 0.10 / 0.15 / 0.20 / 0.30 / 0.50 /
  1.00 mm. Drop heights 12 mm to 122 m. Matches the findings table exactly.
- Strader consistency window: no-cavitation bound ψD ≤ **1.0985e-4 m**; hardening bound
  ψD ≥ **4.7275e-3 m**. Ratio **43.037**, window empty. Closes as (53/850)² = 1/257.
- Cavitation bound alone: ψ ≤ **0.005492** at D = 20 mm, **18.2×** below the assumed 0.1.

**Holds.** Applied: §1 Q4 rewritten from "speculative but the highest upside" to a named aim;
§7.4 rewritten from "worth one afternoon, not worth building the project on" into a
**measurement of ψ** (a departure at V_dep gives ψ = k(p_atm−p_vap)/(η V_dep D); no departure
by 2.5 m/s bounds ψ ≤ 0.078 for 0.2 mm sand, ψ ≤ 0.044 for 0.15 mm); §8 sequencing now
contains a Q4 step; the "honest expectation" paragraph upgraded.

### C4. Percentiles cannot deliver the diameter Kozeny-Carman wants

Verified at d₅₀ = 0.3 mm: a_max ratios d₁₀ vs d₅₀ of 1.1295 / 1.3111 / **1.5888** / 2.0830 /
2.5243 at C_u = 1.2 / 1.5 / **2.0** / 3.0 / 4.0, and Sauter vs d₅₀ of 1.0057 / 1.0283 /
**1.0850** / 1.2274 / 1.3858. d₃₂/d₅₀ = 0.9031 against d₁₀/d₅₀ = 0.5606 at C_u = 2, so the
informative contrast is d₁₀ against d₃₂ (**1.464×**), not either against the median.
Discrimination threshold at 15% scatter: C_u ≈ **1.234** for d₁₀ vs d₅₀, C_u ≈ **1.248** for
d₁₀ vs d₃₂ (recomputed; not in the findings doc). Both independent of d₅₀ — checked at 0.15,
0.30 and 0.60 mm.

**Holds.** Applied: §4.1 now requires the per-sieve retained masses to be archived, gives the
discrete Sauter formula d₃₂ = 1/Σ(wᵢ/dᵢ), specifies √2 sieve spacing, and adds the C_u
admission check. §7.1 makes d₁₀-vs-d₃₂ the primary contrast.

### C5. Sweep ψ, D and the Kozeny-Carman constant, do not quote them

**Holds** — it is a methodological instruction, and it is the one the document most needed,
because §2's single quoted ψ = 0.1 turned out to contradict the project's own founding
observation (see "extra" below).

**Applied:** a "Sweep, do not quote" paragraph closing §2, and a ψ sensitivity table.

---

## 2. The corrected speed range, derived

The constraint is D/V, not V. Two anchors, both from the existing documents:

**Slow end — the field.** A footstep is D/V = 0.10 m ÷ 1.0 m/s = 0.1 s. The matched lab
speed is V = D/(0.1 s). For the 20 mm sphere:

> V = 0.020 m ÷ 0.1 s = **0.20 m/s**,  H = V²/2g = 0.04/19.62 = **2.04 mm**

**Fast end — the literature.** §6 makes reproducing Strader et al. a non-optional acceptance
test, and their reported range runs to 2.5 m/s. That speed is therefore already required:

> V = **2.5 m/s**,  H = 6.25/19.62 = **0.319 m**

### **Corrected range: 0.20 to 2.5 m/s for the 20 mm sphere; drop heights 2.0 mm to 0.32 m.**

D/V spans 0.100 s down to 0.008 s, a factor of 12.5, with the footstep at the top end. The
old range spanned 0.040 to 0.005 s and never touched it.

Three independent checks that this range is well placed, all computed:

| check | value | inside range? |
|---|---|---|
| drained/hardening crossover, 1.0 mm sand, D = 20 mm | 0.764 m/s | yes |
| cavitation ceiling, 0.2 mm sand, ψ = 0.1, D = 20 mm | 1.96 m/s | yes |
| Strader's own published window | 1 – 2.5 m/s | yes (top half) |

For other spheres, V = D/(0.1 s) and H = (10D)²/2g ≈ 5.1 D², so the matched drop height grows
as D²: 0.51 mm at D = 10 mm, 2.04 mm at 20 mm, 8.16 mm at 40 mm. The largest sphere is the
easiest to release cleanly at the matched speed.

---

## 3. Everything changed in PROTOCOL.md

| where | change |
|---|---|
| header | provenance line recording this reconciliation and its five items |
| new **Notation** block | φ = packing fraction, n = porosity, n = 1 − φ, and the 7.6× cost of confusing them |
| §1 Q1 | added the discrimination numbers and the d₁₀-vs-d₃₂ framing (C4) |
| §1 Q4 | "speculative" → "a named aim"; V_cav table; the 43× empty window; ψ ≤ 0.0055; reframed as measuring ψ(φ₀, d) (C3) |
| §1 honest expectation | Q4 named as the result that cannot come out merely confirmatory |
| §2 formula | `φ ≈ 0.4` → `n = 0.40`, full Kozeny-Carman written out, prefactor 9.9 × 10⁻⁴ (C1) |
| §2 table | ρ_s V² 8 → **10.6 kPa** (ρ_s = 2650, not an implied 2000); ratios 12 → **9.6** and 0.5 → **0.4**; k for 1.0 mm 1 × 10⁻⁹ → 9.9 × 10⁻¹⁰; d₅₀ column relabelled d₃₂ |
| §2, new | "Only the group D/V matters" (C2) |
| §2, new | ψ sensitivity table (C5) |
| §2, new | "Sweep, do not quote" (C5) |
| §3.6 | field rod: record diameter and drop height so its D/V can be computed (its D/V is 0.0071 s, 14× short of a footstep) |
| §4.1 | archive the grading curve; Sauter formula; √2 sieve spacing; C_u ≥ 1.25 admission check (C4) |
| §5 impact | **0.5–4 m/s → 0.2–2.5 m/s**, with the derivation, the D/V sweep instruction, the D² scaling of matched drop height, and the V·D caveat (C2) |
| §5 imaging | frame-rate numbers rederived for the new range; slow-end sub-pixel problem and an fps-selection rule |
| §6 | control runs over Strader's 1–2.5 m/s; slow end is ours, so disagreement there is not a rig failure |
| §7.1 | d₃₂ named; d₁₀-vs-d₃₂ made the primary contrast (C4) |
| §7.4 | rewritten as a ψ measurement with bounds, sample ordering, the d² discriminator, and the >2.5 m/s extension rule (C3) |
| §8 | new step 2 (D/V sweep on beads) and new step 6 (Q4), renumbered (C2, C3) |

### Also changed in `dilatancy_design.py` (docstrings only, no computed value)

The wrong numbers had propagated here too, which is why this file was checked:

- notation warning: "cost a factor of **7.5**" → **7.6** (the comment two lines below already
  computed 7.6)
- `consistency_window`: "EMPTY by a factor of **~40**" → **43.0**
- `match_field_regime`: docstring pointed only at the 200 mm sphere, the expensive fix the
  findings doc says *not* to buy. Now gives both options, names the 2 mm drop, and carries
  the V·D caveat.
- module docstring Q4: "roughly 2 m/s … 0.2 m drop" → 1.96 m/s, 0.19 m; notes the reframing
  to a ψ measurement.

`DESIGN_FINDINGS_2026-08-07.md` was left alone. It is a dated record, its numbers all check
out, and rewriting it would destroy the audit trail this reconciliation depends on.

---

## 4. Found during reconciliation, not one of the five

**§2's footstep sentence was wrong.** It read: "The same estimate at footstep scale
(D ~ 0.1 m, V ~ 1 m/s) puts fine beach sand within a factor of a few of the crossover."
Computed at D = 0.1 m, V = 1 m/s: the ratio is **95** for 0.2 mm sand and **3.8** for 1.0 mm,
with the crossover at **d ≈ 1.95 mm**. So it is *coarse* sand that sits within a factor of a
few, not fine, and the sentence was out by a factor of ~30 on the case it named.

Correcting it exposes something better. At ψ = 0.1 the model says every sand finer than
granule gravel hardens underfoot — which **contradicts the field observation in this
document's own Origin line** ("did not harden at all on steeper, coarser stretches"). The
contradiction is a statement about ψ:

| ψ | grain size at which a footstep stops hardening |
|---|---|
| 0.10 | 1.95 mm |
| 0.03 | 1.07 mm |
| 0.010 | 0.62 mm |
| 0.0055 (Strader bound) | 0.46 mm |
| 0.003 | 0.34 mm |

At ψ ≈ 0.005 the crossover lands in fine-to-medium sand, where the beach puts it. So the
barefoot observation and Strader's unbroken power law point at the same ψ, roughly 20× below
the textbook 0.1, from two completely independent directions. That is a second, independent
argument for C3, and it is now in §2.

**The §2 table's fine-sand row sits on the cavitation ceiling.** The demanded Δp is 101.2 kPa
but p_atm − p_vap = 98.985 kPa, so the deliverable suction is 99 kPa and the true ratio is
**9.3**, not 9.6. The table now says so rather than quoting an unreachable demand as if it
were a stress.

**A symbol collision was nearly introduced.** The natural way to specify sieve spacing is
"half-φ intervals" — the sedimentological grain-size scale, a *third* meaning of φ. §4.1 now
specifies "successive openings in a ratio of √2" and says explicitly why.

---

## 5. Nothing failed verification

All five corrections reproduce from the code to the precision quoted. The one qualification
is on C2, recorded above and now in the protocol: matching D/V matches the drainage regime
but under-drives the cavitation demand by 25×, and no bench-scale rig can match both groups
at once.

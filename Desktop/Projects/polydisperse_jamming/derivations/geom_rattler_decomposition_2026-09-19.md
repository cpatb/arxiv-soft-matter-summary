# Item B — The rattler decomposition: can depopulation alone carry the CV effect?

**VERDICT: UNDETERMINED — and the reason is a specific, cheap, missing measurement.**

**What this anchor can and cannot decide.** It CAN decide whether backbone invariance is
*impossible*, because the exact rearrangement bound on the rattler volume closes that
question with no free parameter (it did not fire once the test was run seed-paired). It
CANNOT decide whether backbone invariance is *true*: that needs the rattler VOLUME
fraction ψ, and **no shard records it — every rattler field on disk is a COUNT.** It says
nothing about *why* a backbone would be size-blind.

Every number below is printed by `geom_rattler_decomposition_2026-09-19.py` →
`geom_rattler_decomposition_2026-09-19.txt`.

---

## 1. The right decomposition is not the obvious one

The brief suggests `φ_J = φ_backbone + φ_rattler`. That is fine as bookkeeping but it
hides the physics, because the two terms are not independent: rattlers sit *inside* the
backbone's pore space. The useful exact form is a ratio.

Let `V` be the box, `V_B` the backbone solid volume, `V_R` the rattler solid volume. Define

    φ   = (V_B + V_R)/V        the measured packing fraction (rattlers included)
    φ_B = V_B/V                the backbone solid per unit TOTAL volume
    ψ   = V_R/(V_B + V_R)      the rattlers' share of the SOLID

Then, and SymPy asserts the difference to exactly zero (§1):

> **φ = φ_B / (1 − ψ)**   *(identity, no model)*

`φ_B` is the right variable because it is what stays fixed if the load-bearing skeleton
is geometrically unchanged: it measures backbone solid against the *same* box.

### Which moment of the PSD enters

With N grains and rattler set R of size qN,

    ψ = Σ_{i∈R} d_i³ / Σ_i d_i³ = q · ⟨d³⟩_R / ⟨d³⟩.

**The functional the depopulation picture needs is the third moment CONDITIONED ON THE
RATTLER SET** — not ⟨d³⟩, not d₃₂, not d₄₃. That answers the brief's question directly,
and it is the reason the picture is not testable from what is on disk: the shards record
`rattler_frac`, a *number* fraction, and the conditional moment is exactly the part that
is not a number.

At CV = 0 every grain is identical, so `⟨d³⟩_R = ⟨d³⟩` and **ψ = q exactly**. The
monodisperse cell is a zero-width anchor.

### The exact bound

By the rearrangement inequality, for fixed q the rattler volume is minimised by taking
the q smallest grains and maximised by the q largest. Define the PSD's **volume-share
function**

    Ψ(q) = ∫₀^{F⁻¹(q)} d³ f(d) dd / ⟨d³⟩,

so that

> **Ψ(q) ≤ ψ ≤ 1 − Ψ(1−q).**

Closed forms (§1, derived, plus a 10⁶-draw Monte-Carlo check agreeing to ≤ 3e−3):

| family | Ψ(q) |
|---|---|
| lognormal | `Φ(Φ⁻¹(q) − 3σ)`, `σ² = ln(1+CV²)` |
| uniform on [1−w, 1+w], `w = CV√3` | `(t⁴ − (1−w)⁴)/((1+w)⁴ − (1−w)⁴)`, `t = 1−w+2wq` |
| bidisperse, atoms `1∓CV`, 50:50 by number | `2q·w_s` for `q ≤ ½`; `w_s + (2q−1)(1−w_s)` for `q > ½`, with `w_s = (1−CV)³/[(1−CV)³+(1+CV)³]` |

Limiting checks asserted: Ψ(0) = 0, Ψ(1) = 1, Ψ(q) ≤ q with equality iff CV = 0.

---

## 2. The zero-parameter gate

**Hypothesis H1 (backbone invariance).** `φ_B` depends on friction and protocol but not
on the size distribution. Motivation: the project's fact 3 — `Z_nr` and the fabric
invariant T are flat in CV.

Because the CV = 0 interval collapses to a point, the gate is one-sided and sharp:

> **H1 is REFUTED at a cell iff  φ(CV)·(1 − Ψ(q(CV)))  <  φ(0)·(1 − q(0)),**
> i.e. iff even making every rattler one of the *smallest* grains leaves the backbone
> less dense than the monodisperse one.

That gate can fire. §4 prints, for every (rung, family), the value of φ(CV = 0.5) below
which it *would* have fired — margins run +0.0150 to +0.0588.

**Result on seed means (§3):** it fires in 1 of 12 (rung, family) blocks — bidisperse at
p̂ = 1e−3, CV = 0.2 and 0.3, by −0.0011 and −0.0034.

**Result seed-paired (§3b):** those margins are 0.0011–0.0034 in φ against this
protocol's own monodisperse run-to-run scatter of 2.1e−3 (`paper/numbers.tex
\phiScale`). Re-run per seed:

| rung | family | CV | M (margin) | se | t | verdict |
|---|---|---|---|---|---|---|
| 1e−3 | bidisperse | 0.2 | −0.00105 | 0.00283 | −0.37 | negative, NOT resolved |
| 1e−3 | bidisperse | 0.3 | −0.00335 | 0.00321 | −1.04 | negative, NOT resolved |

**0 of 60 cells has a resolved negative margin at t_crit = 4.303.** So H1 is not
impossible. It is also not confirmable from disk. **UNDETERMINED.**

---

## 3. What H1 requires, if it is true — the number a new run must hit

Invert for the required conditional volume ratio `R = ψ_req/q = ⟨d³⟩_R/⟨d³⟩` (§5). At
p̂ = 3e−5:

| family | CV 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|
| lognormal | 0.8046 | 0.5723 | 0.4446 | 0.3828 | **0.3274** |
| uniform | 0.8138 | 0.5732 | 0.3799 | 0.3003 | **0.2620** |
| bidisperse | 0.7981 | 0.5444 | 0.3590 | 0.2690 | **0.2074** |

Equivalently, rattlers must have an effective diameter `R^{1/3}` of 0.59–0.69 of the mean
at CV = 0.5. For the bidisperse family the PSD has two atoms, so R fixes the composition
of the rattler set *exactly*: at CV = 0.5, p̂ = 3e−5, H1 requires **92.4 % of the small
grains and 7.3 % of the large grains to be rattlers**. That is a directly checkable
statement.

### Size selectivity — an emergent near-constant

Define `S = (q − ψ_req)/(q − Ψ(q))`, which interpolates between

* `S = 1` — *size-ordered*: the rattlers are exactly the q smallest grains, which is what
  a pure coordination threshold `z(d) < 2` gives, since expected neighbour count grows
  monotonically with d;
* `S = 0` — *size-blind*: rattling is independent of size.

`S` outside [0, 1] refutes H1. Measured requirement over 60 cells (§5b):

| p̂ | family | CV 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| 1e−3 | lognormal | 0.553 | 0.689 | 0.836 | 0.864 | 0.848 |
| 1e−3 | uniform | 0.588 | 0.807 | 0.861 | 0.893 | 0.926 |
| 1e−3 | bidisperse | 0.547 | **1.024** | **1.028** | 0.966 | 0.944 |
| 3e−4 | lognormal | 0.401 | 0.674 | 0.778 | 0.803 | 0.790 |
| 3e−4 | uniform | 0.414 | 0.725 | 0.817 | 0.876 | 0.900 |
| 3e−4 | bidisperse | 0.748 | 0.882 | 0.938 | 0.909 | 0.885 |
| 1e−4 | lognormal | 0.125 | 0.464 | 0.766 | 0.790 | 0.815 |
| 1e−4 | uniform | 0.524 | 0.659 | 0.783 | 0.832 | 0.873 |
| 1e−4 | bidisperse | 0.638 | 0.816 | 0.858 | 0.875 | 0.857 |
| 3e−5 | lognormal | 0.460 | 0.671 | 0.750 | 0.765 | 0.806 |
| 3e−5 | uniform | 0.450 | 0.635 | 0.791 | 0.835 | 0.856 |
| 3e−5 | bidisperse | 0.691 | 0.839 | 0.878 | 0.856 | 0.854 |

Range 0.125 to 1.028, median 0.811; 2 cells (the bidisperse pair above) are the §3 near-
refutation. **For CV ≥ 0.3 every one of the twelve (family, rung) blocks requires
S = 0.75–1.03.** If H1 is right, rattling is strongly but not perfectly size-selective,
and the degree of selectivity is nearly a constant across three families, five CVs and
four pressure rungs.

That near-constancy is a **prediction, not a measurement**. It is what H1 requires.

---

## 4. The bidisperse non-monotonicity is the sharpest stress test — and H1 passes it

φ_J for the bidisperse family *peaks near CV = 0.3 and falls*, while its rattler number
fraction rises monotonically 0.04 → 0.50. Under H1 the required ψ must therefore be
NON-monotone even though q is monotone. §6 confirms that it is, at every rung:

| rung | q monotone? | ψ_req monotone? | ψ_req peaks at |
|---|---|---|---|
| 1e−3 | True | **False** | CV = 0.4 |
| 3e−4 | True | **False** | CV = 0.3 |
| 1e−4 | True | **False** | CV = 0.3 |
| 3e−5 | True | **False** | CV = 0.3 |

H1 does not merely tolerate the non-monotonicity: it **requires** a non-monotone
conditional rattler volume, produced by `R(CV)` falling faster than `q(CV)` rises and then
reversing. Nothing in the data contradicts that, and nothing on disk confirms it.

---

## 5. The measurement that settles it

Field inventory of a shard record (§7, printed from
`out/phase2_ln0_shard00of06.jsonl`):

* top level: `converged, cv, cv_realised, d10_realised, d32_realised, d43_realised,
  d_max_over_L, d_span, family, force, iters, key, mu_p, n_contacts, n_placed,
  n_requested, notes, p_hat, p_hat_target, phi, protocol, rattler, rattler_frac, seed,
  stress_dev_hat, unbalanced, wall_s, z_all, z_nr`
* `rattler` sub-dict: `largest_frac_k2/k3/k4, rattler_frac_k2/k3/k4, z_nr_k2/k3/k4`

Every rattler field is a **count**. The three size fields (`d10/d32/d43_realised`) are
whole-pack moments, none conditioned on the rattler set. **ψ is not recoverable from these
files.**

> **The cheapest decisive measurement in this whole item is one extra scalar per cell:**
> `Σ_{i∈rattler} d_i³ / Σ_i d_i³`. The rattler set is already constructed
> (`packing.iterated_rattlers`), the radii are already in memory, and the k = 2/3/4
> convention already exists — so it is `rattler_vol_frac_k2/k3/k4` beside the three
> `rattler_frac_k*` entries, at zero simulation cost on any future run.

With ψ in hand, `φ·(1−ψ)` becomes directly comparable across CV and the item closes in one
analysis pass, with no new DEM cells for cells that are re-run.

### Prior art — the distinction is published, the curve is not

A citation pass on 2026-09-19 corrects an overclaim I was about to make. **The
number-versus-volume distinction for rattlers IS in the literature**, and two results
support the size-selectivity picture directly:

* **Biazzo, Caltagirone, Parisi & Zamponi**, *J. Chem. Phys.* **132**, 176101 (2010),
  doi:10.1063/1.3428387 — Fig. 2 plots the per-species rattler *number* fraction with an
  **inset plotting the rattler VOLUME fraction**. Text: *"the fraction of rattlers among
  the small particles is equal to 1 … the total volume fraction of rattlers is always
  quite small and rarely exceeds 5 % of the total."* Bidisperse; axes are (size ratio,
  composition), not a CV.
* **Monti, Clemmer, Srivastava, Silbert, Grest & Lechman**, *Phys. Rev. E* **106**,
  034901 (2022), doi:10.1103/PhysRevE.106.034901 — states the disparity outright:
  *"while rattlers may constitute a large fraction of the total NUMBER of particles, they
  typically only account for a small fraction of the total particle VOLUME"* (number
  fraction to 97 %, yet φ − φ_nr < 0.1). **This is independent external support for a
  high size selectivity S.**
* **Cantor, Azéma, Sornay & Radjaï**, *Phys. Rev. E* **98**, 052910 (2018),
  doi:10.1103/PhysRevE.98.052910 — floater proportion resolved *per size class* against
  the span S, the cleanest published demonstration that small particles preferentially
  rattle as the distribution widens.
* **Petit, Kumar, Luding & Sperl**, *Phys. Rev. E* **106**, 054903 (2022),
  doi:10.1103/PhysRevE.106.054903 — the mirror image: per-species non-rattler *volume*
  across the (δ, X_S) plane, and no rattler number fraction anywhere.

What was **not** found published is the specific curve this item needs: a *normalised
rattler volume fraction* `ψ = V_R/V_solid` plotted against the CV (or Δ, or span) of a
**continuous** size distribution — the curve that makes the number/volume divergence
visible as a function of *width*. Every width-axis paper located (Hoy 2025
doi:10.1103/3m8z-qxst, Zhang & Scheffold 2015 doi:10.1103/PhysRevE.91.032302, Cantor 2018,
Göncü/Durán/Luding doi:10.1016/j.crme.2010.10.004) reports the **number** fraction only.
Desmond & Weeks contains zero occurrences of "rattler".

So the correct statement is not *"nobody has measured ψ"* but **"ψ against distribution
width, for a continuous PSD at finite friction, was not found published"**, with Biazzo
(2010), Monti (2022) and Petit (2022) cited as prior art establishing that the
distinction matters. Two residual risks are **UNRESOLVED** (full texts unreachable):
Oquendo-Patiño & Estrada, doi:10.1007/s10035-020-01043-9 and
doi:10.1103/PhysRevE.105.064901; and Wang, Shi & Ling, doi:10.1631/jzus.A2500203. An ILL
check on those three is required before any novelty claim is made in print.

---

## 6. A note on a number quoted elsewhere

`MECHANISM_2026-09-08.md` §M3 writes the rattler fraction as "0.040 → 0.234". On the
N = 1500, μ_p = 0.3 table printed in §2 of the script, (0.040, 0.234) is the p̂ = 1e−3
rung at **CV = 0 and CV = 0.3**, not CV 0 → 0.5; the CV 0 → 0.5 range at that rung is
0.0393 → 0.4476 (lognormal), 0.0404 → 0.4576 (bidisperse), 0.0404 → 0.4231 (uniform).
`paper/main.tex:162` uses the correct endpoints (0.039 → 0.448, "about elevenfold"), so
this is a stale line in an older document, not a manuscript error.

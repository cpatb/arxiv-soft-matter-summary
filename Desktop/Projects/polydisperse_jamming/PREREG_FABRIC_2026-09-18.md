# PREREG — the fabric rerun, 2026-09-18

Registered **before** `grade_fabric.py` exists and **before** any cell of the campaign is
run. Sizing numbers are printed by `diag_fabric_null.py` → `out/diag_fabric_null.txt`
(committed with this file; it runs no cell).

---

## 0. What this is, and the correction to the framing it inherited

The Status card says fabric "tests directly" the second half of the rattler result —
polydispersity **removes** grains from the load-bearing network (rattler fraction
0.040 → 0.234 across CV) **without restructuring** it (backbone Z barely moves). Two
things had to be fixed before that sentence could be graded, and both change the design:

**(a) The deviator itself carries no CV signal that symmetry permits.** These packs are
compressed under `PeriTriaxController` with `goal = (g, g, g)`, i.e. isotropically. The
ensemble-mean deviatoric fabric is then zero at every CV by symmetry. What is *not* fixed
by symmetry is the deviator's **magnitude relative to the finite-sample floor**, which is a
rotation-invariant scalar. So the registered observable is

    T = n_contacts · A²,        A = fabric_tensor(...)["aniso"]

For `n` **independent** unit normals E[A²] = 1/n exactly (derivation, every step, in
`diag_fabric_null.py`; Monte Carlo reproduces E[T] = 1 to 2.0 %, sd(T) = 0.6375 against the
predicted √(2/5) = 0.6325). T > 1 is orientational **correlation** in the contact network.

**(b) A raw `aniso` vs CV trend is manufactured by the contact count.** `n_contacts` falls
3683 → 2050 (lognormal), 3684 → 2184 (uniform), 3684 → 2039 (bidisperse) across CV at
`p_hat` = 1e-3, so the floor alone rises **0.0157 → 0.0210**, a slope of **+0.0107 per unit
CV**, *in the same direction a restructuring signal would*. T removes that by construction.
`aniso` must never be read against CV.

**(c) An all-contact fabric cannot answer the question, so the code was changed before the
run.** The rattler fraction is the thing that moves; a fabric over every real interaction
therefore mixes a change in the rattler **population** with a change in the backbone's
**structure**. `packing.contact_normals(O, with_ids=True)` + `packing.backbone_mask` now
also record `fabric_nr`, over contacts whose **both** ends are non-rattlers — **at every
threshold in `RATTLER_THRESHOLDS = (2, 3, 4)`**, as `aniso_k2/k3/k4`, with the primary
`RATTLER_MIN_CONTACTS = 2` aliased to the bare keys. Recording all three is this project's
own convention for the rattler ladder and it is not decoration here: **at k = 2 a rattler
is a grain with ≤ 1 contact**, so in the F1 gate cell (bidisperse CV 0.3) removing 404 of
1500 grains removed **2 of 2762 contacts**, and the k = 2 backbone fabric is all but
identical to the all-contact one. The conflation this split exists to avoid is
**negligible at k = 2 and real at k = 4**, which is a fact about the definition, not about
the packing. **`fabric_nr` at k = 2 is the primary channel; k = 3 and k = 4 are reported
beside it and a verdict that holds at only one k is labelled a statement about the
definition.**

**Honest limit, registered up front.** E[T] = 1 is the null for *independent* normals, and
a contact network is not independent: contacts sharing a grain are geometrically
constrained, and `z_all` itself falls with CV (≈ 4.5 → 3.1). So **the absolute value of T
carries no claim**, and a T drift within a family is degenerate with a z drift. The
campaign is a **comparative** measurement; F7 below is the only thing that can break that
degeneracy, and it can only be read if F2 finds a drift.

---

## 1. The grid, and the tag guard

| | |
|---|---|
| runner | `phase2_phiJ.py`, **unchanged** apart from the `--seeds` fix in §5 |
| families | lognormal, uniform, bidisperse |
| CV | 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 |
| friction | `mu_p` = 0.3, **matched** |
| N | 1500 |
| primary rung | `p_hat` = 1e-3, **12 seeds** (1–12) → 3 × 6 × 12 = **216 cells** |
| secondary rung | `p_hat` = 1e-4, **6 seeds** (1–6) → 3 × 6 × 6 = **108 cells** |
| tag | **`fab`** (primary), **`fabgate`** (the F1 cells) |

**This is not a φ_J result and must never be quoted as one.** A single pressure rung
cannot be extrapolated to `p_hat` → 0; the campaign reads fabric against CV at matched
pressure and friction, and nothing else.

**12 seeds is not a round number, it is the design.** `out/diag_fabric_null.txt` §3: at
sd(T) = 0.6375 and six CV rungs, 2·se(slope) expressed as a change across CV 0 → 0.5 is
0.880 / 0.622 / 0.508 / **0.440** / 0.394 at 3 / 6 / 9 / **12** / 15 seeds. 12 is the
smallest count that resolves the registered effect size.

**Registered effect size: ΔT = 0.5 across CV 0 → 0.5.** A declared convention, not a
measurement: T = 1 is what an orientationally structureless network returns, so this is
"half the whole floor's worth of extra correlation, appearing across the ladder".

**TAG GUARD (structural, asserted in the grader, not by eye).** Every record must come
from a file whose basename starts with `phase2_fab`; the grader exits non-zero if any input
file matches `mixr*`, `extr*`, `tern*`, `kl*`, `ug*`, or if any admitted record lacks a
non-empty `fabric` **and** `fabric_nr`. Pre-2026-09-18 records have neither and cannot
pool with these by construction. N is printed before and after every filter.

---

## 2. Gates

Each gate is a **sentence** first. The grader prints the sentence beside the number.

### F0 — reproduction, runs FIRST and exits non-zero on mismatch
*Every admitted record's `fabric["n"]` equals its own `n_contacts`, and at **each**
threshold k ∈ {2, 3, 4} its `fabric_nr["n_k{k}"]` equals the backbone contact count implied
by observables recorded independently of the fabric code:*

    n_nr_contacts(k) = z_nr_k{k} · (1 − rattler_frac_k{k}) · n_placed / 2

*to within 0.5 of an integer (both sides are counts; the right-hand side is a ratio of
integers and the rounding is the only slack allowed).*
**k = 4 is DEGENERATE at `mu_p` = 0.3 and is registered as such before the run**: mean
`z_all` here is 3.1–3.9, below the threshold, so iterated removal consumes the whole pack
and the backbone is EMPTY. All three F1 gate cells returned `n_k4 = 0` with
`z_nr_k4 = NaN`. The gate therefore requires the two sides to be empty *together* — a
count on one side and a NaN on the other is a failure — and **k = 4 carries no F2
verdict**. That leaves k = 2 (primary) and k = 3.
This is not circular: `z_nr_k*`, `rattler_frac_k*` and `n_placed` are written by
`rattler_report` and the PSD insertion, none of which sees a contact normal. It can fail —
a `backbone_mask` that used `or` instead of `and`, or that misaligned ids against normals,
breaks it immediately. Checked at all three k, so a threshold wired to the wrong rattler
set cannot hide behind the one where the two sets nearly coincide.

### F1 — the solver must not move (bit identity) — ALREADY RUN, before this file
F1 was run before this prereg was committed, and that is deliberate: it is a **wiring
gate on three cells that already existed on disk**, it decides nothing about CV, and a
campaign launched on a code path that moves the solver would be worthless. Its result is
recorded here with the rest. **No production cell had been run when this file was
committed.**

*Three cells that already exist on disk — bidisperse CV 0.3, uniform CV 0.3, lognormal
CV 0.4, all seed 1, `mu_p` = 0.3, `p_hat` = 1e-3, N = 1500 — re-run under the fabric code
return `phi` **bit-identical** to the committed record.*
**`phi` is the gate.** It is the continuous channel: `z_all` (2·k/N), `rattler_frac` (k/N)
and `z_nr` are grid-quantized and a real perturbation can hide inside one quantum
(lessons.md 2026-08-31). They are printed beside it, not read as the gate.
Reference values, read from the committed shards now so the comparison cannot drift:

| cell | file | `phi` |
|---|---|---|
| bidisperse CV 0.3 | `phase2_bidisperse_shard00of12.jsonl` | `0.6120403164804106` |
| uniform CV 0.3 | `phase2_uniform_shard00of12.jsonl` | `0.6101301588958694` |
| lognormal CV 0.4 | `phase2_shard00of12.jsonl` | `0.6339122825028567` |

(The lognormal value is independently confirmed: `phase2_capfix_lg3_shard00of07.jsonl` ran
the same cell at a different iteration cap and returned the identical 16 digits.)

**OUTCOME, already in hand: F1 PASSES 3/3.** `phi` came back bit-identical on all three
cells; `z_all`, `z_nr` and `rattler_frac` did too, and are printed but not read as the
gate. Shards: `out/fabric/phase2_fabgate*.jsonl`. `grade_fabric.py` re-checks this from
the committed shards rather than taking this sentence for it.

### F2 — PRIMARY. The registered prediction is the NULL
*If polydispersity only depopulates the load-bearing network, the orientational
correlation of the surviving backbone does not change with CV: the per-family OLS slope of
`T_nr` against CV is consistent with zero at 2·se, on all three families.*
Read at **k = 2** (primary) and reported at **k = 3**; k = 4 is empty at this friction
(F0) and carries no verdict. A verdict that changes between k = 2 and k = 3 is reported as
a statement about the rattler definition.

Verdict grammar, registered in full so no reading is chosen after the fact:

- `|slope| > 2·se` on any family → **DRIFT**, the null is refuted on that family.
- `|slope| ≤ 2·se` **and** the *achieved* 2·se·0.5 ≤ 0.5 → **NULL, INFORMATIVE**.
- `|slope| ≤ 2·se` **and** the achieved 2·se·0.5 > 0.5 → **UNINFORMATIVE**. Not a pass.

`se` is recomputed from the **measured** scatter, never from the null's 0.6375, because
real orientational structure would make the scatter larger and a design sized on the null
would then claim a resolving power it does not have.

### F3 — the all-contact channel, reported beside it
*The same slope on `T_all`, and the paired per-cell contrast `T_all − T_nr`.*
No prediction registered. It is printed because `T_all` is the quantity that conflates the
rattler population with the backbone, and the contrast is the size of that conflation.

### F4 — the confound must be demonstrated, not asserted
*The raw slope of `aniso` against CV is statistically indistinguishable from the slope the
finite-sample floor produces on its own (0.95/√n at each cell's own `n_contacts`), on all
three families, at 2·se of the difference.*
If the raw slope **exceeds** the floor slope significantly, there is a real deviatoric
signal under isotropic compression and the symmetry argument in §0(a) is wrong — a finding,
and one that would have to be reported before anything else in this file is believed.

### F5 — rung consistency (cannot decide F2)
*On the 6 seeds shared by both rungs, the paired difference `T_nr(1e-3) − T_nr(1e-4)` is
reported with its own se.* At 6 seeds this rung cannot resolve the registered effect
(2·se·0.5 = 0.622 > 0.5) and is registered as a **consistency channel only**.

### F6 — housekeeping that has bitten this project before
*Cells with `converged = False` are excluded, with N printed before and after. The
`d_max/L > 0.25` finite-size flag is reported **broken down by family**, never pooled
(memory: it fires on 100 % of lognormal cells and 0 % of zero-skew ones, i.e. entirely on
one family).* Any family whose cells are wholly flagged has its F2 verdict labelled as
resting on flagged cells.

### F7 — CONDITIONAL, readable only if F2 finds a drift
*Does `T_nr` collapse on `z_nr` across the three families at matched CV?* Within a family
z and CV are collinear by construction, so a drift in T cannot be attributed to either;
the three families have **different** z at matched CV (the 4.5× shape result), which is the
only leverage available here. A collapse on z says the drift is depopulation seen again;
a failure to collapse says it is not.

---

## 3. What a drift would and would not license

**Would.** That the load-bearing network's orientational statistics change with CV under
isotropic compression — i.e. that the "removes grains without restructuring" reading of the
rattler result is **false as stated**, subject to F7's z degeneracy.

**Would not.** Anything about φ_J (one pressure rung, no extrapolation). Anything about
sheared or impacted packs — fabric under shear is a different quantity with a different
symmetry, and this campaign says nothing about it. Anything about the drop-impact side.
Anything about *which* mechanism moves it. And nothing about the deviator's **direction**,
which the ensemble symmetry pins at zero regardless of what one realisation shows.

---

## 4. What is pre-computed, so it cannot be tuned afterwards

From `out/diag_fabric_null.txt`, committed with this file:

- closed form E[T] = 1, MC-confirmed to 2.0 %, sd(T) = 0.6375 (predicted 0.6325);
- floor slope with CV: +0.0107 (lognormal), +0.0094 (uniform), +0.0108 (bidisperse) at
  `p_hat` = 1e-3; +0.0140 / +0.0109 / +0.0115 at 1e-4;
- 2·se(slope)·0.5 by seed count: 0.880 / 0.622 / 0.508 / 0.440 / 0.394 at 3 / 6 / 9 / 12 / 15.

---

## 5. Code changed before the run, and why each change was necessary first

1. **`packing.contact_normals(O, with_ids=True)` + `packing.backbone_mask` +
   `PackingState.fabric_nr`** — without them the campaign's primary channel does not exist
   and 324 cells would have measured the conflated quantity (§0c).
2. **`phase2_phiJ.py`: the `--seeds` override moved AFTER the default.** Until now the two
   lines were the other way round, so `--seeds` was accepted and silently overwritten.
   **No record on disk is affected** — every campaign that passed the flag passed `1,2,3`,
   which is `SEEDS` — but this campaign needs 12 seeds and would silently have run 3.

*Note on the date: this file follows the repo's running date convention (the Status card's
2026-09-18 block); the host clock reads 2026-09-16.*

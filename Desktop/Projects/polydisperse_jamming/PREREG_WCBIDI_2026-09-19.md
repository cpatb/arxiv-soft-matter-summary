# PREREG — the frictionless bidisperse block, and whether the Wyart-Cates translation holds for the non-monotone family

**Written 2026-09-19, BEFORE any cell of this campaign is launched.** Every number quoted
as a prediction below is computed from data already on disk and is recorded here so the
comparison cannot be made after the fact. The immediate cause is referee pass 3 on
`derivations/wc_endpoints_move_together_2026-09-19.md`, and my own §3 caveat (d) there:
*"TWO families, not three. Bidisperse — the non-monotone one — has no frictionless block, so
the family whose behaviour is least like the others is the one not tested."*

---

## 1. What is being run

| | |
|---|---|
| runner | `phase2_phiJ.py`, **unchanged** |
| family | `bidisperse` |
| friction | `mu_p = 0.0` **only** (the frictional block already exists) |
| CV | 0.0, 0.1, 0.2, 0.3, 0.4, 0.5 — the full production ladder |
| seeds | 1, 2, 3 |
| `p_hat` | the production ladder 1e-3, 3e-4, 1e-4, 3e-5 |
| N | 1500 |
| tag | **`wcbidi`** — matches no existing tag (`capfix`, `shpcap`, `fs_n3k`, `commonlad`, `ratem`, `mixr*`, `extr*`, `tern*`, `kl*`, `ug*`) |
| records | **72** = 6 CV × 3 seeds × 4 rungs |
| invocation | `scripts/launch_wcbidi.sh`, invoked **by file path from PowerShell**, held open in a background Bash call for the whole run |

**Cost, from the wall times already on disk.** Frictionless cells are the expensive ones:
median `wall_s` per record at `mu_p = 0` is 487–724 s by CV rung, with a tail to 4587 s.
72 records × ~600 s median ≈ **12 CPU-hours**; at 12 single-threaded workers that is
**≈ 1–2.5 h wall**, tail-dominated. This is the reason the block was never run.

**Launch discipline** (the failure mode that killed three campaigns on this machine).
`scripts/launch_phase2.sh`'s own docstring: *"`wsl.exe -e bash -lc '... &'` returns as soon
as the launcher shell exits, and WSL then kills the process group… Keeping this script in
the foreground of a held-open background Bash call keeps the whole group alive."* The
launcher therefore **waits** on its workers; the Bash call that starts it is backgrounded
and not allowed to return.

---

## 2. Why this campaign and not another

`wc_endpoints_move_together_2026-09-19.md` reports that the two Wyart-Cates endpoints move
together under polydispersity, so the `φ` axis translates rather than deforming. It runs on
**lognormal** and **uniform**. Both of those families are **monotone** in CV. The bidisperse
family is **not** — `main.tex` states it peaks near CV = 0.3 — and it is the family with no
frictionless block.

So the translation claim has never been tested where it is most likely to fail.

---

## 3. Gates, in the order they are graded

### W1 — WIRING, and it can fail

At CV = 0 the generator ignores the family label: all three families draw the **same**
monodisperse packing (established 2026-09-19 on the `mu_p = 0.3` block, where every rung
with a single record is bit-identical across labels). The same must hold at `mu_p = 0`.

**Registered:** the 12 new `(bidisperse, CV = 0, mu_p = 0)` records must be **bit-identical**
on `phi`, `z_all`, `z_nr`, `rattler_frac`, `p_hat` and `n_contacts` to the existing
`(lognormal, CV = 0, mu_p = 0)` and `(uniform, CV = 0, mu_p = 0)` records at the same
(seed, rung).

⚠ **This is a gate, not a formality.** If it fails, the family label is doing something at
zero width, and the shared-cell finding that this week's manuscript fix rests on is wrong.

Existing values it must reproduce (`out/phase2_fits.csv`, extrapolated):
`phi_J_lin(mu_p = 0, CV = 0) = 0.63237` for both lognormal and uniform.

### G1 — ladder health

**Registered prediction: 72 of 72 records converge** (`converged = true`) under the default
`MAX_ITER`. The existing `mu_p = 0` blocks for the other two families converged on the same
grid, so a failure here is a property of the bidisperse PSD, not of the protocol.
If any cell caps, it is reported as an amendment and **not** re-run at a raised cap inside
this campaign — a raised resource limit is a change of protocol (standing rule, S6).

### G2 — the shift comparison, **registered NOT INFORMATIVE IN ADVANCE**

The natural test is "is `Δφ_0` resolved apart from `Δφ_m` for this family?". It **cannot
fail**, and that is registered here rather than discovered afterwards:

- `Δφ_m^bidi` = **+0.01389 ± 0.00248** (on disk).
- The `mu_p = 0` standard errors for the other two families are 0.00268 and 0.00162, so a
  propagated `se_diff ≈ 0.00329`.
- `t_crit(2 dof) × se_diff = 0.01414`, which is **102 % of the effect itself**.

So any outcome passes. **G2 will be printed and labelled NOT INFORMATIVE regardless of its
value.** It is listed so that a pass is not read as support.

### G3 — **PRIMARY, and it can fail: does the frictionless branch inherit the non-monotonicity?**

`φ_m^bidi(CV)` on disk, mean over 3 seeds:

| CV | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| `φ_m` | 0.57514 | 0.57653 | 0.58717 | **0.59661** | 0.59409 | 0.58902 |
| s.e. | 0.00081 | 0.00116 | 0.00087 | 0.00173 | 0.00066 | 0.00201 |

argmax at **CV = 0.3**; peak-to-endpoint drop **0.00759**.

**Two candidate outcomes, both stated before the run:**

- **(T) TRANSLATION.** `φ_0^bidi(CV)` is also non-monotone with an interior maximum, and the
  window stays flat. Registered test: **argmax of `φ_0^bidi(CV)` ∈ {0.2, 0.3, 0.4}**.
- **(D) DEFORMATION.** `φ_0^bidi(CV)` is monotone increasing, like the lognormal and uniform
  frictionless branches. Registered test: **argmax = 0.5**.

**Why an argmax and not a t-test on the drop.** The two-point test is underpowered here by
the same arithmetic as G2 (`t_crit × se ≈ 0.013` against a 0.00759 drop). The argmax is a
**six-point shape** statement and the two candidates put it at opposite ends of the ladder,
so it discriminates where the two-point test cannot. Its weakness is registered too: with
three seeds per rung an argmax can move by noise, so **the runner-up rung is reported beside
the argmax**, and if the top two rungs are within 1 s.e. the gate returns **UNDETERMINED**
rather than a verdict.

**Resolving power, computed in advance.** Under (D) the expected argmax is 0.5 and under (T)
it is 0.3. On the frictional block the same statistic sits at 0.3 with the 0.5 rung 0.00759
(≈ 3.8 s.e.) below it, so the shape is visible at this seed count on the channel that has
already been measured.

### G4 — a blind point prediction on the shift

If the endpoint ratio measured on the other two families carries to this one:

`Δφ_0^bidi / Δφ_m^bidi` = **0.9349** (lognormal) or **0.9854** (uniform), so

> **REGISTERED: `Δφ_0^bidi` lands in [+0.0125, +0.0146].**

Predicted `se ≈ 0.0033`, so this is roughly a ±1 s.e. band and it can miss. A miss is
reported as a miss; the band is not widened afterwards.

### G5 — the window

`W(CV) = φ_0(CV) − φ_m(CV)` for the bidisperse family.

> **REGISTERED: `W^bidi(CV = 0) = 0.05723`, EXACTLY**, because at CV = 0 the frictionless and
> frictional cells are both shared with the uniform family (this is W1 restated and is not an
> independent result).
>
> **REGISTERED: the range of `W^bidi` across the six CV rungs is below 2 × the propagated
> paired s.e.**, which is the bar the other two families met (1.75 and 1.06).

⚠ A range over six correlated points against a two-point standard error is not obviously the
right yardstick — referee pass 3 raised exactly this on the companion page. **Both** are
therefore reported: the range against 2 s.e. as above (for comparability with the published
two families), **and** a per-rung paired comparison `W(CV) − W(0)` with its own s.e., which
is the statistic a referee would ask for. The second is registered as the one to quote if
they disagree.

---

---

## 3b. AMENDMENT 1, written while the campaign is running and BEFORE any of its records were read

The `mu_p` convergence analysis promised in §4 below was done from data already on disk while
the 72 cells were still in flight. It changes what G5 can mean, so it is recorded here rather
than used to reinterpret the result afterwards. **No `wcbidi` record had been opened when this
was written**; the campaign's own output files are untouched by it.

**`mu_p = 0.3` is NOT Wyart-Cates's frictional endpoint.** Paired on seed at CV = 0,
lognormal:

| comparison | difference | s.e. | \|t\| vs 4.303 |
|---|---|---|---|
| `φ_J(0.3) − φ_J(3.0)` | **+0.01711** | 0.00120 | **14.29 — RESOLVED** |
| `φ_J(1.0) − φ_J(3.0)` | +0.00161 | 0.00058 | 2.78 — not resolved |

`mu_p = 0.3` captures only **77.2 %** of the total `mu_p = 0 → 3.0` drop; **22.8 % remains
above it**, at 16.5 s.e. `mu_p = 1.0` is converged within this design's power.

**Consequences for this campaign, stated now:**

1. **G5's window is a `mu_p = 0.3` window, not a Wyart-Cates window.** It is still a
   well-defined quantity and the gate stands as written, but the page it feeds must not call
   `φ_J(mu_p = 0.3)` "`φ_m`". This campaign cannot fix that, because there is no
   `mu_p = 1.0` bidisperse block either.
2. **G3, the primary gate, is UNAFFECTED.** It asks whether the *frictionless* branch
   inherits the bidisperse non-monotonicity. That question does not reference `φ_m` at all.
3. **G4 is unaffected** for the same reason: it is a prediction about `Δφ_0`.
4. ⚠ **A new prediction, registered now.** The CV shift is measured to be friction-independent
   on the lognormal family: `Δφ_J(CV 0→0.5)` = **+0.05803 / +0.06207 / +0.05621** at
   `mu_p` = 0 / 0.3 / 1.0, a spread of 0.0059 against s.e. ≈ 0.003 each. If that carries to
   the bidisperse family, **G4's band should hold for the reason stated rather than by
   coincidence**, and a G4 miss becomes more informative than it was when registered.

⚠ **And it retracts a statistic in the companion page, `wc_endpoints_move_together`.** That
page tested the window's flatness as a RANGE over six rungs against twice a two-point
propagated s.e. On the per-rung PAIRED statistic — which cancels the shared seed draw and is
what a referee asks for — the `mu_p = 0.3` window is **not** flat:
`W(0.5) − W(0) = −0.00404 ± 0.00015`, **t = 27.70**. With the converged `mu_p = 1.0`
endpoint it *is* flat: `−0.00058 ± 0.00172` and `+0.00182 ± 0.00065`, t = 0.34 and 2.78.
Two errors were cancelling: the wrong endpoint created a deformation and the wrong yardstick
hid it. See that page's own §0 for the full correction.

---

## 4. What is NOT being run, and why

- **No raised iteration cap.** If cells fail G1 they are reported, not rescued.
- **No extra seeds.** The bound on G2 is a property of n = 3 and is registered as such;
  adding seeds after seeing the answer would be choosing the design from the result.
- **No `mu_p = 1.0` bidisperse block.** Whether `mu_p = 0.3` stands in for Wyart-Cates's
  frictional limit is a separate question, answerable from the lognormal `mu_p = 1.0` ladder
  already on disk, and is being done as analysis rather than as cells.
- **Nothing is re-run at `mu_p = 0.3`.** That block exists and this campaign must not
  perturb it; the new tag guarantees it cannot pool.

## 5. What would make this campaign worthless

- W1 fails → the shared-cell premise is wrong and every CV = 0 number in the week's work
  needs re-deriving before anything here is read.
- More than 2 of 72 cells cap → the ladder is not comparable to the other two families' and
  G3's argmax is being read off a censored curve.
- The two top rungs of `φ_0^bidi(CV)` within 1 s.e. → G3 returns UNDETERMINED, which is a
  real outcome and closes the question at "not decidable at three seeds".

# PRE-REGISTRATION — the constant-pressure production ladder

**Written 2026-09-12, after reading `PREREG_CONSTP_2026-09-12.md` §7 and before any
production cell launched.** The shakeout cleared the two deciding gates (`φ_ss` moves
0.0023 across a 5× gain range and 0.0032 across preparations 0.09 apart), so a production
run is allowed. Every bar below is fixed here. Grader: `grade_constp_prod.py` →
`out/constp_prod.txt`.

---

## 1. The one question

**`Δφ_c`(CV 0 → 0.5), measured by a route that takes no input from the isotropic `φ_J`.**

The project's headline is `Δφ_c` = **+0.049 to +0.055** on the common ladder, and its
standing caveat is that every estimator behind it (E2, E4, E4m) needs a ladder placement,
a fitting window or an isostatic contact target, and at `μ_p` = 0.3 all three routes end
up reading the isotropic `φ_J` (`PREREG_ROOT_PINNING_2026-09-09.md` §7,
`PREREG_ZISO_2026-09-12.md` §8). Constant-pressure shear needs none of them.

## 2. The grid — 18 cells, and it contains no `φ_J`

`family` lognormal, `μ_p` = 0.3, **N = 1500**, `p*` = **3e-04** (an isotropic ladder rung),
`K` = **0.15** (inside the usable window [0.1, 0.25] the shakeout measured),
`γ_max` = 3, `shear_steps` = 200000.

- `CV ∈ {0.0, 0.2, 0.5}` — the common ladder's own levels.
- `φ_start ∈ {0.55, 0.64}` — **the same two absolute values at every CV**, which is what
  makes this free of `φ_J`: no ladder is centred on anything. The pair brackets every
  expected `φ_c` (E4 reads 0.5745 / 0.5850 / 0.6299 across the three CVs), so at every CV
  one cell dilates to the critical state and one contracts to it.
- `seeds ∈ {1, 2, 3}`.

## 2a. AMENDMENT, 2026-09-12, written before the amended run launched

**`N` = 800 for this run, not 1500, and the reason is the harness rather than the
physics.** Measured today: a process spawned by a `wsl.exe` invocation does **not**
survive that invocation exiting — not under `nohup`, not under `setsid nohup`, and not
while a second WSL session is held open (tested with a bare `sleep 300`, which was gone
by the next call). The agent harness caps one call at 600 s. An `N` = 1500 cell at
200 000 steps costs ~900–1000 s, so **no `N` = 1500 cell can complete in this
environment at all**; the first attempt was killed at 600 s having produced nothing.

**First attempt at `N` = 600 (the shakeout size), and it was WRONG — the level carrying
the effect cannot run there.** All six CV = 0.5 cells refused, from `packing`'s own
periodic-collider guard: *"largest grain d_max = 3.868 needs a final cell edge ≥ 9.671 but
phi = 0.72 gives 9.045. Use N ≥ 734."* This is not a surprise the repo had to discover —
`phase8_criticalstate` sets `N_PROD = 1500` and says "N is bounded below by the periodic
collider at CV = 0.5" — and it should have been read before choosing 600. The records were
deleted rather than kept, because a set whose CV = 0.5 arm is empty cannot answer P1 and
keeping it invites a later pooling across N.

**`N` = 800 then refused CV 0.5 SEED 3** ("needs a final cell edge >= 10.745 but phi = 0.72 gives 10.520. Use N >= 853"), because `d_max` is a random DRAW and the guard's requirement is therefore seed-dependent: seeds 1 and 2 clear at 800 and seed 3 does not. A refusal that fires only at CV 0.5 **censors the axis under test** -- the 2026-09-01 lesson -- and 2 surviving seeds would have met the letter of §3 while quietly making the high-CV arm thinner than the others.

**The run therefore goes at `N` = 1000.** The guard's requirement grows like `d_max` (the maximum of N draws from a truncated lognormal, i.e. slowly) while the available cell edge grows like N^(1/3), so the margin improves with N: needed/available was 1.22 at N = 600 and 1.07 at N = 800. `truncate` is deliberately NOT tightened, because that would change the distribution being compared. ~370 s/cell expected, inside the 600 s ceiling. Everything else in §2 stands.

⚠ **If a seed still refuses at N = 1000, the censoring is REPORTED, not chased**: the refusal count per CV goes in the outcome and P1 is quoted on the surviving seeds with that stated, because the alternative is picking N by which value makes the arms match.

Consequences, stated rather than buried:

- **P1 is still answerable**, because it is a CONTRAST at fixed `N`: whatever finite size
  does to `φ_c`, it does to all three CV levels, and the common-ladder reference it is
  compared against is itself a contrast. ⚠ It is not *guaranteed* to cancel — a
  finite-size effect that depends on CV would not — so P1's verdict carries "at N = 800"
  and the `N` = 1500 repeat is the check.
- **P2 cannot be read** and is recorded as DEFERRED, not failed: it compares N = 600 with
  N = 1500 and this run is neither. ⚠ Its bar (0.005 against the shakeout's 0.5935 at
  N = 600) is now a 600 -> 800 comparison rather than 600 -> 1500, so the grader reports
  the number and calls the gate DEFERRED.
- **P3 compares against N = 1500 isotropic values**, so any offset it reports now mixes
  the protocol difference with a finite-size difference (N = 800 vs 1500). Its SIGN is still informative
  (the shakeout's +0.0102 at CV 0 is 14 isotropic seed scatters, far above plausible
  finite-size drift), but the magnitude is not quotable.
- **P4 and P5 are unaffected** — both are internal to this run.

**What unblocks `N` = 1500:** either `vmIdleTimeout = -1` in `%USERPROFILE%\.wslconfig`
plus a genuine service-style launch (an author decision, not an agent's), or running the
campaign from a terminal the author keeps open. Recorded in
`NEXT_STEPS_2026-09-13.md` §3.

## 3. Admission, fixed in advance

A cell is **refused** — printed, never dropped silently — unless all of:

1. **stayed jammed**: `z_win_min ≥ 2.0` and `p_win_min ≥ 0.1 p*` over the averaging
   window. This is the shakeout's post-hoc diagnostic promoted to an admission rule,
   because a cell can hold `⟨p⟩` on average while visiting an unjammed state (measured:
   `p_hat` to 1e-08 against a 3e-04 set point at K ≥ 0.5).
2. **servo held**: `|⟨p⟩/p* − 1| ≤ 0.05`. ⚠ The shakeout held on only 4 of 10, so this is
   expected to bite. If fewer than 2 seeds survive at any (CV, `φ_start`), that cell block
   is reported UNRESOLVED rather than averaged over one seed.
3. **steady** on both channels and **quasi-static** (`I ≤ 1e-2`).

Every number below is quoted on the admitted set, with the refused count beside it.

## 4. Predictions and bars

**P1 — 🔑 the contrast agrees with the common ladder.**
`Δφ_c`(CV 0 → 0.5) within **±0.015** of the common ladder's midpoint +0.0525 — the same
band `PREREG_COMMON_LADDER_2026-09-09.md` §3 used for H_real. *Predicted: agrees.*
A miss is a real disagreement between two structurally different routes, and it would put
the headline back in doubt rather than the new route.

**P2 — the CV = 0 level reproduces the shakeout across a factor 2.5 in N.**
`φ_ss`(CV 0) within **0.005** of the shakeout's 0.5935 (N = 600 → 1500). *Predicted:
within.* A larger move is a finite-size effect on this protocol and must be reported as
one, because the shakeout's S7 corroboration of E4m rests on `N` = 600.

**P3 — steady shear stays denser than isotropic compression, at every CV.**
`φ_ss` − `φ_iso`(p = 3e-04, same CV) **> 0** on 3 of 3 CVs. At CV 0 the shakeout measured
+0.0102. *Predicted: positive at all three.* A sign flip at high CV would be a genuinely
new statement about polydispersity and the critical state, so it is worth being able to
see.

**P4 — start-independence holds at every CV, now with seeds.**
`|φ_ss(0.55) − φ_ss(0.64)|` ≤ **0.005** at each CV, on the admitted cells. *Predicted:
holds.* ⚠ At CV 0.5 the 0.55 preparation starts ~0.08 below `φ_c`, the furthest travel in
the set, so this is where `γ_max` = 3 is most likely to be too little strain — and the
diagnosis is in the `φ` series, not in the gap alone.

**P5 — `α` from both branches, reported not gated.** Sign must be negative on every
admitted cell (it was 10/10 in the shakeout). The dilating/contracting ratio is
**reported per CV**; the shakeout saw up to 4.41× and pre-registered that a clean
disagreement is physical. With 3 seeds the question becomes answerable: *predicted, the
branch ratio exceeds the seed scatter of `α` within a branch* — i.e. the asymmetry is
real rather than noise. This is the first prediction here that could go either way.

## 5. What would make the run uninformative

If admission (§3) leaves fewer than 2 seeds at any CV on either preparation, the contrast
is quoted on whatever survives **with the per-cell table printed**, and P1 is reported as
UNRESOLVED rather than passed or failed. Ten cells decided nothing about `α` and 18 may
not either; saying so is the outcome.

## 6. OUTCOME

*(empty until the cells land)*

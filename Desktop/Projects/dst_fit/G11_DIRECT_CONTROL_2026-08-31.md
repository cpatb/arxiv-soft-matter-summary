# The strain clock is definition-dependent at the 41 % level

**2026-08-31.** `scripts/probe_g11_direct_control.py`, exit 0. This is
`PLAN_2026-08-29.md` item 5 — the control on G11's method, and the one step G11 has been
taking on faith.

---

## 1. What was unchecked

G11 measures the strain clock `gamma_0` by inverting a post-reversal viscosity transient
through Wyart–Cates to get the frictional fraction `f(gamma)`, and fitting an exponential
to that. `f` is never observed; it is inferred. And `eta` is strongly nonlinear in `f`
near jamming, so an inversion that is wrong in a *strain-dependent* way would move the
clock without leaving a residual to notice.

Lin et al.'s `Fig3` is the control. Those panels are **simulations** of the same reversal
protocol, so the contact and hydrodynamic stresses are evaluated directly rather than
inferred: `Fig_3b` gives both against strain at `gamma_dot = 0.74 gamma_dot_0`. So the
same record can be read two ways, and the clock fitted twice.

**The two channels differ only by the inversion.** Both come from `Fig_3b`; the same
`fit_recovery` code path fits both, with the same binning, the same exact amplitude
elimination and the same log-grid profile scan. The direct channel reaches that code by
being pushed *forward* through the Wyart–Cates map first, so the function inverts it back
to itself — verified to 1.5e-15 before anything is fitted.

---

## 2. The answer

| channel | `gamma_eff` | 68 % interval | chi2_r | bins |
|---|---|---|---|---|
| Wyart–Cates inversion of the total viscosity | **0.4881** | [0.4623, 0.5152] | 553.7 | 91 |
| directly evaluated contact stress fraction | **0.2892** | [0.2764, 0.3025] | 501.5 | 98 |

**A factor of 1.69, or +41 %, on the same record.** The intervals do not come close to
overlapping.

**And the difference is bit-identical across two packing triples** — 40.75 % and 40.75 %,
spread 0.00 % — which is the pre-registered check that it is the *inversion* and not the
gauge. That is the affine-invariance argument G11 relies on, working exactly as claimed:
the two triples are `(0.50, 0.58, 0.66)` and `(0.40, 0.52, 0.70)`, and they are far enough
apart that the inverted `f` leaves the admissible interval `[0, 1]` in **opposite
directions** — running to −1.546 under the first and to +1.266 under the second — and the
clock ratio does not move at all.

---

## 3. What it means, and what it does not

**It does not show the inversion is wrong.** Wyart–Cates's `f` is a fraction of frictional
*contacts*; `s = sigma_contact / sigma_total` is a fraction of *stress*. There is no law
saying they must be equal, or even affinely related. Two honest quantities can disagree.

**What it shows is that `gamma_0` is definition-dependent at the 41 % level**, on the one
record in the literature where both definitions can be evaluated. G11 quotes 0.54 to 1.64
on the inversion definition and does not name it as a definition. It should.

The mechanism is visible in the record rather than inferred:

| strain | `eta_r` | `s` (direct) | `f` (inverted, triple 1) | n |
|---|---|---|---|---|
| 0 – 0.5 | 8.23 | 0.3335 | **−1.3560** | 190 |
| 0.5 – 2 | 26.21 | 0.7615 | 0.3641 | 571 |
| 2 – 8.49 | 36.89 | 0.7907 | 0.7668 | 2475 |

**`s` is 96 % of the way to its plateau by strain 2; `f` is 47 %.** That gap is the whole
of the 41 %. The two are affine to R² = 0.880 with a worst residual of 0.2565 against an
`s` span of 0.6903 — 37 % of the range — so the map between them is substantially
nonlinear, and it is most nonlinear exactly where the clock is set: over the first half
strain unit, where the inversion is outside its own domain.

---

## 4. What follows

1. **G11's range needs the definition named.** "`gamma_0` = 0.54 to 1.64" becomes
   "`gamma_0` = 0.54 to 1.64 **read from the Wyart–Cates inversion**, and a directly
   evaluated contact fraction on a simulated record of the same protocol reads 41 %
   shorter."
2. **The ordering survives; the margin does not.** G11's headline is that every one of 21
   estimates exceeds the `gamma_0 = 0.2` the Basilisk solver runs (attribution corrected 2026-09-27: Richards et al. 2019 used gamma_0 = O(0.1); 0.197 is Han et al. 2018's start-up scale; see PRIOR_ART_STRAIN_CLOCK_2026-09-27.md), by 2.7× to 8.2×.
   Carrying a 41 % downward systematic takes the bottom of that from 2.7× to **1.6×**.
   Still above one, so the ordering claim stands and the solver's 0.2 is still low — but
   the comfortable end of the margin is gone, and 1.6× is close enough that it should be
   stated rather than rounded past.
3. **This is a control on the method, not a second measurement.** Lin's simulation rates
   are in units of its own `gamma_dot_0`, so the clock here is dimensionless. It cannot
   be compared with the experimental 0.54–1.64 as a number; only the *ratio between two
   readings of one record* transfers.

---

## 5. Two smaller things the probe turned up

**`Fig_3a` and `Fig_3b` differ by a constant factor of 1.024300, to five digits**
(standard deviation 1.3e-5 over 3236 shared strains). They are the same run plotted with
two different reductions. Nothing in the deposit's READMEs explains the factor and this
does not either; it is recorded because a 2.4 % offset between two panels of one figure is
the sort of thing that gets discovered later and blamed on the analysis.

**A pre-registered test was replaced because it measured the wrong thing.** P1 was first
written as "the two panels agree to better than 2 % in median" and it failed, at 2.43 %.
That was the wrong test rather than a failed one — a fixed ratio held to five digits is
what "the same run" looks like and is not what "different runs" looks like — so it was
replaced by a constancy test, which is the stricter bar (0.5 % on a ratio against 2 % on a
value). Both numbers are printed by the probe and the original wording is kept in its
docstring.

**A rank correlation is not reported and the reason is worth keeping.** 76 % of this
record is the post-recovery plateau, so Spearman is computed almost entirely on plateau
noise: it comes back at 0.04 where Pearson is 0.94. The same asymmetry is why
`fit_recovery` log-bins before fitting, and it is why the R² above is quoted as a shape
diagnostic and not as a goodness of fit.

---

## 6. Reproducing

```
python scripts/probe_g11_direct_control.py     # exit 0 if P1, P2 and P4 hold
```

P3 prints its verdict either way, deliberately: "the inversion distorts the clock" is as
much a result as "it does not", and gating the exit code on the outcome would be
pre-registering a conclusion.

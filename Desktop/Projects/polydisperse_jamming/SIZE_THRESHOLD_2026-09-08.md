# Is there a critical particle size? — the hypothesis, enumerated and narrowed

**2026-09-08.** Author's hypothesis, in his words: *"maybe there is just a size threshold
after which you transition from solid-like to inertial-dilatant thickening."*

Everything below is printed by `size_threshold.py` into `out/size_threshold.txt`; no
number here is typed by hand. Companion: `MECHANISM_2026-09-08.md` and
`derivations/dilatancy_sign_mechanism.md`.

---

## Result box

> **The hypothesis is not one hypothesis — it is seven**, because a size threshold is
> exactly what a critical value of a size-dependent dimensionless group *is*. Enumerated
> and filtered on two requirements (does it separate 800 nm silica from 10 µm
> cornstarch, and does it have the right *polarity*), **six die and one survives**:
>
> | | candidate | `a*` [µm] | verdict |
> |---|---|---|---|
> | **T1** | **macroscopic drainage, `ID = 6.6e-3`** | **7.66** | **SURVIVES** |
> | T1′ | macroscopic drainage, `ID = 1` | 94.3 | both materials below |
> | T2 | grain-scale drainage | — | predicts *no* threshold |
> | T3 | dilatant strengthening capacity | 1.3–4.0 | **reverse polarity** |
> | T4 | Brownian, `Pe = 1` | 0.113 | both materials above |
> | T6 | Stokes, `St = 1` | 61.2 | both below |
> | T8 | capillary / granulation | 143 | both below |
> | T7 | grain inertia, `I = 1` | 2121 | both below |
>
> And the sharpened version is better than the original: because the two live
> size-dependent effects grow at **different rates and point opposite ways** —
> strengthening capacity as `a¹`, drainage as `a²` — **the solid-like state is a WINDOW
> in particle size, not a half-line**, with a lower edge set by **salt** and an upper
> edge set by **velocity**.

---

## 1. Why "a threshold" has to be turned into "which group"

Any statement of the form "above size `a*` the behaviour changes" is a statement that
some dimensionless group crosses a critical value, and the group is what carries the
content: two thresholds can sit at the same size and still be told apart, because they
move differently when you turn a knob. So the useful question is not *is there a
threshold* but *which group sets it*, and the way to answer it is to enumerate every
group in the problem whose value depends on `a`, compute where each puts `a*`, and ask
which ones the observation can already distinguish.

Seven were found. Each formula is an identity, a group already implemented in this
repo's own modules, or a textbook definition, and every convention that has variants
(Péclet, Stokes) is stated at the point of use because the variants move `a*` by up to
~2×.

## 2. The two filters

**Filter 1 — separation.** The observation to be explained is that 800 nm silica is
solid-like and 10 µm cornstarch reaches the inertial-dilatant mode. A candidate must
therefore put `a*` *between* them. Four candidates fail here: T4 (0.113 µm — both
materials are already non-Brownian at impact strain rates), T6 (61 µm), T8 (143 µm) and
T7 (2121 µm) all leave both materials on the same side, so none of them can be what
distinguishes the two.

**Filter 2 — polarity, and this is the one that does real work.** A candidate must have
`a > a*` on the **inertial** side. **T3 does not.** Its criterion (`PORE_PHYSICS.md`
Eq. 9.1) is

```
   p_dil^max / P_conf = (a / 6λ_D)(Δφ/φ)
```

which *grows* with size, so above `a*_T3` frustrated dilation is *more* able to
strengthen the pack — bigger grains are **more solid-like**, the reverse of the
hypothesis. T3 is a real criterion and it is excluded as *this* threshold.

T2 is excluded by construction: if the drainage length scales with the grain, every
factor of `a` cancels out of `ID` and there is **no threshold at all**. That is a
falsifiable prediction rather than a missing number, and §5 says how to test it.

**One survivor: T1**, macroscopic drainage, at `a* = 7.66 µm`.

## 3. What the surviving candidate costs: the ID boundary

`a*` goes as `√(ID_boundary)`, and at `ID_boundary = 1` it is **94.3 µm** — which puts
cornstarch on the *solid* side and contradicts [M26]'s own observation. For the threshold
to fall between the two materials,

```
   ID_boundary ∈ [7.2e-5, 1.1e-2]
```

which brackets the **6.6e-3** that [M26]'s Eq. 6 returns at their own Regime-2/3
transition, and **excludes `ID = 1` by 89×**.

⚠ **This is an inference *from* the hypothesis, not evidence *for* it.** It says only
that the hypothesis is self-consistent on the smaller calibration and not on the nominal
one. But it does bear on the Status card's standing landmine (§8, "`ID` Eq. 6 does not
evaluate to 1 at the boundary it is drawn at, and a co-author can settle it"): the
size-threshold reading independently requires the smaller value.

## 4. 🔑 The sharpened hypothesis: a window, not a threshold

T1 and T3 are not redundant. Both grow with size and they grow at **different rates**:

| | criterion | scaling | meaning of "above" |
|---|---|---|---|
| **T3** | strengthening capacity | `a¹` | dilatancy *can* raise the particle pressure |
| **T1** | drainage | `a²` | the pore fluid keeps up, so it *cannot* |

Below T3 dilatancy has nothing to give; above T1 the fluid drains and it does not matter
what it had. **Frustrated dilation can hold a drop solid only in the band between them.**

| salt | `U₀` [m/s] | lower edge (T3) | upper edge (T1) | inside the window |
|---|---|---|---|---|
| 0 mM | 0.5 | 3.96 µm | 18.8 µm | **cornstarch** |
| 0 mM | 3.0 | 3.96 µm | 7.66 µm | *neither* |
| 1 mM | 0.5 | 1.32 µm | 18.8 µm | **cornstarch** |
| 1 mM | 3.0 | 1.32 µm | 7.66 µm | *neither* |
| 100 mM | 0.5 | 0.13 µm | 18.8 µm | **silica + cornstarch** |
| 100 mM | 3.0 | 0.13 µm | 7.66 µm | **silica** |

Three consequences, each on an axis the lab already turns:

**(a) [M26]'s re-entrant transition is recovered as a window edge sweeping past a fixed
particle size.** Cornstarch at 10 µm is *inside* the window at 0.5 m/s and *outside the
top* by 3 m/s, because the upper edge falls as `U₀^(-1/2)`. That is their velocity-driven
Regime 2 → Regime 3 transition, obtained here without fitting anything to it.

**(b) ⚠ 800 nm silica sits BELOW the window at 0 and 1 mM salt — so its solid-like branch,
which is observed, is probably not frustrated dilation at all.** It would have to be
frictional / repulsive (Wyart–Cates) jamming. If that holds, **silica and cornstarch are
solid for different reasons**, which matters well beyond this question: it is a caution
against reading silica as a small-particle model of cornstarch. Note the honest limit of
this claim — the criterion says the poroelastic channel is *weak* there
(`p_dil^max/P_conf` = 0.20 at 0 mM), not that it is zero.

**(c) The lower edge crosses 800 nm inside the lab's own salt series.** `a*_T3 = 6λ_D φ/Δφ`
falls with screening: 3.96 µm at nominal zero salt, 1.32 µm at 1 mM, 0.40 µm at 10 mM,
0.13 µm at 100 mM. **So the same silica should enter the dilatancy window somewhere
between 1 and 10 mM** — no new particles, no new apparatus, and Paper 1 already sweeps
this axis for a different reason.

## 5. What separates the candidates — the experiment

This is the payoff of enumerating rather than picking. The thresholds sit at different
sizes, but more usefully they **move differently**:

| candidate | vs `U₀` | vs salt | vs `Δφ` | vs `D₀` |
|---|---|---|---|---|
| T1 drainage (`ℓ = D₀`) | `U₀^(-1/2)` | none | `Δφ^(+1/2)` | `D₀^(+1/2)` |
| T2 drainage (`ℓ ∝ a`) | *no threshold* | — | — | — |
| T3 strengthening | **none** | **`∝ λ_D`** | `Δφ^(-1)` | none |
| T4 Brownian | `U₀^(-1/3)` | none | none | `D₀^(+1/3)` |
| T6 Stokes | `U₀^(-1/2)` | none | none | `D₀^(+1/2)` |
| T8 capillary | `U₀^(-2)` | none | none | none |

> **T3 is the only candidate that moves with salt and the only one that does not move
> with velocity.** Every other finite threshold falls with `U₀` and is salt-blind.

So the discriminating experiment is a **size ladder crossed with a velocity ladder and a
salt ladder**, and the readings are unambiguous:

- `a*` falls as `U₀^(-1/2)` and ignores salt ⇒ **T1**, drainage.
- `a*` ignores `U₀` and tracks `λ_D` ⇒ **T3**, strengthening capacity.
- **no threshold anywhere on the size axis** ⇒ **T2**, grain-scale drainage — which
  would also settle the drainage-length question that `derivations/…` §6 calls the
  largest single uncertainty in the lab's use of `ID`.
- `Δφ` is a second, independent separator with the **opposite sign** in the two live
  candidates: raising `φ` *raises* `a*_T1` (7.66 → 16.5 µm from φ = 0.44 to 0.48) and
  *lowers* `a*_T3` (3.96 → 1.44 µm). A φ ladder alone therefore widens the window from
  one end and narrows it from the other.

## 6. Does a pure size threshold explain the **CV** observation?

On a pure size threshold, width acts only by moving the moment the threshold reads. T1
reads Carman–Kozeny, whose length is the **Sauter** mean, and for a lognormal at fixed
number mean `d₃₂/d₁₀ = (1+CV²)²`. So crossing `a*` from a given `d₁₀` requires a specific
CV — a number, not a story:

| `d₁₀` [µm] | `d₃₂/d₁₀` needed to reach 7.66 µm | CV needed |
|---|---|---|
| 0.8 | 9.58 | **1.45** |
| 1.5 | 5.11 | **1.12** |
| 3.0 | 2.55 | **0.77** |

**This is the sharpest thing on the page for the bench.** At the lab's 800 nm number-mean
a pure size threshold needs **CV ≈ 1.45** to carry a sample across — far broader than the
CV 0.2–0.5 this project has simulated, and broader than most Stöber batches.

> **If the author's two samples have CV below about 1, a pure size threshold acting
> through `d₃₂` cannot carry them across T1 and something else must be moving as well.
> If they are broader than that, it can.**

Which makes the **raw size distributions of the two bench samples** — already Q1 on the
Status card, still open — the measurement that decides between the author's hypothesis
and the mechanism in `MECHANISM_2026-09-08.md`. Not a new experiment: a measurement of
samples that already exist.

## 7. What would kill this

| observation | verdict |
|---|---|
| `a*` measured, and it does **not** move with `U₀` | T1 dead; T3 or something unlisted |
| `a*` measured, and it moves with **salt** | T1 dead |
| no threshold found anywhere over a decade of `a` | T1 dead, T2 alive, and `ID`'s `a²` scaling goes with it |
| the two bench samples have CV ≪ 1 and still switch regime | a *pure* size threshold cannot be the whole story |
| silica shows a pore **suction** signature at 0 mM salt | §4(b) wrong — frustrated dilation is acting below where the criterion says it can |

## 8. Referee

*(pending — `physics-referee` pass in flight at the time of writing. Nothing in §§1–7 is
to be quoted into the Status card before it returns.)*

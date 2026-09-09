# Is there a critical particle size? — the hypothesis, and what survives a referee

**2026-09-08.** Author's hypothesis, in his words: *"maybe there is just a size threshold
after which you transition from solid-like to inertial-dilatant thickening."*

Everything below is printed by `size_threshold.py` into `out/size_threshold.txt`.
Companions: `MECHANISM_2026-09-08.md`, `derivations/dilatancy_sign_mechanism.md`.

> ## 🔴 THE FIRST VERSION OF THIS PAGE WAS WRONG AND IS WITHDRAWN
>
> It claimed "of seven candidate thresholds exactly one survives", and that headline was
> an artefact. It computed a single critical size `a*` **at the lab silica's operating
> point** and then sorted **cornstarch** with it — but `a*` depends on `φ_c`, `D₀` and
> `ρ`, and [M26]'s cornstarch differs in all three. Recomputed at cornstarch's own point,
> `a*` = **15.4 µm** and cornstarch falls **below** its own threshold, so the "survivor"
> does not separate the two materials at all.
>
> This is the project's own standing lesson — *never grade a knob at one setting of the
> knob it interacts with* — committed again, one day after it was written down. Found by
> `physics-referee`, confirmed independently before acting.
>
> Also withdrawn: the ID-boundary "inference" (§3, an identity), and the size **window**
> as a quantitative claim (§4, thinner than its own error bars).

---

## 1. What replaces it

A size threshold is what a critical value of a size-dependent dimensionless group *is*,
so the hypothesis is a family with one member per such group. The repair is to stop
converting to a critical size across materials, because that conversion smuggles in one
material's operating point. Two statements, of which only the second involves an `a*`:

**(A) Cross-material — compare the GROUPS, not a size.** Evaluate every size-dependent
group at each material's own parameters:

| group | silica @ 3 m/s | starch @ 0.5 | starch @ 7.2 | separated? |
|---|---|---|---|---|
| `ID` drainage, `ℓ = D₀` | 7.19e-05 | 4.64e-04 | 6.67e-03 | **yes** |
| `ID` drainage, `ℓ = 10a` | 7.03 | 0.742 | 10.7 | no |
| `p_dil^max/P_conf` capacity | 0.202 | — (Hertzian) | — | n/a |
| `Pe` Brownian | 358 | 7.3e4 | 1.0e6 | yes |
| `St` | 1.71e-04 | 2.15e-03 | 0.031 | yes |
| `ρU₀²/p_capillary` | 5.6e-03 | 1.7e-03 | 0.343 | no |
| `I` inertial number | 3.8e-04 | 2.8e-03 | 2.8e-03 | yes |
| **`De = τ_p/τ_imp`** | 3.6e+05 | — (Hertzian) | — | n/a |
| **`τ_χ/τ_imp`** | **1.06** | 1.5e-03 | 0.021 | yes |
| **`Pe_grav`** | 0.260 | 3489 | 3489 | yes |
| **`Re_p`** | 1.1e-03 | 0.015 | 0.222 | yes |
| **`λ_D/r_throat`** | **0.484** | 0.039 | 0.039 | yes |

The last five were **missing from v1's enumeration** and were added on the referee's
list. Two of them matter:

- **`τ_χ/τ_imp` = 1.06 for the lab's silica** — Jérôme's own validity condition for
  treating the pore fluid as incompressible, and the lab's particle sits *exactly on it*
  while cornstarch is two decades clear. `inertial_dilatant.py` calls this "the
  size-dependent statement that survives review", and v1 did not enumerate it.
- **`λ_D/r_throat` = 0.484 at nominal zero salt** — the double layer fills half the pore
  throat, which is why §5 below retracts the claim that the drainage groups are
  salt-blind.

**(B) Within one material — `a*` is well posed**, because a ladder that holds chemistry,
`φ`, `φ_c`, `D₀` and `U₀` fixed and moves only `a` has no other operating point to
smuggle. **That is exactly the experiment being proposed (E-4)**, so the concept survives
precisely where it is going to be used, and nowhere else.

## 2. What can still be said about silica vs cornstarch — without any `a*`

[M26] observe the regime transition somewhere inside their own velocity range, so the
boundary lies between the `ID` their cornstarch has at the bottom and the top of it:

```
   ID_boundary ∈ [4.64e-04, 6.67e-03]
   lab silica  =  7.19e-05
```

**Silica sits 6.4× below even the bottom of that interval**, at the lab's `φ`, `D₀` and
3 m/s.

> On the drainage reading, the lab's silica **cannot reach the inertial regime by
> velocity alone** at the lab's drop size. It would need a larger particle or a larger
> drop.

That is the defensible version of the author's hypothesis, and it needs no critical size.

⚠ **The v1 claim that the hypothesis "independently requires" an ID boundary near 6.6e-3
is withdrawn as circular.** The bracket was algebraically `ID(silica)` to
`ID(cornstarch)` at the lab's parameters, and 6.6e-3 is itself just `ID(cornstarch)` at
the top of [M26]'s velocity range — computed here as **6.675e-03**. Saying "the bracket
contains 6.6e-3" restates where 6.6e-3 came from.

## 3. The window: structure survives, numbers do not

Two criteria grow with size at different rates and point opposite ways — capacity as
`a¹`, drainage as `a²` — so *in structure* frustrated dilation can hold a drop solid only
in a band. **That structure is sound.** The numbers are not:

- v1 quoted the band as [3.96, 7.66] µm, a factor **1.94**.
- `PORE_PHYSICS.md` §6.3 states the affine `K_p` is an **upper bound** and real packings
  measure **2–5× softer**. `K_p` sets the lower edge, so the lower edge is a *lower
  bound* and moves up by the same factor: **at 2× the window is already empty.**
- The upper edge is calibrated (through the ID boundary) and the lower edge is not, so
  the *ratio* carries the whole O(10)–O(100) prefactor freedom `inertial_dilatant.py`
  declares in its own docstring.

**A 1.94× band cannot be asserted against a ≥2× one-sided uncertainty on one of its
edges.** The window is withdrawn as a quantitative claim and with it the sorting of any
material by it — including v1's §4(b) conclusion that silica's solid branch "has to be"
frictional jamming, which was an elimination from a list that never enumerated the
frictional-onset group.

⚠ A separate error in the same section: v1 computed the lower edge for **cornstarch**
from the double-layer contact law. Cornstarch is **Hertzian** (`PORE_PHYSICS.md` §6.4),
where `K_p ∝ E^(2/3)P^(1/3)` carries **no particle size at all**, so that criterion has no
lower edge for cornstarch and both cornstarch rows of v1's window table were computed
with the wrong law. (This project had already written the caveat down — "`a¹` belongs to
the double layer" — and then violated it two sections later.)

## 4. The salt discriminator is broken

v1's headline experimental design rested on "the capacity criterion is the **only**
candidate that moves with salt". **False.** `PORE_PHYSICS.md` §5.2 — cited by v1 — says
the double layer occupies the pore throat at low salt, and at 800 nm:

| salt | `λ_D` | `λ_D/r_throat` |
|---|---|---|
| 0 mM | 30 nm | **0.484** |
| 1 mM | 10 nm | 0.161 |
| 10 mM | 3 nm | 0.048 |
| 100 mM | 1 nm | 0.016 |

So the permeability that the drainage groups are built on is itself salt-dependent at the
small end, by an order-unity factor. A salt ladder moves the drainage groups too — and,
separately, moves `σ*`, i.e. the frictional mechanism. **The salt axis is not a clean
discriminator.** What survives is that the candidates still scale differently with
**velocity** (`U₀^(-1/2)`, `U₀^(-1/3)`, `U₀^(-2)`, none), so a size × velocity ladder
retains discriminating power.

## 5. The CV question — the one part that got stronger

On a pure size threshold, width acts by moving the moment the threshold reads.
Carman–Kozeny's length is the **Sauter** mean and `d₃₂/d₁₀ = (1+CV²)²` for a lognormal at
fixed number mean, so crossing `a*` needs a specific CV. **The algebra was checked and is
exact, and the Sauter mean is the right moment** (CK's length is *defined* by
`S_v = 6/d₃₂`).

But v1 held `φ_c` fixed while doing it, and this project has now **measured** that `φ_c`
rises with CV (`MECHANISM_2026-09-08.md` §8: +0.058 to +0.062 over CV 0 → 0.5). Since
`a* ∝ √(φ − φ_c)`, broadening moves the *threshold down* as well as moving `d₃₂` up, and
the two push the same way:

| `φ_c` | `Δφ` | `a*` (silica ladder) | CV needed via `d₃₂` |
|---|---|---|---|
| 0.420 | 0.0200 | 7.71 µm | **1.45** (v1's number) |
| 0.425 | 0.0150 | 6.65 µm | 1.37 |
| 0.430 | 0.0100 | 5.40 µm | 1.26 |
| 0.435 | 0.0050 | 3.80 µm | **1.09** |

*(`a*` here targets `ID` = 6.675e-03, cornstarch's own value at the top of [M26]'s
velocity range, rather than the rounded 6.6e-3 — which is why it reads 7.71 and not the
7.66 v1 quoted.)*

A `φ_c` shift of +0.015 — a quarter of what the DEM measured at CV = 0.5 — takes the
requirement from 1.45 to ~1.1. **So v1's bold "CV ≈ 1.45 needed" is an upper bound that
assumed a fixed `φ_c`; the honest range is ~1.0–1.5**, and the two CV routes are not
independent.

⚠ **And the provenance of "800 nm" is unstated.** Whether it is a *number* mean is not on
record; `psd.py` notes laser diffraction reports `d₄₃`. If 800 nm is a volume or intensity
mean then `d₃₂` **falls** as CV rises and this route reverses sign. That is Q1 on the
Status card, and it remains the measurement that decides the most here.

## 6. What would kill what is left

| observation | verdict |
|---|---|
| a critical size found in a single-material ladder, and it does **not** move with `U₀` | the drainage reading is dead |
| no threshold found anywhere over a decade of `a` | drainage dead, grain-scale drainage alive, and `ID`'s `a²` scaling goes with it |
| the two bench samples have CV ≪ 1 and still switch regime | a *pure* size threshold is not the whole story |
| `τ_χ/τ_imp` ≈ 1 turns out to matter | the incompressible-pore-fluid assumption fails for the lab's own particle, and every group above is affected |

## 7. Referee

`physics-referee`, second pass, on this page and its script. **Three FATAL and five
MATERIAL findings, all acted on above**: the operating-point artefact (§0), the window's
error bars (§3), the incomplete enumeration (§1), the circular boundary bracket (§2), the
salt discriminator (§4), the fixed-`φ_c` CV number and the `d₄₃` provenance (§5), the
Hertzian/EDL mix-up (§3), and the fact that the two-filter procedure could not return
"no threshold" because the null was eliminated by the first filter by construction.

Findings recorded as **surviving** the attack: the capacity criterion's inversion and its
reverse polarity; radius-versus-diameter consistency throughout; the Péclet convention;
the Sauter-mean algebra in §5. Minor and not acted on: two polarity labels
(`St`, capillary) are arguable but do not change any verdict, and the candidate numbering
skips T5 with no note.

⚠ The referee could not resolve the window construction to a published DOI and marks it
`[UNVERIFIED]` as new; it is in kind the standard two-sided bound (capacity below,
drainage above) rather than a new result.

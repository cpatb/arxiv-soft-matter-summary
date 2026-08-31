# The airflow axis, and what it does to the containment story

**2026-08-31.** Code in `cfd/`, results in `cfd/cfd_provisioning.json`,
`cfd/cfd_floor_probe.json`, `cfd/cfd_provisioning_fit.json`.

`DESIGN.md` §8 has carried one open gap since the Tier-1 CFD landed:

> *Over-provisioning carries no thermal benefit.* `under = max(0, 1-p)` means `B` is
> identical at `p = 1` and `p = 1.4`, so blowing extra cold air costs `p^3` fan power
> and buys nothing. A two-sided `phi(p)` is the fix and needs the 3D run to calibrate.

It did not need the 3D run. It needed the 2D solver to have a through-flow, which it
did not have, because the cavity was buoyancy-only and therefore had no `p` in it at
all. It has one now, and the answer to the gap arrived with two things attached that
matter more than the gap did.

---

## 1. What was added

Cold supply enters through a floor grille under the cold aisle and hot air leaves
through a ceiling grille over the hot aisle. In a streamfunction formulation a
through-flow is not a new equation, it is a non-zero boundary `psi`: `psi` *is* the
cumulative volumetric flux along the boundary, so the openings are the places where it
changes and everything else is a constant. It enters the existing sparse factorization
as a Dirichlet lift, and the Thom wall vorticity becomes a difference rather than a
value. The partition stands on the floor to the right of the inlet, so its constant is
`psi = 0` and the existing internal-solid handling is already right.

`flux = 0` reproduces the buoyancy-only cell exactly: the three stored points in
`cfd_recirc.json` come back to 1.1e-5, which is the JSON's own rounding to four
decimals. The eta-only calibration is untouched by this change.

**`p` is anchored on the cavity, not on a fan.** `Q_1 = max|psi|` of the buoyancy-only
solution at `eta = 0`, and `p = flux / Q_1`. So `p = 1` means "the CRAH pushes as much
air through the aisle as the unconfined buoyant loop circulates on its own". That is a
definition, chosen at `eta = 0` so the two axes stay independent. It is the weakest
link in everything below and §5 sweeps it.

---

## 2. The sweep: 48 cells, 0 diverged

`phi`, the recirculated fraction at the rack-inlet face, Ra = 4e3, N = 40:

| eta \ p | 0.00 | 0.25 | 0.50 | 1.00 | 1.50 | 2.00 | 3.00 | 4.00 |
|---|---|---|---|---|---|---|---|---|
| **0.0** | 0.2411 | 0.1216 | 0.0862 | 0.0435 | 0.0272 | 0.0226 | 0.0211 | 0.0213 |
| 0.2 | 0.1897 | 0.0772 | 0.0496 | 0.0179 | 0.0063 | 0.0028 | 0.0011 | 0.0006 |
| 0.4 | 0.1382 | 0.0582 | 0.0332 | 0.0091 | 0.0022 | 0.0006 | 0.0001 | 0.0000 |
| 0.6 | 0.1163 | 0.0530 | 0.0272 | 0.0060 | 0.0012 | 0.0003 | 0.0000 | 0.0000 |
| 0.8 | 0.1007 | 0.0385 | 0.0147 | 0.0020 | 0.0003 | 0.0001 | 0.0000 | 0.0000 |
| **0.95** | 0.0729 | 0.0144 | 0.0030 | 0.0002 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

**Over-provisioning buys a lot of thermal benefit.** At `eta = 0`, `phi` falls 5.5x
between `p = 0` and `p = 1` and another 1.9x by `p = 2`, then saturates. The model
credited that whole range with nothing.

A three-parameter dilution law was fitted and **rejected**, which is worth recording
rather than hiding: `phi = [phi0 exp(-eta/tau) - floor] / (1 + p/k) + floor` puts
`floor` on its own lower bound and misses the `eta = 0.95, p = 4` cell by four orders
of magnitude. The floor it is trying to absorb is `eta`-dependent — at `eta = 0` there
is no partition and the hot wall conducts straight to the inlet face, giving a large-p
limit of 0.0213, while at `eta = 0.95` the same limit is 0. One additive constant
cannot be both. `dc_thermal` interpolates the measured grid instead.

---

## 3. The bigger finding: the "7% recirculation floor" is conduction

`CFD_RECIRC` gives `phi(eta = 1) = 0.0713`, and `DESIGN.md` builds "containment ROI
flattens; ~7% floor" on it — the finding it calls the real technical content of the
CFD step. The sweep put a number beside it that the eta-only sweep never had to face:
at `eta = 0.95, p = 0` the temperature channel reads **0.0729** and the passive tracer
reads **0.0002**, a factor of 360.

The two channels differ in one thing: the tracer's diffusivity is 0.02 against the
thermal 1, i.e. Lewis 50. So sweeping it separates advection from conduction, because
**an advected signal is Lewis-independent and a diffused one is proportional to the
diffusivity**. At Ra = 1e4, N = 40 (`cfd/run_floor_probe.py`):

| tracer diffusivity | Le | phi from T | phi from tracer | ratio |
|---|---|---|---|---|
| 1.00 | 1 | 0.0694 | 0.06945 | **1.000** |
| 0.50 | 2 | 0.0694 | 0.05717 | 0.823 |
| 0.20 | 5 | 0.0694 | 0.03001 | 0.432 |
| 0.02 | 50 | 0.0694 | 0.00253 | 0.036 |

and the control at `eta = 0`, where the signal should be advective:

| tracer diffusivity | Le | phi from T | phi from tracer | ratio |
|---|---|---|---|---|
| 1.00 | 1 | 0.2256 | 0.22563 | 1.000 |
| 0.02 | 50 | 0.2256 | 0.22644 | **1.004** |

At perfect containment the signal scales with the diffusivity all the way to Le = 1,
where the passive scalar reproduces the temperature to five decimals. At zero
containment the same 50x sweep moves it by 0.4%. **The floor is diffusion.**

That is a statement about this cell, which runs at Ra = 4e3 to 1e5. A hall is
Ra ~ 1e9 to 1e12, where molecular conduction across a cold aisle is nothing. So the
floor does not transfer, and neither does the flattening of the containment ROI that
rests on it.

**A pre-registered test failed to decide this and it is recorded as failed.** The
first version of `run_floor_probe.py` fixed the rule "conduction if `phi_T` at
eta = 0.95 falls by more than 2x from Ra = 4e3 to 1e5, advection if it changes by less
than 30%". It fell by 1.42x and changed by **29.8%** — inside the advection branch by
two parts in a thousand. A rule that lands on its own boundary has not decided
anything, and the Lewis sweep above is what replaced it. The Ra ladder is kept in the
JSON because the monotone fall of `phi_T` with Ra (0.0729, 0.0694, 0.0599, 0.0512 at
N = 40; 0.0655, 0.0623, 0.0537, 0.0455 at N = 60) is consistent with the verdict, but
it is not what establishes it.

---

## 4. What it does to the model

`dc_thermal` gains `recirc_model = "cfd2d_p"`, interpolating the measured grid. The
mass-balance floor `phi >= 1 - p` stays, because the cavity has no rack fan and
nothing in the CFD knows that a hall short of supply air must make the difference up
from the hot aisle. The `(1 + under)` amplification is **not** applied on top of the
table — the table already measures the p dependence, and applying both would count
under-provisioning twice.

| | free optimum `(T_s, p)` | PUE | total kW | containment ROI at `p = 1` |
|---|---|---|---|---|
| `linear` placeholder | 21.02 C, 0.685 | 1.2351 | 249.00 | 0.0214 |
| `cfd` (eta only) | 21.02 C, 0.685 | 1.2351 | 249.00 | 0.0126 |
| **`cfd2d_p`** | **21.02 C, 0.685** | **1.2351** | **249.00** | **0.0022** |

Two results, and the first is the one the gap was actually about.

**Closing gap (1) moves the free optimum by nothing.** All three models give the same
optimum to five digits, because it sits at `p = 0.685`, below design, where the
mass-balance floor `1 - p = 0.315` dominates every geometric term and `eta` drops out
of `B` entirely. So the one-sidedness of the old `p` axis was never what pinned the
optimizer — mass balance and `p^3` fan power were. And the new two-sided branch does
not tempt it upward either: at `p = 2` the recirculated fraction is half what it is at
`p = 1`, and eight times the fan power is not worth it.

**The containment ROI shrinks by 5.7x**, 0.0126 to 0.0022, and the supply temperature
containment unlocks at `eta = 0` rises from 23.0 C to 26.4 C — because at design
airflow the cold aisle is already flooded and there is little recirculation left for
containment to remove.

---

## 5. The anchor is the weak link, so it was swept

`P_ANCHOR_SCALE` moves the identification between the hall's provisioning ratio and
the cavity's `Q_1`. Containment ROI (PUE at `eta = 0` minus PUE at `eta = 1`):

| p \ scale | x0.25 | x0.5 | x1 | x2 | x4 |
|---|---|---|---|---|---|
| 0.60 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 0.80 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 1.00 | 0.0063 | 0.0046 | 0.0022 | 0.0011 | 0.0011 |
| 1.50 | 0.0054 | 0.0033 | 0.0014 | 0.0010 | 0.0011 |
| 2.00 | 0.0046 | 0.0022 | 0.0011 | 0.0011 | 0.0011 |

Three things survive a 16x move in the anchor:

1. **The free optimum does not move at all** — `T_s = 21.02 C`, `p = 0.685`,
   PUE 1.2351, identical to five digits at every scale. It is on the mass-balance
   floor, where the calibration is irrelevant.
2. **The containment ROI is smaller than the 0.0126 currently quoted under every
   anchor**, by between 2x and 11x.
3. **Below design airflow the containment ROI is exactly zero**, at every anchor,
   because `1 - p` pins `B` and containment has nothing to act on. That is a
   structural property of the model, not of the CFD — and the model's own optimizer
   wants to run there.

---

## 6. Limitations, stated rather than buried

- **Ra = 4e3 for the production sweep**, 4e3 to 1e5 for the floor probe. A hall is
  1e9-1e12. The `eta` and `p` shapes are what is being carried across, not the
  absolute conduction share, and §3 is the reason to say so.
- **Laminar and two-dimensional.** No turbulence model, no rack-to-rack coupling, no
  third dimension. The off-diagonal structure of `B` is still the geometric
  `exp(-D/coupling_len)` placeholder and no CFD has touched it.
- **No rack fan.** The rack is a solid partition, so the cavity cannot know what
  fraction of the supply the servers actually ingest. That is why the mass-balance
  floor has to be imposed on top of the table rather than measured in it, and why the
  `p` anchor cannot be read off the geometry.
- **The grid is clamped, not extrapolated**, above `p = 4` and outside
  `eta` in [0, 0.95]. Above `p = 4` the table is already flat, so holding it is the
  measurement; the `eta` clamp is an assumption on a 5% sliver.
- **`cfd_recirc.json`'s numbers are unchanged and still correct** for what they are:
  the `p = 0` edge of this surface. What changed is the claim that they describe a
  hall at design airflow.

---

## 7. Next

1. **The `eta` axis should be re-measured at high Ra.** §3 disqualifies the floor but
   does not replace it; the advective `phi(eta)` at hall-like Ra is unmeasured, and
   the containment ROI needs it.
2. **Anchor `p` on something physical.** A rack with a fan — a prescribed
   through-partition flow — would give the cavity an airflow demand to normalise
   against, and would turn `P_ANCHOR_SCALE` from a sweep into a measurement.
3. The 3D `B_ij` off-diagonals remain untouched and remain the honest gap in
   `DESIGN.md` §8.

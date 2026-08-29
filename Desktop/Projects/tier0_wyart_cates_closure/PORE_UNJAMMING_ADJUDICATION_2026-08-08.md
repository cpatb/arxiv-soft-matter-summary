# What actually sets the un-jamming clock

**2026-08-08/09.** Adjudication between two independent investigations that were run
without sight of each other and reached incompatible conclusions. One derived the
un-jamming timescale from pore-scale physics; the other went and found what has been
measured. They disagree, and the disagreement is the useful part.

Sources: `PORE_UNJAMMING_DERIVATION_A.md` (this folder) and
`Desktop\Papers\STATE_OF_THE_FIELD\UNJAMMING_TIMESCALE_LITERATURE_2026-08-08.md`.

---

## 1. Why this was worth a night

`wc_drop_impact.c` evolves a frictional-fraction field `fhat` toward its
stress-determined steady value on one of two time constants — `tau_up` while rising,
`tau_relax` while falling. **Both are placeholders in every preset**, and the solver's
own header says why: steady flow curves fix `{eta_s, phi0, phiM, sigma*, delta}` and
nothing else, because a transient cannot be read off a steady curve. The values in the
code — 2.0 s for silica, 0.02 s for cornstarch — were guessed.

They are also the parameters that decide whether the drop is jammed while it matters. So
two hand-set numbers were gating the result the whole DNS exists to produce.

The steer that opened this was: *the theoretical pore size of a jammed suspension affects
how it relaxes out of a jam, and how long it stays jammed.*

---

## 2. What the derivation concluded

To un-jam, the frictional network must dilate. The drop is incompressible and its free
surface is capillary-sealed — the Laplace pressure across a 148 nm pore is ~1 MPa against
a `sigma*` of 44 Pa, a factor of 2×10⁴ — so the new pore volume has to be fed by liquid
dragged through the packing. That Darcy flow is then the clock:

```
tau_relax = eps_dil * eta_f * L_d^2 / ( k(phi,a) * sigma* )
k(phi,a)  = (1-phi)^3 a^2 / (45 phi^2)        Carman-Kozeny
r_h       = (1-phi) a / (3 phi)               hydraulic pore radius
```

It ranked the alternatives and disposed of them cleanly. Brownian relaxation scales as
`a^3`, which *inverts* the observed ordering. Inertial and acoustic times are 5–8 orders
too fast. Lubrication pull-apart at a single contact gives `3*pi*eta_f/(2*sigma*)`, in
which **the particle radius cancels exactly** — 1.2 ms for silica, 3.7 ms for cornstarch.

That cancellation is the derivation's sharpest observation, and it generalises: *every*
contact-scale mechanism reduces to `eta_f/sigma*` times a number that does not contain
`a`, and `eta_f/sigma*` is three times **larger** for cornstarch. So every local mechanism
predicts the wrong sign. Only a mechanism carrying a length that does not scale with the
particle can produce a slow silica and a fast cornstarch — and permeability, `k ∝ a^2`,
is exactly such a mechanism.

The numbers were the persuasive part:

| | derived | placeholder in the solver |
|---|---|---|
| silica | 1.5 s (0.3–5) | 2.0 s |
| cornstarch | 0.013 s (0.004–0.04) | 0.02 s |
| **ratio** | **111** (40–300) | **100** |

The ratio contains no free parameter — it is `(eta_s/eta_c)(sigma_c*/sigma_s*)(k_c/k_s) =
1.205 × 0.257 × 360`, and it reads 111 / 113 / 107 under Carman–Kozeny, Happel and
Beetstra, so it is not an artifact of one permeability model. Better still, fixing the one
free length `L_d = 55.5 um` to reproduce silica's 2.0 s then **predicts cornstarch at
0.0180 s against a placeholder of 0.020 s** — one parameter, one datum, the second datum
predicted to 10% across a factor of 100.

That is a genuinely good piece of physics, and it is why it needed a second opinion rather
than a nod.

---

## 3. What the literature says instead

**The un-jamming time has been measured, repeatedly, and nobody needed to derive it.**

| study | system | measured release time |
|---|---|---|
| Maharjan & Brown, PRF 2, 123301 (2017) | cornstarch, phi_eff 55–61% | 0.02 – ~1 s, `T1 ∝ gammadot_c^-1`, exponent 1.0 ± 0.1 |
| Cho et al., PRF 7, 063302 (2022) | cornstarch 13 um | two-step, `tau_1, tau_2 ≈ 1e-2 – 1e1 s` |
| Richards et al., PRL 123, 038004 (2019) | cornstarch, phi = 0.47–0.48 | **`t_r = 0.24 ± 0.05 s`** |
| Barik & Majumdar, PRL 128, 258002 (2022) | PS 2.65 um in PEG | `tau_f ≈ 0.2 s`, `tau_s ≈ 2 s` |

Three things follow, and each one damages the derivation somewhere different.

**(a) The cornstarch placeholder is wrong, and it is wrong in the direction that
manufactured the agreement.** Every measurement puts cornstarch release at 0.2–2 s. The
solver uses 0.02 s — at or below the bottom of the measured range, ~10× too fast. The
derivation predicted 0.013 s and was scored against 0.02 s. **It agreed with a guess.**
Once the target moves to the measured 0.2–2 s, silica and cornstarch are no longer two
decades apart, and a mechanism whose whole job was to produce two decades of separation is
solving a problem that does not exist.

**(b) The size scaling is not there.** Measured release times cluster at 0.1–2 s across a
5× range of particle diameter, from 2.65 um polystyrene to ~13 um cornstarch. `a^-2`
predicts a 25× separation over that range. Nothing like it appears.

**Honest caveat, because this is the load-bearing point:** the measured times *do* span
about 20×, which is comparable to the 25× being tested for. What is missing is any
*correlation* with size — the scatter runs across studies, materials and protocols, and
**no single study has varied particle size systematically while holding everything else
fixed.** So the evidence against `a^-2` is strong enough to stop us building it into a
solver, and not strong enough to call it refuted. See §6.

**(c) Darcy is in the problem, but on the other axis.** Mobaseri, Driscoll & Cheng
(arXiv:2601.12642) build the same Kozeny–Carman permeability into a drop-impact model — so
the pore-scale picture is not exotic, it is current practice in this exact geometry. But
their dilatancy development time evaluates to **~9 us for 10 um cornstarch**, five decades
below the measured release times. Athani et al. (PRF 10, 043301) reach the same structure
independently: Darcy bounds the *stress*, while the state variable relaxes on a *strain*
clock.

> **The resolution: pore-scale drainage sets the stress, not the clock.**

The steer was pointing at real physics and at the wrong observable. Permeability governs
how large the dilatant stress is; it does not govern how long the jam lasts.

---

## 4. What the clock actually is

While the material is flowing, the state variable relaxes on a **strain**, not a time.
Richards et al. write it as one rate covering both directions:

```
dfhat/dt = -( gammadot/gamma0 + 1/t_r ) [ fhat - fss(sigma) ]
```

`gamma0 ≈ 0.2` — and Han et al. (PRF 3, 073301) measure the same quantity independently at
**0.197 ± 0.002**, which is about as good as agreement gets in this field. `t_r = 0.24 s`
takes over once flow stops.

This is strictly better than what the solver had, for a reason that has nothing to do with
which number is bigger: **the loading/unloading asymmetry stops being an input.** The old
model produced asymmetry by typing in `tau_up != tau_relax`. Here the same rate applies
both ways, and the asymmetry falls out — loading happens while flowing, so the strain term
dominates and it is fast; unloading happens largely at rest, so only `t_r` survives. The
measured asymmetry is a factor ~2000 (`t_c/t_r = 5e-4` in Richards' own sample), and the
law reproduces it without a second parameter.

Implemented as `WC_CLOCK=1` in `wc_drop_impact.c`, commit `9fbad8a`, four gates passing —
including an exact-limit gate in which `gamma0 -> infinity` provably reduces the new law to
the old one and reproduces it byte for byte.

**And it is not a cosmetic change.** Impact shears the drop at `gammadot ~ U0/D0 ~ 670 s^-1`,
so the strain term gives an effective clock of ~0.3 ms against the 2 s the silica preset
assumed — a factor of ~7000. The old law held `fhat` **frozen through the entire impact**.
At L8, peak `jam_frac` moves from 0.116 to 1.000 when the law is swapped.

⚠ **That headline number is not yet trustworthy and is under test.** Both new values sit at
their ceiling, and the solver caps the relaxation rate at `1/(4*Delta)` = 8 at L8 against a
requested strain term of 5–50 — so L8 may have measured the grid rather than the law. A
resolution ladder plus a `gamma0`-sensitivity row (`wcclock_ladder.sh`) is running to
separate them. The *ordering* — the measured law jams the drop far more than the guessed
one — is unlikely to move; the value 1.000 may well be numerical.

---

## 5. On electrostatics

The steer was that electrostatics is not a big parameter. **Half right, and the half that
is wrong is in the experiment rather than the solver.**

*For the solver:* right, and for a stronger reason than "the effect is small". Whatever
electrostatics does, it does through the onset stress `sigma*` — a single scalar the solver
already takes as a calibrated input. There is no separate double-layer physics to model.

*For the material:* not small, and specifically not small for the Paper 1 system.
Ishibashi et al. (Colloids Surf. A 732, 139148, 2026) fit silica 110–2580 nm as
`sigma_c = 5kT a^-3 + 0.01 F* a^-2` with `F* = 0.3 nN` at 0.1 mM and **`F* = 0` at 100 mM** —
i.e. the double-layer term is the whole story at low salt and gone at high salt, above a
critical coagulation concentration. Maranzano & Wagner (J. Rheol. 45, 1205) get a smaller
but real factor 4.4 (750 -> 170 Pa) on charge neutralisation.

Against that, chemistry and roughness are larger where they act: Bourrianne et al.
(PRR 4, 033062) take *identical* 300 nm fumed silica, change only the surface chemistry, and
**extinguish DST entirely**; roughness alone moves `phi_DST` from 18% to 8.8%.

So: minor for cornstarch-in-water and for screened silica above the CCCS; **dominant for
sub-micron charge-stabilised silica at low ionic strength**, which is exactly the
salt-screening axis Paper 1 is built on. Treating it as a small parameter there would
discard the experiment's own independent variable.

---

## 6. What would settle the open question

The size scaling is the one thing genuinely unresolved, and it is unresolved because nobody
has done the obvious experiment.

1. **Stress-drop rheometry at fixed material, varying gap, 100 um to 2 mm.** A drainage
   mechanism predicts `tau ∝ gap^2` — a 400× swing. A contact mechanism predicts none.
   Nothing else separates them so cleanly, and it needs no new material. **2–4 days.**
2. **One material, three particle sizes, same chemistry, same phi/phi_M, same protocol,
   stress-drop release time measured for each.** This is the measurement the entire
   literature is missing, and the reason the `a^-2` question cannot be closed today. A
   null result here is publishable on its own — it would rule out the drainage picture for
   the release clock across the field, not just for us.
3. Cheap and immediate: **run the DNS at `WC_CLOCK=1` with `t_r` swept over the measured
   0.1–2 s bracket** and report the sensitivity, rather than quoting one value. If the
   impact-phase results are insensitive across that whole range — which the strain term
   makes likely, since it dominates during impact — then the clock stops being a threat to
   the programme's conclusions and can be reported as such.

---

## 7. Verdict

- **`tau_relax` is no longer a free parameter, and the replacement is measured, not
  derived.** `gamma0 = 0.197 ± 0.002` and `t_r = 0.24 ± 0.05 s` come from two independent
  experiments; the two numbers they replace came from neither.
- **The pore-size steer identified real physics on the wrong axis.** Permeability sets the
  magnitude of the dilatant stress; it does not set the release time. Keep the mechanism,
  move it to the stress.
- **The derivation's agreement with the placeholders was an artifact of the placeholders
  being wrong**, and it is worth recording how convincing it looked: no free parameter in
  the ratio, robust across three permeability models, one-parameter cross-prediction good
  to 10% across a factor of 100. The check that caught it was not a better derivation, it
  was going and looking up what had been measured.
- **The size scaling is open**, and closing it is a two-to-four-day bench job that no one in
  the field appears to have done.

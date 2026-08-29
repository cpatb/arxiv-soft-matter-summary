# Stress in an impacting drop, and the inertial–dilatant number for 800 nm silica

**Prepared for the Driscoll group, 6 August 2026.**
All numbers regenerate from `results.py`; all equations are checked in `test_impact_stress.py` (40 tests).
Working directory: `C:\Users\compl\Desktop\Projects\impact_stress\`

---

## 0. The short version

Your PI asked two things. Here are the answers, then the reasoning.

**(1) Stress localization and stress values in an impacting drop.** There is a complete, closed-form
answer in the literature, and one of its two key papers has a figure taken in your own lab. The wall
pressure is

$$p(r,t) \;=\; \frac{3}{\pi}\,\frac{\rho U_0^2 R}{\sqrt{d(t)^2 - r^2}}, \qquad d(t)=\sqrt{3RU_0 t}$$

This is a self-similar field whose maximum sits **in a ring at the moving contact line, not at the impact
point**. The wall shear stress has the same structure but a stronger, $1/(\text{distance to the rim})$
divergence, and is smaller than the pressure by a factor $\sim\!\mathrm{Re}^{-1/2}$. For your drops
that factor is about 10, and it is the difference between a gauge that works and one that does not.

**(2) The inertial–dilatant number.** The term exists, it is one paper old, and **Michelle is a
co-author**: Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng, *Phys. Rev. Lett.* **136**,
148201 (2026). They define

$$\mathrm{ID} \equiv \frac{F_i}{F_d} = \frac{\rho\,\kappa\,U_0}{A\,\mu_w\,\Delta\phi\,D_0}$$

For your 800 nm silica, ID comes out around $6\times10^{-5}$, roughly **220× smaller than for the
cornstarch the paper was built on**, because $\mathrm{ID}\propto\kappa\propto a^2$ and your particles
are 12× smaller. Colloidal silica sits far into the dilatant (solid-like) corner, and the re-entrant
liquid regime that is the paper's headline is out of reach on a bench drop tower.

Three things came out of doing this carefully that are worth your attention:

- **The gauge your project has been using, $\rho U_0^2$, cannot measure an onset, and the shear
  gauge can.** $\rho U_0^2$ is 1.9–36 kPa over your velocity range against $\sigma^\ast \sim 20$–100 Pa,
  so every drop of every suspension is 40–800× above threshold. The wall shear stress is ~10× smaller
  and crosses $\sigma^\ast$ **inside your measured velocity window**. Your project already retired the
  inertial gauge on exactly these grounds; this supplies the replacement.
- **Equation 6 of the ID paper does not evaluate to 1 at the boundary it is drawn at.** Two independent
  calculations put it at $6.6\times10^{-3}$ at the paper's own Regime-2/3 conditions. Capping the
  dilatant pore pressure at what the pore fluid can physically sustain (~1 atm of suction) moves it to
  0.35 at the same point, within a factor of 3 of unity instead of a factor of 100. That cap is a
  proposed extension, it is falsifiable, and you can simply ask a co-author whether the Supplemental
  Material carries a prefactor that arXiv does not show.
- **Your fitted $\dot\gamma_c$ is convention-bound by a factor of 4.5.** The $U_0/D_0$ convention your
  front collapse uses is the one Cheng's group refuted in PNAS last year. This does not make your
  number wrong; it makes it a number that has to be quoted with its convention attached.

---

## 1. Orientation: an impact is four events, not one

Take one of your drops: 800 nm silica at $\phi = 0.44$ in water, $D_0 = 2.5$ mm, $U_0 = 2.3$ m/s
(velocity band B), onto a glass slide.

| quantity | symbol | value |
|---|---|---|
| suspension density | $\rho = \phi\rho_p + (1-\phi)\rho_f$ | 1440 kg/m³ |
| low-stress viscosity | $\eta_0$ | 34 mPa·s (bracket 7–34, see §10) |
| sound speed (Wood's law) | $c$ | 1629 m/s |
| Reynolds number | $\rho U_0 D_0/\eta_0$ | 241 |
| Weber number | $\rho U_0^2 D_0/\gamma$ | 265 |
| Ohnesorge | $\eta/\sqrt{\rho\gamma D_0}$ | 0.067 |
| Mach number | $U_0/c$ | $1.4\times10^{-3}$ |
| gas Stokes number | $\eta_{\rm air}/(\rho U_0 D_0)$ | $2.2\times10^{-6}$ |
| inertial time | $R/U_0$ | 544 µs |
| impact time | $D_0/U_0$ | 1087 µs |

Now the four events, in order. Every number below is for that drop.

### Stage 0: the air cushion (before contact)
The drop never touches the glass cleanly. Air has to be squeezed out of the closing gap, and it cannot
get out fast enough, so it lifts the drop's underside into a dimple. Scales:
film thickness $h_d \sim R\,\mathrm{St}^{2/3} \approx$ **210 nm**, dimple radius
$\ell \sim R\,\mathrm{St}^{1/3} \approx$ **16 µm**, gas pressure
$p_g \sim \rho U_0^2\,\mathrm{St}^{-1/3} \approx$ **0.6 MPa**. The important consequence is not the
number: it is that the gas film **screens the water hammer** and delays true contact.

### Stage 1: the water hammer (nanoseconds)
At first touch the liquid cannot get out of the way, so it is compressed. One-dimensional acoustics
gives $p = \rho c U_0 =$ **5.4 MPa**, the largest stress anywhere in this problem and 700× the dynamic
pressure. And it is almost entirely irrelevant, because of where and when it acts:

- it lasts $t_d = RU_0/2c^2 =$ **0.54 ns**
- over a contact patch of radius $RU_0/c =$ **1.8 µm**
- touching a fraction $\tfrac{3}{8}(U_0/c)^3 = 1.1\times10^{-9}$ of the drop's volume.

One part in a billion. Say this out loud in a talk: *the biggest stress in a drop impact is a billion
times too small to matter, and the air film screens most of what's left.*

### Stage 2: Wagner (the one that matters)
Once the shock detaches, the liquid is effectively incompressible, and the flow is set by an
expanding wetted disc. This is where the force, the pressure field, and the stress localization live.
Derived in §2.

### Stage 3: the lamella
The drop is now a spreading sheet. Pressure has collapsed; what survives is **viscous shear at the
wall**. This is the stage that matters for a shear-thickening drop, and it is §5.

---

## 2. The pressure field, from scratch

You should be able to reproduce this on a whiteboard. It takes four lines.

**Step 1. How big is the contact patch?** A sphere of radius $R$ descending at $U_0$ intersects the
plane in a circle of radius $\sqrt{2RU_0t}$. But the free surface piles up ahead of the geometric
intersection, so the *wetted* radius is larger. Wagner's condition fixes the correction:

$$\boxed{\;d(t) = \sqrt{3\,R\,U_0\,t}\;}$$

Cheng, Sun & Gordillo write this as $r_0 = \sqrt{6D V t}/2$; with $D=2R$ it is the same equation.
It has been confirmed experimentally by Riboux & Gordillo and by Gordillo, Sun & Cheng.

**Step 2. Replace the drop by a disc.** Near the wall, the flow is that of a flat disc of radius
$d(t)$ pushed into a half-space of liquid at speed $U_0$. The velocity potential on such a disc is
classical:

$$\Phi(r,t) = -\frac{2}{\pi}U_0\sqrt{d(t)^2-r^2}$$

**Step 3. Unsteady Bernoulli.** In an impact, the pressure comes from the *acceleration* term, not
the kinetic term:

$$p = -\rho\frac{\partial\Phi}{\partial t} = \frac{2}{\pi}\rho U_0\,\frac{d\,\dot d}{\sqrt{d^2-r^2}}$$

**Step 4. The trick.** Because $d^2 = 3RU_0t$ is linear in $t$, the product $d\dot d$ is a *constant*:

$$d\dot d = \tfrac12\frac{\mathrm d(d^2)}{\mathrm dt} = \tfrac32 RU_0$$

Substituting,

$$\boxed{\;p(r,t) = \frac{3}{\pi}\frac{\rho U_0^2R}{\sqrt{d^2-r^2}}
= \underbrace{\frac{\sqrt3}{\pi}\rho U_0^2}_{\text{scale}}\;
\underbrace{t^{*-1/2}}_{\text{amplitude}}\;
\underbrace{(1-\xi^2)^{-1/2}}_{\text{shape}}\;}$$

with $t^*=U_0t/R$ and $\xi = r/d(t)$.

**This is exactly Eq. [10] of Cheng, Sun & Gordillo (2022), prefactor included.** The test suite
checks it algebraically, and Philippi, Lagrée & Antkowiak (2016) get there independently, obtaining
the same closed form $P(\xi,0)=3/(\pi\sqrt{3-\xi^2})$ from the full self-similar solution.

### How to read it

Three factors, three separate ideas.

- $\rho U_0^2$ sets the stress scale.
- $t^{*-1/2}$ is the amplitude. **The earlier you look, the bigger it is.** It equals $\rho U_0^2$ at
  $t^*=3/\pi^2 = 0.30$, i.e. about a third of the way through the impact.
- $(1-\xi^2)^{-1/2}$ is **the shape, and this is the stress localization.** Flat across the middle of
  the contact patch, divergent at its rim.

The maximum is therefore **off-axis**, in an annulus racing outward with the contact line. That is the
counter-intuitive bit and it is worth a slide: a steady jet hitting a wall has its maximum pressure at
the stagnation point, because for a steady flow $\partial\Phi/\partial t=0$ and only the (negative)
kinetic term survives. An impact is the opposite: the acceleration term dominates, and it is largest
where the contact line is sweeping fastest.

The ring is not a theoretical curiosity. It has been seen three ways: as ring cracks in zinc sulphide
from 340 m/s drop impacts; as an off-axis pressure peak on piezo sensors under a 2D drop; and directly,
by Cheng's group, imaging the deformation of an elastic gel under an impacting drop
(Sun *et al.*, *Nat. Commun.* **13**, 1703, 2022).

![**Figure 1.** The self-similar wall stress. Left: pressure, shape $(1-\xi^2)^{-1/2}$, flat across the patch and divergent at the rim. Right: shear stress, shape $\xi/(1-\xi^2)$, exactly zero on the axis and more sharply peaked at the rim; the dashed line is $\sigma^*$ for the silica. Each curve is one instant. The shapes are identical and only the amplitude changes.](figures/fig1_selfsimilar_fields.png)

### Where the compressible stage plugs back in

Equation above blows up as $t^*\to0$. What stops it? Set the Wagner pressure equal to the water-hammer
pressure:

$$\frac{\sqrt3}{\pi}\rho U_0^2\,t^{*-1/2} = \rho c U_0 \;\;\Longrightarrow\;\; t^*_\times = \frac{3}{\pi^2}\mathrm{Ma}^2 = 0.30\,\mathrm{Ma}^2$$

Independently, the geometric shock-detachment argument gives $t^*_d = \tfrac12\mathrm{Ma}^2$. Two
arguments sharing no algebra, agreeing to a factor of 1.6. So the compressible stage is not a separate
phenomenon bolted onto $t=0$. **It is what makes the incompressible solution finite.** That is a nice
line for a talk.

---

## 3. Stress localization, quantified

Everything at $U_0 = 2.3$ m/s, evaluated at the peak-force time.

| stage | stress | ÷ $\rho U_0^2$ | localized over | lasts |
|---|---:|---:|---:|---:|
| gas cushion | 0.59 MPa | 77 | 16 µm | 0.09 µs |
| water hammer (rigid) | 5.4 MPa | 708 | 1.8 µm | 0.0005 µs |
| Wagner pressure, on axis | 7.5 kPa | 0.98 | 1.2 mm | 171 µs |
| **Wagner pressure, at the rim** | **42 kPa** | **5.5** | **20 µm** | 171 µs |
| wall shear, area-averaged | 1.4 kPa | 0.18 | 1.2 mm | 171 µs |
| wall shear, at the rim | 13 kPa | 1.8 | 20 µm | 171 µs |
| — reference: $\rho U_0^2$ | 7.6 kPa | 1 | — | — |
| — reference: $\sigma^\ast$ (silica) | **44 Pa** | 0.006 | 0.4 µm | — |

Peak force $F_{\max} = 0.79\,\rho D_0^2U_0^2 =$ **37 mN** at $t=171$ µs. That is measurable with a
piezo sensor, which is exactly what the 2026 paper did and something your rig does not currently do.

![**Figure 2.** Every stress scale of one impact on one axis, with $\sigma^*$ in red for comparison. The water hammer is the largest and the least relevant.](figures/fig2_stress_ladder.png)

**The localization factor.** The pressure at the rim runs 5.5× the on-axis value, at all times. That
number comes from cutting the singularity off at the screening length
$\Delta = 0.02\sqrt{VDt}$ that Philippi *et al.* measured in their simulations, i.e. at 1.6% of the
contact radius. **Treat 5.5 as a soft number.** Philippi *et al.* say plainly that the physical origin
of $\Delta$ is unknown, and offer two candidates that give different answers, one bounded and one not.
The area-averaged shear stress does not have this problem, because its integral converges; that is why
§5 uses it as the gauge.

---

## 4. What changes on a fluid or soft substrate

Your bath work and the Carbopol ladder both live here. There are three separate laws, and the mistake
to avoid is applying the wrong one.

### (a) The compressible stage: acoustic impedance sharing
Two half-spaces meeting at relative speed $U_0$ share the load in proportion to their acoustic
impedance $Z = \rho c$:

$$p = U_0\,\frac{Z_1Z_2}{Z_1+Z_2} \qquad\Longrightarrow\qquad \frac{p}{p_{\rm rigid}} = \frac{Z_2}{Z_1+Z_2}$$

Rigid target ($Z_2\to\infty$) recovers $\rho c U_0$. A bath of the *same* liquid gives **exactly one
half**. This is de Haller's 1933 result, still used in the droplet-erosion literature.

### (b) The incompressible stage: inertia sharing
Redo the Wagner derivation with a target that also moves. Both sides develop the same interfacial
pressure through their own potential flow, so $\rho_1u_1 = \rho_2u_2$ with $u_1+u_2=U_0$, giving

$$\frac{p}{p_{\rm rigid}} = \frac{\rho_2}{\rho_1+\rho_2}$$

**The same law with $Z\to\rho$.** Two media sharing an interface share the load in proportion to their
inertia: acoustic impedance when waves matter, mass density when they do not. Matched liquids: one
half again.

This is not a new result in disguise; it is the same statement as two things already in print.
Bisighini *et al.* report that the drop/target interface penetrates a deep pool at $U_0/2$, and Langley,
Castrejón-Pita & Thoroddsen unify soft- and rigid-substrate air-film scalings by defining an
*effective impact velocity* $V_e = V/2$ for a deformable surface. The contribution here is putting
(a) and (b) side by side and noticing they are one law.

*Caveat, stated so it is not lost:* the factor 3 in $d^2 = 3RU_0t$ comes from free-surface pile-up on
the **drop** side. With a deformable target there is pile-up on both sides, so that 3 becomes some
other order-one number. Equation (b) is the leading-order amplitude change; the geometric prefactor is
an open $O(1)$.

### (c) An elastic solid: a deflection criterion, not a sharing law
A solid does not flow; it deflects, by $w\sim pd/E^\ast$. That matters once $w$ is comparable to the
gap the drop is closing at the same radius, $\sim d^2/2R$. With the Wagner pressure this reduces to a
single group:

$$\mathcal{C} = \frac{\rho U_0^2}{E^\ast} \qquad\text{(the Cauchy number)}$$

$\mathcal{C}\gtrsim 1$: soft for the whole impact. $\mathcal{C}\ll1$: rigid except for the first
$\mathcal{C}$-th of it. This is the same group Yokoyama *et al.* (*Nat. Commun.* **17**, 607, 2026) use
to organise the crossover from inertial to Hertzian peak-force scaling, and the same comparison
Howland *et al.* (*PRL* **117**, 184502, 2016) make when they show splashing is suppressed once the
peak Wagner pressure exceeds $E$.

### (d) A yield-stress bath: $J$ is the switch between (b) and (c)
Your Carbopol arm is already indexed on $J = \tau_y/\rho U_0^2$. In this framework $J$ is not a
stiffness. It is **the switch that decides which of the two laws applies**. $J\ll1$: the bath yields
everywhere the impact touches, behaves as a liquid, use (b). $J\gtrsim1$: it never yields, behaves as
an elastic solid, use (c). That is a cleaner statement of what the Carbopol ladder scans than
"stiffness", and it is worth having in the substrate-grant framing.

### Numbers for your targets, at $U_0 = 2.3$ m/s

| target | which law | reduction vs rigid | $\mathcal{C}$ | $J$ |
|---|---|---:|---:|---:|
| glass slide | rigid | 1.00 | $2\times10^{-6}$ | — |
| deep water bath | (b) inertia | **0.41** | — | — |
| deep bath of the same suspension | (b) inertia | **0.50** | — | — |
| Carbopol, $\tau_y=1$ Pa | (b), $J\ll1$ | 0.41 | — | $1.3\times10^{-4}$ |
| Carbopol, $\tau_y=30$ Pa | (b), $J\ll1$ | 0.41 | — | $3.9\times10^{-3}$ |
| PDMS gel, $E=100$ kPa | (c) elastic | — | **0.057** | — |
| acrylic, $E=3$ GPa | (c) elastic | — | $1.9\times10^{-6}$ | — |

Both your Carbopol concentrations sit at $J\lesssim4\times10^{-3}$, deep in the liquid-like corner.
The bath yields under the impact throughout the planned ladder. If you want to reach $J\sim1$ you need
$\tau_y \sim \rho U_0^2 \approx$ 7.6 kPa, which no Carbopol concentration in your design reaches.
That is worth knowing *before* the 20-run grid, not after.

**One caution the literature is loud about and worth repeating.** Compliance changes the *magnitude*
of the stress but, in three independent measurements, **not the localization length**: Howland finds
$a=\sqrt{3RUt}$ independent of $E$; Sun *et al.* find the turning point and the shear maximum
independent of $E$ and of wettability; Langley finds no clear stiffness dependence below a critical
value. So "a softer substrate spreads the load out" is not what happens. What happens is that the load
gets smaller and the ejecta slow down.

---

## 5. The suspension question: which stress does dilatancy see?

### The threshold
A dense suspension thickens when the local stress exceeds an onset stress $\sigma^\ast$, at which the
repulsion holding particles apart is overcome and contacts turn frictional (Wyart & Cates 2014). For
your material, four independent estimates:

| source | $\sigma^\ast$ (Pa) |
|---|---:|
| Ishibashi (2026), fitted law with a measured $F^\ast = 0.3$ nN | 19 |
| Lin *et al.*, *PRL* **115**, 228304 (2015), scaled by $a^{-2}$ | 63 |
| Shah, Arora & Driscoll (2022), **your own 830 nm particles** | ~100 |
| your repo's flow-curve fit, `params_fitted.yaml` | 44 |

Geometric mean of the three external ones: **55 Pa**. The repo's own 44 Pa sits inside their range,
which is the main reason to trust it as an order-of-magnitude anchor despite the fact that the fit
failed its own acceptance gate.

### The gauge problem, and why $\rho U_0^2$ fails
Your project already established this and it is worth restating precisely, because the fix follows from
the diagnosis. Over your accessible velocities:

$$\rho U_0^2 = 1.9\ \text{to}\ 36\ \text{kPa} \qquad\text{against}\qquad \sigma^\ast \approx 20\text{–}100\ \text{Pa}$$

That is $\Lambda_P = \rho U_0^2/\sigma^\ast = 43$ to 814. **Every drop of every suspension is far above
threshold at the slowest speed anyone can drop from.** A gauge with no reachable threshold cannot
measure a threshold.

### The fix, and why it is physics rather than fitting
Dilatancy is a *shear* phenomenon: $\sigma^\ast$ in Wyart–Cates is a shear stress, and isotropic
pressure does not by itself build the anisotropic force chains that thicken a suspension. So use the
shear stress the impact actually imposes. Philippi *et al.* solved for it, and Cheng, Sun & Gordillo
print it as their Eq. [13]:

$$\tau(r,t) = \frac{2}{\pi^{3/2}}\,\rho U_0^2\,\mathrm{Re}_R^{-1/2}\,t^{*-1/2}\,\frac{\xi}{1-\xi^2},
\qquad \mathrm{Re}_R = \frac{\rho U_0 R}{\eta}$$

Three things to notice, and they carry the whole argument:

1. **The prefactor is $\rho U_0^2\,\mathrm{Re}_R^{-1/2}$.** The shear stress is smaller than the pressure
   by $\sim1/\sqrt{\mathrm{Re}}$. At your $\mathrm{Re}_R\approx120$ that is a factor of 11, and once
   the geometric factor is included, the area-averaged shear stress is ~100× below $\rho U_0^2$.
2. **The shape is $\xi/(1-\xi^2)$**, a $1/(\text{distance to rim})$ divergence, sharper than the
   pressure's inverse square root, and **exactly zero on the axis** by symmetry. So the shear maximum
   stays at the drop's edge for the whole impact, while the pressure maximum migrates to the centre.
   (Sun *et al.* found the same $\xi/(1-\xi^2)$ form collapsed their gel measurements best.)
3. **$\tau\propto\sqrt\eta$, not $\eta$.** A factor-5 uncertainty in your drop's viscosity is only a
   factor-2.2 uncertainty in the stress. That is what makes this usable despite $\eta_0$ being unresolved.

For a robust single number, integrate rather than pick a point. Philippi's own drag force divided by
the contact area gives

$$\langle\tau\rangle = \frac{10.7}{3\pi}\,\rho U_0^2\,\mathrm{Re}_R^{-1/2}\,t^{*-1/2}
= 1.135\,t^{*-1/2}\sqrt{\frac{\rho\,\eta\,U_0^3}{R}}$$

That is **the proposed gauge**, and its ratio to $\sigma^\ast$ is the group that answers "did the impact
trigger thickening?":

$$\boxed{\;\Lambda_F = \frac{\langle\tau\rangle}{\sigma^\ast}
= 1.135\,t^{*-1/2}\,\frac{1}{\sigma^\ast}\sqrt{\frac{\rho\eta U_0^3}{R}}\;}$$

### Does it work?

| $U_0$ (m/s) | $\rho U_0^2$ (Pa) | $\Lambda_P$ | $\langle\tau\rangle$ (Pa) | $\Lambda_F$ | frictional area fraction |
|---:|---:|---:|---:|---:|---:|
| 1.15 | 1904 | 43 | 278 | 6.3 | 0.82 |
| 2.3 | 7618 | 172 | 788 | 17.8 | 0.97 |
| 5.0 | 36000 | 814 | 2524 | 57 | 1.00 |

Setting $\Lambda_F = 1$ gives an onset velocity

$$U_c = \left[\frac{\sigma^{\ast2}R\,t^*}{1.29\,\rho\,\eta}\right]^{1/3} \;\propto\; \sigma^{\ast\,2/3}$$

Across the full $\sigma^\ast$ (19–100 Pa), $\eta_0$ (7–34 mPa·s) and criterion brackets, $U_c$ lands
between **0.19 and 1.75 m/s**. Your measured onsets are **1.473 / 1.120 / 1.138 / 0.455 m/s** at
0 / 10 / 50 / 100 mM LiCl.

So: the shear gauge puts the onset in the right place, and the inertial gauge cannot put it anywhere.
**But be honest about what that is worth.** A prediction spanning a factor of 9 that brackets a
measurement spanning a factor of 3 is weak evidence on its own. Of the factor-9 spread, $\sigma^\ast$
contributes 2.9×, the criterion 1.8×, and $\eta_0$ 1.7×. See Figure 3.

![**Figure 3.** The inertial gauge $\rho U_0^2$ (grey) never approaches $\sigma^*$ over the accessible velocity range. The shear gauge $\langle\tau\rangle$ (blue, three viscosity choices) crosses it right where the jamming onset is measured (green lines).](figures/fig3_gauge_comparison.png)

### The direction that is actually useful: invert it
The *scaling* $U_c\propto\sigma^{\ast2/3}$ carries no fitted constant, so ratios between salt levels
are free of every uncertainty above, because viscosity, drop size, criterion and prefactor all cancel:

$$\frac{\sigma^\ast(\text{salt})}{\sigma^\ast(0\ \mathrm{mM})} = \left[\frac{u_0^\ast(\text{salt})}{u_0^\ast(0\ \mathrm{mM})}\right]^{3/2}$$

| salt | $u_0^\ast$ measured (m/s) | **predicted $\sigma^\ast/\sigma^\ast_{0\rm mM}$** | $\sigma^\ast$ implied (Pa) |
|---|---:|---:|---:|
| 0 mM | 1.473 | 1.00 | 77–404 |
| 10 mM | 1.120 | **0.66** | 51–268 |
| 50 mM | 1.138 | **0.68** | 53–274 |
| 100 mM | 0.455 | **0.17** | 13–69 |

**These four ratios are the falsifiable output of this whole analysis.** They need no rheometer
calibration, no absolute stress, and no assumption about your $\phi$ axis. Measure $\sigma^\ast$ on a
rheometer at the four salt levels and check them. Two features stand out and are testable: 10 mM and
50 mM are predicted to be indistinguishable (0.66 vs 0.68), and the drop to 100 mM is a factor of four.
A screening curve that plateaus between 10 and 50 mM and then falls is a specific, checkable shape.

The absolute values are consistent too, in a way worth noticing: the 100 mM number (13–69 Pa) sits
squarely inside the published bracket for this material (19–100 Pa), while the unsalted number
(77–404 Pa) sits above it. That is the right direction, because your nominal "0 mM" sample has a residual
ionic strength of only ~$2.5\times10^{-6}$ M from dissolved CO₂, so its double layer is unusually
unscreened and its $\sigma^\ast$ should be unusually high.

### Where in the drop is it frictional?
Because the shear field is known everywhere, you can locate the $\tau = \sigma^\ast$ contour. Solving
$\xi/(1-\xi^2) = K$:

$$\xi_c = \frac{-1+\sqrt{1+4K^2}}{2K}, \qquad K = \frac{\sigma^\ast\sqrt{t^*}}{(2/\pi^{3/2})\rho U_0^2\mathrm{Re}_R^{-1/2}}$$

Material at $\xi > \xi_c$ is above onset; material at $\xi < \xi_c$ is still lubricated. **So the
thickened region is an annulus hugging the rim, not a disc at the centre** (Figure 4). That is a
qualitative prediction your under-view footage can check directly, and as far as I can establish
nobody has computed this contour for a suspension drop before. Shah, Arora & Driscoll say in words
that the solidification front must be where "the stress falls below the critical stress"; this makes
that sentence into a calculable curve.

![**Figure 4.** Where in the contact patch the suspension is above $\sigma^*$ (shaded). The frictional region is an annulus at the rim that fills inward as the impact proceeds, and it appears earlier and reaches further in at higher speed.](figures/fig4_frictional_annulus.png)

---

## 6. The inertial–dilatant number

### What it is
From Mobaseri *et al.*, *PRL* **136**, 148201 (2026), Eq. 6. They set up two competing forces:

- **Inertial** (Wagner / Gordillo–Sun–Cheng): $F_i \sim \rho U_0^2D_0^2$
- **Dilatant** (Reynolds–Darcy, after Jérôme, Vandenberghe & Forterre): the suspension must expand to
  shear, expansion sucks solvent inward through the pores, and Darcy resistance to that flow generates
  a pore suction that presses the particles together. $F_d\sim A\mu_w U_0D_0^3\Delta\phi/\kappa$

Their ratio:

$$\mathrm{ID} \equiv \frac{F_i}{F_d} = \frac{\rho\kappa U_0}{A\mu_w\Delta\phi D_0} = \frac{t_d}{t_i}$$

with $\kappa = (1-\phi)^3a^2/180\phi^2$ (Carman–Kozeny, $a$ = particle **diameter**),
$\Delta\phi = \phi - \phi_c$, $\mu_w = 1$ mPa·s, and
$A = \frac{1-\sin\beta}{1+\sin\beta}e^{\pi\tan\beta}$ a granular bearing-capacity coefficient
($\beta_0 = 44\pm2°$, cornstarch's repose angle). $\mathrm{ID}\ll1$ means dilatancy wins and the impact
is solid-like; $\mathrm{ID}\gg1$ means inertia wins and you get their re-entrant liquid regime.

### What it gives for your silica
The one input nobody has measured for your material is $\phi_c$, so here it is as a function of $\phi_c$,
at $U_0 = 2.3$ m/s, $D_0 = 2.5$ mm, $\phi = 0.44$, $a = 800$ nm:

| $\phi_c$ | $\Delta\phi$ | $\kappa$ (m²) | **ID** | $t_d$ | Darcy $P_f$ (Pa) |
|---:|---:|---:|---:|---:|---:|
| 0.36 | 0.08 | $3.2\times10^{-15}$ | $1.2\times10^{-5}$ | 13 ns | $1.4\times10^{8}$ |
| 0.40 | 0.04 | $3.2\times10^{-15}$ | $2.7\times10^{-5}$ | 29 ns | $7.1\times10^{7}$ |
| 0.42 | 0.02 | $3.2\times10^{-15}$ | $5.5\times10^{-5}$ | 60 ns | $3.6\times10^{7}$ |
| 0.43 | 0.01 | $3.2\times10^{-15}$ | $1.1\times10^{-4}$ | 122 ns | $1.8\times10^{7}$ |
| 0.47 (Shah's value) | **−0.03** | — | **no dilatant branch at all** | — | 0 |

**Headline: ID $\approx 6\times10^{-5}$, about 220× smaller than cornstarch at matched conditions.**
The whole difference is $\kappa\propto a^2$: your particles are 12.5× smaller than cornstarch, so
$\kappa$ is 156× smaller, so the Darcy resistance to dilation is 156× larger. Colloidal silica is
*more* dilatant, not less, and by a lot.

Two consequences for the group:

- **Regime 3 is unreachable for you.** ID $=1$ would need $U_0\sim2\times10^5$ m/s on the literal
  formula. Whatever you see at high speed is not the re-entrant liquid the paper describes.
- **The last row is a real problem and it is informative.** If $\phi_c\approx0.47$, which is what
  Shah, Arora & Driscoll measured on these same 830 nm particles, then at $\phi = 0.44$ your
  suspension is *below* the dilatancy onset, $\Delta\phi<0$, and there is no dilatant branch for
  inertia to compete with. But your drops demonstrably arrest and grow solidification fronts at
  $\phi = 0.42$–0.45. So either your $\phi$ axis is offset from Shah's (which your own §8 landmine
  already says it is), or $\phi_c$ for your particles is below 0.42. **ID is only defined once
  $\phi_c$ is measured, and running it backwards from the fact that fronts exist at $\phi = 0.42$
  says $\phi_c \lesssim 0.42$.** That is a number the project does not currently have and can get.

The same tension shows up in the Wyart–Cates parameters independently: the fitted $\phi_M = 0.5176$
means that at $\phi = 0.44$ the model cannot shear-jam at any stress, only thicken continuously. The
repo's own closure class refuses to instantiate at the paper's operating point for exactly this reason
(pinned as a test).

![**Figure 5.** ID against particle size at fixed $\Delta\phi$. Solid: Eq. 6 as published, $\mathrm{ID}\propto a^2$. Dashed: capped at the pressure the pore fluid can sustain, which removes the size dependence entirely. The two vertical lines are your silica and the paper's cornstarch.](figures/fig5_ID_vs_particle_size.png)

### Sensitivities
$\beta_0$ enters exponentially and is unmeasured for silica (44° is cornstarch's repose angle):
ID varies from $1.0\times10^{-4}$ at $\beta_0=30°$ to $3.4\times10^{-5}$ at $51°$, a factor of 3.
$\phi_c$ contributes a factor of 9 over the plausible range. Neither is fatal for the *ratio* between
systems, which is what §6's headline uses.

---

## 7. A problem with Eq. 6, and a proposed fix

### The problem
Evaluate ID with the paper's own constants at the paper's own Regime-2/3 boundary conditions
(cornstarch, $\phi = 0.41$, $U_0 = 5.1$ m/s, $D_0 = 4$ mm, $a = 10$ µm, $\phi_c = 0.37$, $\beta_0=44°$):

$$\kappa = 6.8\times10^{-13}\ \mathrm{m^2},\quad \rho = 1226,\quad A = 4.0,\quad \Delta\phi = 0.04
\;\;\Longrightarrow\;\; \mathrm{ID} = 6.6\times10^{-3}$$

Not 1. This was computed twice, independently, and agrees to three figures. At the paper's top speed of
7.2 m/s it reaches only $9.3\times10^{-3}$. Getting ID $=1$ anywhere in their range would need
$\Delta\phi\approx3\times10^{-4}$, i.e. the boundary would lie on top of the $\phi_c$ line, but Fig. 1(d)
draws them as clearly separate curves.

The likely cause is stated in the paper itself: $F_d$ is written with "$\sim$" after setting
$r_c\sim h\sim D_0$ and $|dh/dt|\sim U_0$. Those substitutions inflate $F_d$ badly. Carried through,
they give $F_d\approx 240$ N against a measured peak force of order 0.5 N. **Their *fits* are
unaffected**, because Eq. 4 uses the measured $r_c(t)$ and $|dh/dt|$; it is only the scaling group that
is in question. The Supplemental Material, which is not on arXiv, may carry a prefactor.

**This is worth one email.** You have direct access to a co-author.

### The proposed fix, and why it is more than bookkeeping
Look at the last column of the table in §6: Eq. 3 predicts Darcy pore suctions of $10^7$–$10^8$ Pa for
800 nm particles. A pore fluid cannot do that. Two ceilings bind first:

- **cavitation** of the pore water at about one atmosphere of suction, $\sim1.0\times10^5$ Pa;
- **interface pull-in.** The free surface is dragged into the pores once suction exceeds $2\gamma/a$,
  which is $1.8\times10^5$ Pa at 800 nm. (This is Brown & Jaeger's "frustrated dilatancy against a
  stiff boundary", where the stiff boundary is the liquid–air interface.)

So for a submicron colloid the Darcy formula overshoots the physically available pressure by
**350×**, and even for cornstarch it overshoots by 83×. The empirical anchor is direct: Jérôme,
Vandenberghe & Forterre put a hydrophone into a dense bath under impact and measured pore suctions of
$-10$ to $-100$ kPa, **saturating at exactly one atmosphere.**

Cap the dilatant pressure at $P_{\rm ceil} = \min(1\ \mathrm{atm},\,2\gamma/a)$ and the number changes
character:

$$\mathrm{ID_{cap}} = \frac{\rho U_0^2}{\pi A\,P_{\rm ceil}}$$

Two structural consequences: it scales as $U_0^2$ rather than $U_0$, and **it stops depending on
particle size at all** once the cap binds. The $a^{-2}$ sensitivity that dominates the uncapped number
turns out to be an artefact of extrapolating Darcy past the cavitation limit.

### Does capping help?

| system | $U_0$ | ID (Eq. 6) | **ID capped** |
|---|---:|---:|---:|
| cornstarch, $\phi=0.41$, $D_0=4$ mm | 5.1 m/s | 0.0066 | **0.18** |
| cornstarch, same | 7.2 m/s (their top speed) | 0.0093 | **0.35** |
| cornstarch, ID $=1$ at | — | 771 m/s | **12.2 m/s** |
| your silica, $\phi_c=0.42$, $D_0=2.5$ mm | 2.3 m/s | $5.5\times10^{-5}$ | 0.0062 |
| your silica, ID $=1$ at | — | $1.9\times10^5$ m/s | **30 m/s** |

Regime 3 *is* observed in their experiments at the top of their velocity range. The capped number is
0.35 there, order unity and off by a factor of 3. The uncapped number is 0.0093, off by a factor of 100.
**The cap improves the location of the regime boundary by a factor of 30, on the one system where the
boundary has been measured.**

Be careful about how much weight this bears: it is one system, and a factor-3 agreement obtained from
a one-parameter fix is suggestive, not proof. The honest claim is that Eq. 6 as printed does not land
the boundary, that the Darcy pressure it contains is physically unreachable, and that fixing the second
problem largely fixes the first. The falsifiable prediction is that **800 nm silica should show
re-entrant liquid-like impact behaviour somewhere around 30 m/s**, reachable with a pressurised jet or
a tall tower, and out of reach entirely on the uncapped formula.

---

## 8. The shear-rate convention problem

This one affects a number your project already quotes, so it needs saying plainly.

Four conventions are in print for "the shear rate of an impact", three of them from Cheng's or
Michelle's groups. For your drop at 2.3 m/s:

| convention | source | value (1/s) |
|---|---|---:|
| (A) $U_0/D_0$ | Shah, Arora & Driscoll 2022, **your fitted axis** | 920 |
| (B) $U_0/2D_0$ | same paper, for $\mathrm{Re_{eff}}$ | 460 |
| (C) $U_0/\ell$, $\ell\approx100$ µm | same paper, for $\mathrm{Re_{eff}}$ | 23 000 |
| (D) $U_0/L_c$, $L_c = 0.986\,\mathrm{Re}^{\alpha(We)}D_0$ | Mobaseri, Kumar & Cheng, *PNAS* **122**, e2500163122 (2025) | 4 160 |

(D) is the only one that was validated rather than asserted: it is defined by matching the *measured*
viscous dissipation integral in >1200 simulations, and the paper explicitly refutes (A), showing
$L_c<0.45\,D_0$ always. For your conditions $L_c = 0.22\,D_0$, so **(D) is 4.5× (A)**.

What that does to your fitted onsets:

| salt | $\dot\gamma_c$ as fitted, on (A) | rescaled to (D) |
|---|---:|---:|
| 0 mM | 632 | 2860 |
| 10 mM | 479 | 2170 |
| 50 mM | 471 | 2130 |
| 100 mM | 196 | 885 |

And on the rheology side, from your own calibrated Wyart–Cates parameters:
$\sigma^\ast/\eta_0 = 1290$ 1/s, or 936 1/s using the self-consistent criterion $f=0.5$.

Which is awkward, and I am not going to smooth it over: **as fitted, on convention (A), your 632 1/s
agrees with the rheology-derived 936–1290 1/s to within a factor of 1.5–2. Rescaled onto the more
defensible convention (D), it disagrees by a factor of 2–3 in the other direction.** The convention
correction makes the agreement worse.

Three readings, and they are distinguishable:
1. $L_c$ was calibrated for Carreau fluids on superhydrophobic surfaces. Transferring it to a jamming
   colloid on hydrophilic glass may not be legitimate.
2. Your $\sigma^\ast$ and $\eta_0$ come from a flow-curve fit on a $\phi$ axis that is offset from your
   drops'. This is the §8 landmine again.
3. The two quantities are not the same object anyway: $L_c$ is a dissipation-matched average over the
   whole impact, while your $\dot\gamma_c$ is the intercept of a front-speed fit.

The actionable part is small and clear: **whenever the fitted $\dot\gamma_c$ is compared with anything,
state the convention.** A factor of 4.5 is larger than the salt effect you are measuring.

---

## 9. What this predicts, and how to test it

Ranked by cost.

1. **Free, from data you have.** The $\sigma^\ast$ ratios in §5: 1.00 / 0.66 / 0.68 / 0.17. Contains
   no fitted constant. Measure $\sigma^\ast$ at four salt levels on the rheometer and check.
2. **Free, from footage you have.** The frictional region should be an **annulus at the contact-line
   rim**, not a central disc, and it should fill inward as $t^*$ grows (Figure 4). Your under-view
   channel can see this.
3. **Cheap bench.** $U_c\propto\sigma^{\ast2/3}R^{1/3}\eta^{-1/3}$. The $R^{1/3}$ is a second needle:
   a 1.5× larger drop should raise the onset velocity by 14%. Small, but it is a sign test on a
   dependence nothing else in the project probes.
4. **One measurement that unlocks ID.** $\phi_c$ for your silica. Without it, ID is undefined for
   your system. And if $\phi_c$ really is 0.47, then ID says your suspension should not dilate at all
   at $\phi = 0.44$, which contradicts every front you have tracked.
5. **New capability, moderate cost.** A piezo force sensor. Predicted $F_{\max}=37$ mN at $t=171$ µs
   for a band-B drop, with the whole $F(t)$ curve predicted parameter-free. The 2026 paper found a
   regime with it that high-speed photography could not see.
6. **The big one, if you want it.** Re-entrant liquid-like behaviour in silica at $U_0\sim30$ m/s under
   the capped model, and never under Eq. 6 as printed. This discriminates between them cleanly.

---

## 10. What is uncertain, and by how much

Stated so you are never caught out in a talk.

| quantity | value used | range | effect on conclusions |
|---|---|---|---|
| $\sigma^\ast$ | 44 Pa | 19–100 Pa | $U_c$ moves 2.9× |
| $\eta_0$ | 34 mPa·s | 7–34 mPa·s | $U_c$ moves 1.7× ($\tau\propto\sqrt\eta$ saves you) |
| criterion ($\Lambda_F=1$ vs half-area) | $\Lambda_F=1$ | — | $U_c$ moves 1.8× |
| rim cut-off $\Delta$ | 1.6% of $d$ | unknown physically | rim stresses only; the area average is immune |
| $\phi_c$ (silica) | 0.42 | 0.36–0.47 | ID moves 9×; at 0.47, ID undefined |
| $\beta_0$ (silica) | 44° | 30–51° | ID moves 3× |
| $\rho_p$ | 2000 kg/m³ | 1900–2600 (**contested in your own files**) | $\rho$, hence every $\phi$ |
| Wagner prefactor on a deformable target | 3 | unknown $O(1)$ | law (b) amplitude only |

Two structural caveats that are not error bars:

- **Wyart–Cates is a steady-state mean-field model.** Your impact runs for ~1 ms with strain rates
  changing over orders of magnitude within that window. Comparing $\sigma$ to $\sigma^\ast$
  instantaneously assumes the frictional network equilibrates faster than the impact, and nobody has
  checked that for a drop. The correction would be a finite start-up strain $\gamma_0$, which appears in
  the shear-front literature and in no drop-impact paper.
- **Two frameworks are being mixed and should not be silently merged.** ID is a Reynolds–Darcy
  construct built from permeability, water viscosity, a friction angle and $\Delta\phi$. It contains no $\sigma^\ast$,
  no $F^\ast$, no $\phi_0/\phi_M$ interpolation. $\Lambda_F$ is a Wyart–Cates construct and contains no
  permeability. They answer different questions ("is the packing above the dilatancy onset and does
  Darcy have time to act?" versus "is the local stress above the frictional transition?") and there is
  currently no published bridge between them for a drop. Royer, Blair & Hudson found for silica that
  the dilatancy onset ($\phi\approx0.53$) and the DST onset ($\phi\approx0.58$) are *different*, and
  concluded both are consequences of friction rather than one driving the other, which is in direct
  tension with the ID paper's premise. Say this out loud rather than papering over it.

---

## 11. Presenting it

A twelve-slide arc that follows the logic above:

1. The question, in one line: what stress does an impact impose, where, and is it enough to thicken?
2. The four stages, with your drop's numbers.
3. **The water hammer is a red herring.** 5.4 MPa, 0.5 ns, one part in a billion of the drop.
4. Derive the Wagner pressure. Four lines. Emphasise $d\dot d = \tfrac32 RU_0$ being constant.
5. Read the three factors. **The maximum is a ring, not a spot.** Contrast with a steady jet.
   (Figure 1.)
6. The stress ladder. (Figure 2.)
7. Non-rigid substrates: two laws with the same structure, $Z_2/(Z_1+Z_2)$ and $\rho_2/(\rho_1+\rho_2)$,
   both giving one half for matched media. $\mathcal{C}$ for solids, $J$ as the switch.
   Your Carbopol ladder sits at $J\lesssim4\times10^{-3}$.
8. **Why $\rho U_0^2$ can't work**, and why the shear stress can. (Figure 3.) This is the pivot of
   the talk.
9. The frictional annulus, and the claim that it is the solidification front. (Figure 4.)
10. ID: what it is, whose paper it is, what it gives for silica, and the $a^2$ sensitivity.
    (Figure 5.)
11. The two problems: Eq. 6 doesn't land at 1, and Darcy predicts unreachable pressures. The cap, and
    what it buys.
12. What to test, ranked. Lead with the four prefactor-free $\sigma^\ast$ ratios.

**Three sentences to have ready.**

> *"The stress in a drop impact is concentrated in a ring at the contact line, and it decays as
> $t^{-1/2}$, so 'the impact stress' isn't a number, it's a field with a known closed form."*

> *"$\rho U^2$ is the wrong gauge for shear thickening by about two orders of magnitude, because
> dilatancy responds to shear stress and $\rho U^2$ is a pressure; the right gauge is smaller by
> $\sqrt{\mathrm{Re}}$ and it crosses $\sigma^\ast$ right where we measure the onset."*

> *"ID for our silica is about $6\times10^{-5}$, roughly 220× more dilatant than the cornstarch the
> number was built for, purely because permeability goes as $a^2$. But the Darcy pressure behind it
> exceeds the cavitation limit by two orders of magnitude, so I think the number needs a cap."*

**If asked "is this published?"** The stress field, yes, thoroughly. The impedance/inertia pairing is
a repackaging of things in print. The shear-stress gauge $\Lambda_F$, the $\tau=\sigma^\ast$ contour,
and the capped ID are new here and have not been refereed.

**If asked "how sure are you?"** Point at the table in §10. The phrase to use is *"the scaling is solid, the
prefactors are order-unity uncertain, and the ratios are the part that carries no fitted constants."*

---

## Appendix A: symbols

| symbol | meaning | units |
|---|---|---|
| $R$, $D_0=2R$ | drop radius, diameter | m |
| $U_0$ | impact speed | m/s |
| $t^* = U_0t/R$ | dimensionless time, **radius clock** | — |
| $\tilde t = U_0t/D_0 = t^*/2$ | dimensionless time, **diameter clock** (Cheng's) | — |
| $d(t)=\sqrt{3RU_0t}$ | wetted / turning-line radius | m |
| $\xi = r/d(t)$ | self-similar radial coordinate | — |
| $Z=\rho c$ | acoustic impedance | Pa·s/m |
| $\mathrm{Re}_R = \rho U_0R/\eta$ | radius-based Reynolds number ($=\mathrm{Re}/2$) | — |
| $\sigma^\ast$ | dilatancy / frictional onset stress | Pa |
| $\Lambda_F$ | $\langle\tau\rangle/\sigma^\ast$, the proposed shear gauge | — |
| $\kappa$ | permeability, $(1-\phi)^3a^2/180\phi^2$ | m² |
| $\Delta\phi = \phi-\phi_c$ | distance past the dilatancy onset | — |
| ID | inertial–dilatant number | — |
| $\mathcal{C} = \rho U_0^2/E^\ast$ | Cauchy number (substrate compliance) | — |
| $J = \tau_y/\rho U_0^2$ | yield number | — |

**Watch the two traps.** (i) Every factor-of-2 confusion in this subject is the $t^*$ / $\tilde t$
clock. (ii) $a$ is a *radius* in the Wyart–Cates literature and a *diameter* in Carman–Kozeny; the code
keeps them in separate classes for that reason.

## Appendix B: sources actually read

Primary, full text retrieved and used:
Cheng, Sun & Gordillo, *Annu. Rev. Fluid Mech.* **54**, 57 (2022) ·
Philippi, Lagrée & Antkowiak, *JFM* **795**, 96 (2016), arXiv:1504.05847 ·
Mobaseri, Gordillo, Burton, Yoon, Lee, Kumar, Driscoll & Cheng, *PRL* **136**, 148201 (2026),
arXiv:2601.12642 ·
Mobaseri, Kumar & Cheng, *PNAS* **122**, e2500163122 (2025), arXiv:2503.23555 ·
Shah, Arora & Driscoll, *Commun. Phys.* **5**, 275 (2022) ·
Sun, Álvarez-Novoa, Andrade, Gutiérrez, Gordillo & Cheng, *Nat. Commun.* **13**, 1703 (2022) ·
Gordillo, Sun & Cheng, *JFM* **840**, 190 (2018) ·
Howland *et al.*, *PRL* **117**, 184502 (2016) ·
Riboux & Gordillo, *PRL* **113**, 024507 (2014) ·
Mandre & Brenner, *JFM* **690**, 148 (2012) ·
Jérôme, Vandenberghe & Forterre, *PRL* **117**, 098003 (2016) ·
Jørgensen, Forterre & Lhuissier, *JFM* **896**, A28 (2020) ·
Pegg, Purvis & Korobkin, *JFM* **839**, 561 (2018) ·
Langley, Castrejón-Pita & Thoroddsen, *Soft Matter* **16**, 5702 (2020) ·
Bisighini, Cossali, Tropea & Roisman, *PRE* **82**, 036319 (2010) ·
Royer, Blair & Hudson, *PRL* **116**, 188301 (2016) ·
Guy, Hermes & Poon, *PRL* **115**, 088304 (2015) ·
Brown & Jaeger, *Rep. Prog. Phys.* **77**, 046602 (2014) ·
Ness, Seto & Mari, *Annu. Rev. Condens. Matter Phys.* **13**, 97 (2022) ·
Athani, Metzger, Forterre & Mari, *PRF* **10**, 043301 (2025) ·
Yokoyama, Maruoka, Matsunuma & Tagawa, *Nat. Commun.* **17**, 607 (2026).

Two citation corrections worth carrying: the short-time self-similarity paper is
Philippi, Lagrée & **Antkowiak** (not Josserand); and the term is **"inertial–dilatant"**, not
"inertio-dilatant", which appears nowhere in the literature.

## Appendix C: how to rerun

```bash
cd C:\Users\compl\Desktop\Projects\impact_stress
python -m pytest test_impact_stress.py -q     # 40 tests
python results.py                             # out/results.json + figures/
```

| file | holds |
|---|---|
| `materials.py` | fluids, Wyart–Cates, your calibrated silica, provenance and caveats |
| `stress.py` | the four stages; Wagner pressure, wall shear, force, compliance laws |
| `dilatancy.py` | $\Lambda_P$, $\Lambda_F$, the $\tau=\sigma^\ast$ contour, onset velocity |
| `inertial_dilatant.py` | ID as published, plus the capped extension |
| `conventions.py` | the four shear-rate conventions and $L_c$ |
| `results.py` | every number and figure in this document |
| `out/results.json` | the numbers, machine-readable |

Validation, in increasing order of what it proves: limits (rigid to $\rho cU$, matched to $\tfrac12\rho cU$);
identities (dimensional and self-similar forms agree); **literature** (Cheng's Eqs. [10], [13] coded
independently from the paper and asserted equal to ours); **cross-checks** (numerical integration of the
shear field reproduces the published drag force to 11%, with the residual traced to the empirical
cut-off; the compressible crossover from two unrelated arguments agrees to 1.6×; our Wyart–Cates
matches the repo's own closure to $10^{-12}$).

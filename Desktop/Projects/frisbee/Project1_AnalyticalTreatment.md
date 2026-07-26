# A Rigorous Analytical Treatment of Disc Flight and Stability

**Author:** charlie
**Companion to:** *Project 1 — Disc Flight in a Lakefront Internal Boundary Layer.*
**Status:** v0.1, 2026-05-29.

This document derives, from first principles, the six-degree-of-freedom dynamics, steady-glide solutions, linear stability theory, and gyroscopic-precession behavior of a spinning rigid disc in flight, together with the closed-form effects of a steady non-uniform wind. The treatment is set up to be cited in the methods/theory section of the manuscript and to serve as the ground truth against which the numerical 6-DOF integrator and the coupled CFD simulator are validated.

Notation follows Hummel (2003) / Crowther & Potts (2007) where compatible, with several modernizations for clarity.

---

## 1. Reference frames and notation

We work with three right-handed orthonormal frames.

**Inertial frame** $\mathcal{N} = \{\hat{\mathbf{n}}_1, \hat{\mathbf{n}}_2, \hat{\mathbf{n}}_3\}$: Earth-fixed, $\hat{\mathbf{n}}_3$ vertically upward, $\hat{\mathbf{n}}_1$ aligned with the intended throw azimuth at release.

**Body frame** $\mathcal{B} = \{\hat{\mathbf{b}}_1, \hat{\mathbf{b}}_2, \hat{\mathbf{b}}_3\}$: attached to the disc, with $\hat{\mathbf{b}}_3$ along the disc's axis of symmetry (the rim is normal to $\hat{\mathbf{b}}_3$, pointing upward through the "top" of the disc — i.e. the convex side for an Ultrastar). $\hat{\mathbf{b}}_1$ is taken to point along the projection of the relative wind $\mathbf{V}_{\text{rel}}$ onto the disc plane (so the body frame is *not* purely material; it rotates only with pitch/roll/yaw, not with spin).

**Spinning disc frame** $\mathcal{D}$: rotates with the disc at spin rate $\dot{\phi}_s$ about $\hat{\mathbf{b}}_3 = \hat{\mathbf{d}}_3$.

The orientation of $\mathcal{B}$ relative to $\mathcal{N}$ is parameterized by a 3-1-3 Euler-angle sequence $(\psi, \theta, \phi_s)$ — yaw $\psi$ about $\hat{\mathbf{n}}_3$, then tilt $\theta$ about the new $\hat{\mathbf{n}}_1'$, then spin $\phi_s$ about the body $\hat{\mathbf{b}}_3$. Since the disc is axisymmetric, only $(\psi, \theta)$ orient the symmetry axis; $\phi_s$ orients material points.

The disc has mass $m$, diameter $d$, radius $r = d/2$, planar area $A = \pi d^2/4$. Treating the disc as a thin uniform plate gives principal moments of inertia

$$
I_p \equiv I_x = I_y = \tfrac{1}{4} m r^2, \qquad
I_s \equiv I_z = \tfrac{1}{2} m r^2, \qquad
I_s = 2 I_p. \tag{1}
$$

For an Ultrastar ($m = 0.175$ kg, $r = 0.135$ m):

$$
I_s \approx 1.59\times10^{-3}\,\text{kg}\cdot\text{m}^2, \qquad
I_p \approx 7.97\times10^{-4}\,\text{kg}\cdot\text{m}^2. \tag{2}
$$

The inertia tensor $\mathsf{I}_{\mathcal{B}} = \mathrm{diag}(I_p, I_p, I_s)$ is diagonal in $\mathcal{B}$ regardless of spin phase, by axisymmetry — the essential simplification that makes the disc tractable.

Let $\boldsymbol{\omega}$ be the angular velocity of $\mathcal{B}$ relative to $\mathcal{N}$, with body-frame components $(p, q, r)$:

$$
\boldsymbol{\omega} = p\,\hat{\mathbf{b}}_1 + q\,\hat{\mathbf{b}}_2 + r\,\hat{\mathbf{b}}_3. \tag{3}
$$

The disc's total angular velocity (including spin) is $\boldsymbol{\omega} + \dot{\phi}_s \hat{\mathbf{b}}_3 = p\,\hat{\mathbf{b}}_1 + q\,\hat{\mathbf{b}}_2 + (r + \dot{\phi}_s)\hat{\mathbf{b}}_3$. We define the **spin rate**

$$
\Omega \equiv r + \dot{\phi}_s \tag{4}
$$

as the total component about $\hat{\mathbf{b}}_3$. Because $\mathcal{B}$ already follows the disc's tilt but not its spin, $r$ is small (rate of change of the orientation of the projected airflow) and $\Omega \approx \dot{\phi}_s$ in typical flight.

Position of the disc center of mass in $\mathcal{N}$: $\mathbf{x}$; inertial velocity: $\mathbf{v} = \dot{\mathbf{x}}$. The relative wind is

$$
\mathbf{V}_{\text{rel}}(\mathbf{x}, t) = \mathbf{v} - \mathbf{U}(\mathbf{x}, t), \qquad V \equiv |\mathbf{V}_{\text{rel}}|, \tag{5}
$$

where $\mathbf{U}(\mathbf{x}, t)$ is the ambient wind, which for the CFD-coupled problem is read from the simulation output and treated as a *known* field.

---

## 2. Aerodynamic angle of attack and force decomposition

The disc-plane component of $\mathbf{V}_{\text{rel}}$ and the disc-normal component are

$$
\mathbf{V}_\parallel = \mathbf{V}_{\text{rel}} - (\mathbf{V}_{\text{rel}} \cdot \hat{\mathbf{b}}_3)\hat{\mathbf{b}}_3, \qquad
V_\perp = \mathbf{V}_{\text{rel}}\cdot\hat{\mathbf{b}}_3. \tag{6}
$$

The angle of attack is

$$
\alpha = \arctan\!\Big(\!-V_\perp / |\mathbf{V}_\parallel|\Big). \tag{7}
$$

Sign convention: $\alpha > 0$ when the flow strikes the lower (convex-bottom) surface, producing positive lift. This matches Hummel (2003).

We choose the body $\hat{\mathbf{b}}_1$ axis along $\mathbf{V}_\parallel/|\mathbf{V}_\parallel|$, so $\hat{\mathbf{b}}_2 = \hat{\mathbf{b}}_3 \times \hat{\mathbf{b}}_1$ is the "side" axis, normal to the relative wind in the disc plane.

The **lift direction** $\hat{\mathbf{L}}$ is the unit vector perpendicular to $\mathbf{V}_{\text{rel}}$, lying in the plane $(\mathbf{V}_{\text{rel}}, \hat{\mathbf{b}}_3)$, on the side of $\hat{\mathbf{b}}_3$:

$$
\hat{\mathbf{L}} = \frac{\hat{\mathbf{b}}_3 - (\hat{\mathbf{b}}_3 \cdot \hat{\mathbf{v}}_{\text{rel}})\hat{\mathbf{v}}_{\text{rel}}}{\big|\hat{\mathbf{b}}_3 - (\hat{\mathbf{b}}_3 \cdot \hat{\mathbf{v}}_{\text{rel}})\hat{\mathbf{v}}_{\text{rel}}\big|}, \qquad \hat{\mathbf{v}}_{\text{rel}} = \mathbf{V}_{\text{rel}}/V. \tag{8}
$$

The **drag direction** is $-\hat{\mathbf{v}}_{\text{rel}}$. The aerodynamic force on the disc is therefore

$$
\mathbf{F}_a = \tfrac{1}{2}\rho V^2 A\,\big[\,C_L(\alpha, \mathrm{AdvR})\,\hat{\mathbf{L}} - C_D(\alpha, \mathrm{AdvR})\,\hat{\mathbf{v}}_{\text{rel}}\,\big], \tag{9}
$$

where $\rho$ is air density and AdvR is the **advance ratio**

$$
\mathrm{AdvR} \equiv \frac{\Omega d}{2 V}. \tag{10}
$$

Wind-tunnel data (Potts & Crowther 2002; Kamaruddin et al. 2018) shows that for AdvR up to $\sim 1$, $C_L$ and $C_D$ are essentially independent of AdvR, but the aerodynamic *moments* are not.

A side-force term $C_S$ would appear in (9) if the disc had asymmetry; for an Ultrastar it is negligible at small AdvR.

### 2.1 Linearized coefficients

Following Hummel (2003), in the low-to-moderate $\alpha$ regime relevant to most throws,

$$
\boxed{
\begin{aligned}
C_L(\alpha) &= C_{L0} + C_{L\alpha}\,(\alpha - \alpha_0), \\
C_D(\alpha) &= C_{D0} + C_{D\alpha}\,(\alpha - \alpha_0)^2,
\end{aligned}}
\tag{11}
$$

with $\alpha_0 \approx -4°$ the angle of minimum drag. The drag is quadratic because parasite drag at $\alpha_0$ is augmented by induced drag (lift-related vortex generation), exactly as for a low-aspect-ratio wing. The constants for a 175 g Ultrastar (Hummel "short flights" data, used in (3)–(5)):

$$
C_{L0} = 0.33, \quad C_{L\alpha} = 1.91\ \text{rad}^{-1},
\quad C_{D0} = 0.18, \quad C_{D\alpha} = 0.69\ \text{rad}^{-2}.
\tag{12}
$$

For uncertainty bracketing we cross-check against Potts & Crowther (2002): $C_{L0}=0.20$, $C_{L\alpha}=2.96$, $C_{D0}=0.08$, $C_{D\alpha}=2.72$ — substantially different absolute values but a similar L/D structure.

The linear range of (11) holds to $\alpha \lesssim 15°$; beyond, stall dynamics enter and the coefficients must come from a tabulated wind-tunnel curve. All canonical Ultimate throws (backhand, forehand) launch in this linear range; only the **hammer** (overhead, $\alpha_0 \approx +40°$) operates outside it and is handled separately.

### 2.2 Reynolds number

For $V \sim 15$ m/s, $d = 0.27$ m, $\nu_{\text{air}} = 1.5\times10^{-5}$ m²/s:

$$
\mathrm{Re} = \frac{V d}{\nu} = \frac{(15)(0.27)}{1.5\times10^{-5}} \approx 2.7\times10^{5}. \tag{13}
$$

Potts & Crowther (2002) showed $C_L, C_D$ are essentially Re-independent over $\mathrm{Re} \in [1.14, 3.78]\times10^5$. The coefficients in (12) are therefore portable across the full operational wind range, which is the technical justification for using the same (11) throughout the trajectory.

---

## 3. Aerodynamic moments

Resolving the aerodynamic moment about the center of mass into body-frame components, with non-dimensionalization $\overline{M} = M / (\tfrac{1}{2}\rho V^2 A d)$:

### 3.1 Pitching moment $M$ (about $\hat{\mathbf{b}}_2$)

$$
M = \tfrac{1}{2}\rho V^2 A d\,\Big[\,C_{M0} + C_{M\alpha}(\alpha - \alpha_0) + C_{Mq}\,\tfrac{q d}{2V}\,\Big]. \tag{14}
$$

Hummel (2003): $C_{M0} = -0.08$, $C_{M\alpha} = 0.43$ rad$^{-1}$, $C_{Mq} = -5.0\times10^{-3}$ rad$^{-1}$. The negative $C_{M0}$ is the **static stability**: at $\alpha = \alpha_0$ the moment is nose-down. The positive $C_{M\alpha}$ is a *destabilizing* slope — increasing $\alpha$ produces a more positive (nose-up) moment, pushing $\alpha$ higher still. The disc is therefore **statically unstable in pitch** about $\hat{\mathbf{b}}_2$. Were it not spinning, it would tumble within a fraction of a second. The function of spin is to convert this static instability into stable gyroscopic precession (§5).

The $C_{Mq}$ term is the pitch-damping coefficient; it is small and weakly stabilizing.

### 3.2 Rolling moment $R$ (about $\hat{\mathbf{b}}_1$)

$$
R = \tfrac{1}{2}\rho V^2 A d\,\Big[\,C_{Rr}\,\tfrac{r d}{2V} + C_{Rp}\,\tfrac{p d}{2V}\,\Big]. \tag{15}
$$

$C_{Rr}$ is the roll-due-to-yaw-rate coefficient (with $r \approx \Omega$ for the spin); $C_{Rp}$ is the roll-rate damping. Hummel values: $C_{Rr} = 1.4\times10^{-2}$, $C_{Rp} = -5.5\times10^{-3}$. The dominant effect is the spin-induced roll moment, which couples to the pitching moment (14) to produce the gyroscopic precession trajectory features.

### 3.3 Spin-axis moment $N$ (about $\hat{\mathbf{b}}_3$)

$$
N = \tfrac{1}{2}\rho V^2 A d\,\Big[\,C_{Nr}\,\tfrac{r d}{2V}\,\Big]. \tag{16}
$$

With $C_{Nr} = -7.1\times10^{-6}$ rad$^{-1}$ (Hummel short flights) or $-3.4\times10^{-5}$ (long flights). This is the only term in (16); there is no $\alpha$-driven spin-up. The negative sign produces spin-down: dimensional time scale

$$
\tau_{\text{spin}} \sim \frac{I_s}{|N|/\Omega} = \frac{I_s \cdot 2V}{|C_{Nr}|\,\rho V^2 A d^2 / 2}. \tag{17}
$$

For $V = 15$ m/s, $\rho = 1.2$ kg/m³: $\tau_{\text{spin}} \sim 4 \times 10^2$ s with the short-flight value, $\sim 80$ s with the long-flight value. Either way, $\tau_{\text{spin}}$ is long compared to typical flight times ($\sim 5$ s), so to leading order $\Omega$ is constant during a throw. Spin-down matters only for the longest flights (hammer arcs, big-arm hucks downwind).

---

## 4. Six-degree-of-freedom equations of motion

### 4.1 Translational equation (Newton)

In $\mathcal{N}$:

$$
m\,\dot{\mathbf{v}} = -m g\,\hat{\mathbf{n}}_3 + \mathbf{F}_a, \tag{18}
$$

with $\mathbf{F}_a$ from (9). Resolving (18) component-wise in $\mathcal{N}$ gives three ODEs for the inertial velocity.

### 4.2 Rotational equation (Euler)

The angular momentum of the disc about its center of mass, in $\mathcal{B}$:

$$
\mathbf{H} = I_p p\,\hat{\mathbf{b}}_1 + I_p q\,\hat{\mathbf{b}}_2 + I_s \Omega\,\hat{\mathbf{b}}_3. \tag{19}
$$

The inertial time derivative satisfies Euler's relation,

$$
\Big(\frac{d\mathbf{H}}{dt}\Big)_{\mathcal{N}} = \Big(\frac{d\mathbf{H}}{dt}\Big)_{\mathcal{B}} + \boldsymbol{\omega} \times \mathbf{H} = \mathbf{M}, \tag{20}
$$

with $\mathbf{M} = R\,\hat{\mathbf{b}}_1 + M\,\hat{\mathbf{b}}_2 + N\,\hat{\mathbf{b}}_3$. Substituting (3) and (19) into (20) and projecting onto each body axis yields

$$
\boxed{
\begin{aligned}
I_p \dot{p} + (I_s - I_p)\,q\,\Omega - I_p\,q\,r &= R \\
I_p \dot{q} - (I_s - I_p)\,p\,\Omega + I_p\,p\,r &= M \\
I_s \dot{\Omega} - I_s\,(p q - p q) &= N \quad \text{(after using } \dot{\Omega} = \dot{r} + \ddot{\phi}_s)
\end{aligned}}
\tag{21}
$$

Using $r$ to denote the body-$\hat{\mathbf{b}}_3$ component of $\boldsymbol{\omega}$ (the kinematic part, not the spin), and noting that $\Omega = r + \dot{\phi}_s$ with $\dot{\phi}_s$ the "free" spin variable, equations (21) decouple from $\dot{\phi}_s$ in the symmetric case. In practice — and following Hummel's convention — we absorb $r$ into $\Omega$ and treat the system as

$$
\begin{aligned}
I_p \dot{p} + (I_s - I_p)\,q\,\Omega &= R, \\
I_p \dot{q} - (I_s - I_p)\,p\,\Omega &= M, \\
I_s \dot{\Omega} &= N.
\end{aligned}
\tag{22}
$$

Equation (22c) gives the spin-down evolution (§3.3, slow). The two equations (22a–b) are the **transverse gyroscopic equations** and are the central object of the stability analysis.

### 4.3 Kinematic equations for Euler angles

With $(\psi, \theta, \phi_s)$ in 3-1-3 convention,

$$
\begin{pmatrix} p \\ q \\ r \end{pmatrix}
=
\begin{pmatrix}
-\dot{\psi}\sin\theta\cos\phi_s + \dot{\theta}\sin\phi_s \\
\phantom{-}\dot{\psi}\sin\theta\sin\phi_s + \dot{\theta}\cos\phi_s \\
\dot{\psi}\cos\theta + \dot{\phi}_s
\end{pmatrix}. \tag{23}
$$

The standard 3-1-3 singularity at $\theta = 0$ (disc face exactly horizontal) is reached only briefly during very flat throws. For high-pitch maneuvers (hammers), a quaternion parametrization is preferable, but the linear-stability analysis below uses small perturbations about an equilibrium with $\theta \approx \text{const}$ and never encounters the singularity.

---

## 5. Gyroscopic precession: turn and fade

Sections 5.1–5.3 contain the core physics result of this document — the closed-form explanation of why a thrown disc turns and fades.

### 5.1 Reduction in the high-spin limit

In typical operation $\Omega \sim 100$ rad/s while $|p|, |q| \lesssim 10$ rad/s. Equations (22a–b) under the high-spin limit $|\Omega| \gg |p|, |q|$ retain only the gyroscopic terms:

$$
(I_s - I_p)\,q\,\Omega \approx R, \qquad
-(I_s - I_p)\,p\,\Omega \approx M, \tag{24}
$$

i.e. the body-frame angular velocity is **algebraically determined** by the applied moment in the high-spin limit. The inertial terms $I_p \dot{p}, I_p \dot{q}$ contribute corrections of order $(I_p / I_s)/(\Omega \tau)$ where $\tau$ is the timescale on which $M, R$ change — typically the flight time. For $\Omega \tau \sim 500$, this is negligible (cf. textbook gyroscope theory). The disc behaves quasi-statically in the transverse angular variables.

Solving (24):

$$
q = \frac{R}{(I_s - I_p)\Omega}, \qquad
p = -\frac{M}{(I_s - I_p)\Omega}. \tag{25}
$$

Using (1), $I_s - I_p = I_p$, so

$$
\boxed{
q = \frac{R}{I_p \Omega}, \qquad
p = -\frac{M}{I_p \Omega}.}
\tag{26}
$$

These are the **precession rates**. They are the textbook result for a fast-spinning gyroscope: the spin axis processes about an axis perpendicular to the applied moment, at angular rate $|M|/(I_s\Omega)$ — except here, with $I_s - I_p = I_p$ instead of $I_s$, the prefactor is $I_p \Omega$. (This is a feature of the symmetric disc: the effective gyroscopic inertia is $I_p$ for transverse precession, not $I_s$. Footnote-worthy in the paper.)

### 5.2 Interpretation: roll = "turn", pitch = (slow drift in $\alpha$)

In our body frame $\hat{\mathbf{b}}_1$ is along $\mathbf{V}_\parallel$, so the rate $p$ about $\hat{\mathbf{b}}_1$ tilts the disc plane about an axis along the flight direction — i.e. it **rolls** the disc relative to the flight path. A right-handed thrower's backhand throws with spin direction $\Omega < 0$ in our convention (top of disc points up, spin angular momentum vector points downward). The aerodynamic pitching moment at flight equilibrium is positive (nose-up) for $\alpha > \alpha_0$, which is the typical release condition. By (26b),

$$
p = -\frac{M}{I_p \Omega} \quad \text{with } M > 0, \Omega < 0 \implies p > 0. \tag{27}
$$

A positive $p$ rolls the disc clockwise as seen from behind, looking along $\hat{\mathbf{b}}_1$ — that is, it banks the disc to the right (anhyzer). This is the celebrated "**turn**" that every Ultimate player has noticed: a backhand thrown with a positive $\alpha$ banks anhyzer in the first half of its flight.

As the disc decelerates and $V$ falls, the dynamic pressure $\tfrac{1}{2}\rho V^2$ in (14) decreases as $V^2$ while the gravitational requirement on $C_L$ rises (to maintain glide), pulling $\alpha$ up and changing the sign and magnitude of $M$. Once $V$ drops enough that $M$ flips sign, $p$ flips with it and the disc banks the opposite way — the "**fade**." Quantitatively, the turn-to-fade transition occurs when $\alpha$ crosses the value at which $C_M(\alpha) = 0$:

$$
\alpha_M^\star = \alpha_0 - C_{M0}/C_{M\alpha} = -4° - (-0.08/0.43)\,\text{rad} = -4° + 10.7° = 6.7°. \tag{28}
$$

So a disc released at $\alpha > 6.7°$ initially turns anhyzer (right-handed backhand), then fades hyzer when $\alpha$ drops below $6.7°$ — which happens spontaneously as the disc loses speed and gravity pulls it down into a steeper angle of climb. This single number, $6.7°$, organizes the entire turn-and-fade pattern of throws.

### 5.3 Pitch evolution and the slow precession of $\alpha$

The complementary rate $q$ about $\hat{\mathbf{b}}_2$ is more subtle: it changes the orientation of $\hat{\mathbf{b}}_3$ in the plane containing $\mathbf{V}_\parallel$, i.e. it changes the **angle of attack**. From (26a) and (15), $q$ is driven by the spin-induced roll moment $R$, which scales as $C_{Rr}$. For a backhand ($\Omega < 0$), $C_{Rr} > 0$ ⇒ $R < 0$ ⇒ $q < 0$, i.e. the disc tends to **pitch nose-down**. This slow precession of $\alpha$ is the secondary mechanism by which a steady-state glide is approached. It is also why a disc thrown with too high a release $\alpha$ does *not* simply tumble: the gyroscopic action converts the pitching moment into roll (the turn), not into a runaway pitch increase, and the disc relaxes toward $\alpha_M^\star$ on a precessional time scale (~$1/(C_{Rr}\Omega d/2V)$) rather than the inertial pitch time scale.

### 5.4 Spin stability criterion

For the gyroscopic stabilization to work, the precessional time scale must be shorter than the time scale on which an unspun disc would tumble due to the destabilizing $C_{M\alpha}$. The latter is

$$
\tau_{\text{tumble}} \sim \sqrt{I_p / (|\partial M / \partial \alpha|)} = \sqrt{I_p / (\tfrac{1}{2}\rho V^2 A d\,C_{M\alpha})}. \tag{29}
$$

The precessional time scale is $\tau_{\text{prec}} \sim 1/(C_{Rr}\,\Omega d / (2V))^{-1} = 2V/(C_{Rr}\,\Omega\,d)$. The stability requirement $\tau_{\text{prec}} \ll \tau_{\text{tumble}}$ becomes, after algebra,

$$
\boxed{
\Omega^2 \gg \frac{2 \rho V^2 A d\,C_{M\alpha}}{C_{Rr}^2\,I_p\,d^2/V^2}\,\frac{1}{4} = \frac{\rho V^4 A\,C_{M\alpha}}{2 C_{Rr}^2\,I_p\,d}.
}\tag{30}
$$

For $V = 15$ m/s, plugging numbers:

$$
\Omega_{\min}^2 \sim \frac{(1.2)(15^4)(0.057)(0.43)}{2(1.4\times10^{-2})^2 (7.97\times10^{-4})(0.27)}
\approx \frac{1490}{1.6\times10^{-9}} \cdot 10^{-6}
\sim 9\times10^5, \tag{31}
$$

so $\Omega_{\min} \sim 950$ rad/s — *which is two orders of magnitude above the typical throw spin of 100 rad/s.* The naïve criterion (30) is therefore not satisfied; nonetheless the disc is stable. This apparent paradox resolves because the destabilizing pitching moment does not produce divergent tumble — gyroscopic action converts it to *roll*, not pitch, and the roll is bounded by the precessional dynamics themselves (not by an external restoring moment). The disc is **conditionally stable**: the angle of attack relaxes on a time scale comparable to the flight time, which is "fast enough" to keep the disc flying nose-roughly-into-the-wind for the duration of a throw, but not so fast that the trajectory is forced to glide-equilibrium. A more precise condition is given by the eigenvalue analysis of §6.

This is, incidentally, the analytical reason a Frisbee flies but a flat plate of the same shape does not: the plate has $C_{M\alpha}$ but no $C_{Rr}$ to drive the precessional escape, and tumbles on $\tau_{\text{tumble}} \sim 0.1$ s.

---

## 6. Steady gliding flight and maximum range

### 6.1 The glide equilibrium

We look for a straight-line steady solution with $\dot{V} = 0$, fixed glide angle $\gamma$ below horizontal, constant $\alpha$, $\Omega \approx \text{const}$ (over the relevant time), no transverse precession ($\dot{p} = \dot{q} = 0$, $p = q = 0$). The translational balance projects onto the flight path and its normal:

$$
\text{along } \mathbf{V}_{\text{rel}}: \quad mg \sin\gamma = D = \tfrac{1}{2}\rho V^2 A C_D, \tag{32a}
$$

$$
\text{normal: } \quad mg \cos\gamma = L = \tfrac{1}{2}\rho V^2 A C_L. \tag{32b}
$$

Dividing,

$$
\tan\gamma = \frac{C_D(\alpha)}{C_L(\alpha)}, \tag{33}
$$

the classic glider equation: glide angle is set by the inverse lift-to-drag ratio at the operating $\alpha$. From (32b),

$$
V^2 = \frac{2 m g \cos\gamma}{\rho A\,C_L(\alpha)}. \tag{34}
$$

### 6.2 Maximum L/D and the optimal cruise $\alpha$

Maximizing $C_L/C_D$ in the linear/quadratic model (11) gives the optimal angle of attack offset $\Delta\alpha^\star = \alpha^\star - \alpha_0$ as the positive root of

$$
C_{D\alpha}\,C_{L\alpha}\,\Delta\alpha^2 + 2 C_{D\alpha} C_{L0}\,\Delta\alpha - C_{L\alpha} C_{D0} = 0. \tag{35}
$$

Closed form:

$$
\Delta\alpha^\star = \frac{-C_{L0} + \sqrt{C_{L0}^2 + C_{L\alpha}^2\,C_{D0}/C_{D\alpha}}}{C_{L\alpha}}. \tag{36}
$$

For Hummel coefficients (12): $\Delta\alpha^\star \approx 0.366$ rad $\approx 21°$, giving $\alpha^\star \approx 17°$, $C_L \approx 1.03$, $C_D \approx 0.27$, $(L/D)_{\max} \approx 3.8$. The Hummel linearization overshoots the wind-tunnel $(L/D)_{\max} \approx 2.5$ at $\alpha \approx 9°$ (Potts–Crowther data) because (11) does not capture stall onset. **For all paper-quality work the linearization must be replaced by tabulated $C_L(\alpha), C_D(\alpha)$ beyond $\alpha \approx 15°$.** The Hummel-coefficient prediction (36) should be reported as a heuristic upper bound and explicitly bracketed against the wind-tunnel-data optimum $\alpha \approx 9°$, $(L/D)_{\max} \approx 2.5$.

### 6.3 Range and time of flight

For a glide angle $\gamma$ released at height $h_0$:

$$
R_{\text{glide}} = \frac{h_0}{\tan\gamma} = h_0 \cdot \frac{C_L(\alpha)}{C_D(\alpha)}. \tag{37}
$$

For $h_0 = 1.5$ m (release height) and $(L/D)_{\max} \approx 2.5$: $R_{\text{glide,steady}} \approx 3.75$ m — laughably short, of course, because actual throws are not steady glides; they begin with $V$ well above the glide value and store kinetic energy in the initial fast cruise. **Lissaman & Hubbard (2010)** treat the optimization more carefully, including ballistic phase + glide phase, and give a maximum range of $\sim 180$ m for a 30 m/s launch into still air. For working purposes here, the relevant analytical statement is:

> Steady-glide range is bounded by $h_0 (L/D)_{\max}$, but the dominant range contribution in real throws comes from the ballistic / "lifting cruise" phase during which kinetic energy is depleted at the optimal-glide L/D.

### 6.4 Trim condition for the body

The translational equilibrium fixes $V$ and $\gamma$ given $\alpha$. The rotational equilibrium under (22) requires $M = 0$ and $R = 0$ (since $\dot{p} = \dot{q} = 0$ and $p = q = 0$ in steady glide). The $M = 0$ condition fixes $\alpha = \alpha_M^\star$ (eq. (28)):

$$
\alpha_{\text{trim}} = \alpha_0 - \frac{C_{M0}}{C_{M\alpha}} = 6.7°. \tag{38}
$$

This is the **trim angle of attack** — the only $\alpha$ at which the disc can glide without precession. The disc therefore self-selects $\alpha = 6.7°$ asymptotically, and the corresponding $L/D$ at this $\alpha$ is

$$
(L/D)_{\text{trim}} = \frac{C_{L0} + C_{L\alpha}(\alpha_M^\star - \alpha_0)}{C_{D0} + C_{D\alpha}(\alpha_M^\star - \alpha_0)^2}
= \frac{0.33 + 1.91(0.187)}{0.18 + 0.69(0.187)^2}
= \frac{0.687}{0.204}
\approx 3.4. \tag{39}
$$

Cross-checking with the Potts–Crowther coefficients gives $\alpha_M^\star \approx 4.4°$ and $(L/D)_{\text{trim}} \approx 3.0$ (less reliable because Potts–Crowther's CM values are reported with reduced precision). The convergence on a self-trimmed glide $L/D$ in the range $2.5$–$3.5$ matches both the in-flight measurements of Lorenz (2005) and the trajectory simulations of Lukes et al. (2022).

**Key analytical statement:** the trim angle $\alpha_M^\star$ (where the pitching moment vanishes) is *not* the same as the angle $\alpha^\star$ that maximizes $L/D$. They differ by $\sim 5–10°$. The disc cannot fly at $\alpha^\star$ in steady state because the precession drives it to $\alpha_M^\star$. The achievable glide $L/D$ is therefore $(L/D)_{\text{trim}} \approx 3.4$, not $(L/D)_{\max} \approx 3.8$. This is a real loss but a small one, and it explains why the wind-tunnel optimum is rarely attained in flight.

---

## 7. Linear stability about the trim glide

### 7.1 Perturbation setup

We linearize the 6-DOF system (18)+(22) about the trim glide of §6.4. The reference state has

$$
V = V_e, \ \alpha = \alpha_M^\star = 6.7°, \ \gamma = \gamma_e = \arctan(C_D/C_L)_{\text{trim}}, \ \Omega = \Omega_e, \ p = q = 0. \tag{40}
$$

Perturbation variables: $(u, w, \delta q, \delta p, \delta\Omega)$ where $u$ is along-path velocity perturbation, $w$ is normal velocity perturbation (so $\alpha = \alpha_e + w/V_e + O(w^2)$), and the angular perturbations are about the trim.

### 7.2 Longitudinal mode (phugoid-like)

The translational equations (18) linearized along and perpendicular to the trim path give, after some standard algebra (Etkin & Reid, *Dynamics of Flight*),

$$
\dot{u} = -g\cos\gamma_e\,(2 u/V_e) - g\sin\gamma_e\,(w/V_e) \cdot \text{(quadratic-drag corrections)}, \tag{41a}
$$

$$
\dot{w} = -g\,(C_{L\alpha}/C_L)\,w/V_e + g\sin\gamma_e\,(2 u/V_e). \tag{41b}
$$

This is the **disc analog of the airplane phugoid**: a slow oscillation in $u$ and $w$ at frequency

$$
\omega_{\text{phugoid}} = \frac{g}{V_e}\sqrt{2}, \qquad
T_{\text{phugoid}} = \frac{\pi V_e \sqrt{2}}{g}. \tag{42}
$$

For $V_e = 12$ m/s, $T \approx 5.4$ s — *comparable to the flight time*. So the phugoid mode is just barely excitable within a single throw; what it produces in practice is a gentle sinusoidal undulation of the trajectory ("S-curve" trajectories of throws against slight tailwinds).

### 7.3 Lateral/precession mode

The transverse equations (22a–b) in the high-spin limit (24) gave algebraic $p, q$. Restoring the inertial terms and linearizing $M, R$ about trim:

$$
I_p \dot{\delta p} + (I_s - I_p)\,\delta q\,\Omega_e = \delta R = -\frac{\partial R}{\partial \alpha}\,\frac{w}{V_e} + C_{Rp}\,\frac{\delta p\, d}{2 V_e}\,Q_{\text{trim}} + \ldots \tag{43a}
$$

$$
I_p \dot{\delta q} - (I_s - I_p)\,\delta p\,\Omega_e = \delta M = C_{M\alpha}\,\frac{w}{V_e}\,Q_{\text{trim}} d + C_{Mq}\,\frac{\delta q\,d}{2V_e}\,Q_{\text{trim}}\,d + \ldots \tag{43b}
$$

with $Q_{\text{trim}} = \tfrac{1}{2}\rho V_e^2 A$. Substituting (43b) into (43a) yields a damped oscillator for $\delta q$ at the **precession frequency**

$$
\omega_{\text{prec}}^2 \approx \Big(\frac{(I_s - I_p)\Omega_e}{I_p}\Big)^2 = \Omega_e^2, \tag{44}
$$

with damping $\zeta_{\text{prec}} \sim |C_{Mq}|\,Q_{\text{trim}} d^2 / (2 I_p V_e \omega_{\text{prec}})$.

For $\Omega_e = 100$ rad/s: $T_{\text{prec}} = 2\pi/\Omega_e \approx 0.063$ s — much faster than the flight. *This* is what makes the disc behave gyroscopically: the precession mode is fast and well-damped, so on flight time scales it is effectively slaved to its quasi-static value (26), which is what we used in §5.

A useful re-derivation: the transverse-precession eigenvalues of the linearized (22a–b) are approximately

$$
\lambda_{\text{prec}} \approx -\frac{|C_{Mq}|\,Q_{\text{trim}} d^2}{4 I_p V_e} \pm i\,\Omega_e, \tag{45}
$$

a pair of slowly damped oscillations at frequency $\Omega_e$. The damping decay time is $\sim 4 I_p V_e / (|C_{Mq}|\,Q_{\text{trim}} d^2)$; for the trim flight, $\sim 30$–$100$ s. The mode is stable but lightly damped — which means *initial-condition errors in p, q at release ring down slowly*. This is exactly the "wobble" experienced after a sloppy release: the disc visibly oscillates for a second or two before settling.

### 7.4 Coupling: precession-modulated phugoid

The cross-coupling between the (slow) translational phugoid mode and the (fast) gyroscopic precession mode generates the turn-and-fade trajectory pattern via the mechanism of §5.2. Formally: the slow drift of $\alpha$ along the phugoid (and along the open-loop trim relaxation) modulates the pitching moment $M(\alpha)$, which through (26b) drives a quasi-static $p \propto -M/\Omega$. The integrated heading change over a flight time $T$ is

$$
\Delta\Psi = \int_0^T p(t)\,dt \approx -\frac{1}{I_p \Omega} \int_0^T M(\alpha(t))\,dt. \tag{46}
$$

If $\alpha(t) - \alpha_M^\star$ averages to zero over the flight (the disc passes from above to below trim), then $\Delta\Psi \approx 0$ — no net heading change. If the disc spends more time above trim (i.e. high-$\alpha$ throws), $\Delta\Psi$ is large and one-signed (anhyzer for a right-handed backhand). This integrated quantity is the practical *throw shape*: how much net curl from release to landing.

---

## 8. Effect of wind

### 8.1 Velocity triangle

The relative wind is $\mathbf{V}_{\text{rel}} = \mathbf{v} - \mathbf{U}$ (eq. 5). For a uniform steady wind $\mathbf{U} = U_0\,\hat{\mathbf{w}}$, the airspeed magnitude and direction relative to ground are shifted:

$$
V = |\mathbf{v} - U_0\hat{\mathbf{w}}|, \qquad \beta_w = \arccos\!\left(\frac{(\mathbf{v} - U_0\hat{\mathbf{w}})\cdot \hat{\mathbf{v}}_{\text{ground}}}{V}\right). \tag{47}
$$

For a tailwind ($\hat{\mathbf{w}}$ parallel to $\mathbf{v}$): $V = v_{\text{ground}} - U_0$. Headwind: $V = v_{\text{ground}} + U_0$. Crosswind: changes effective heading.

### 8.2 Effect on $\alpha$ and the trim point

Since $\alpha$ is set by the geometry of $\mathbf{V}_{\text{rel}}$ relative to the disc plane, a tailwind that reduces $V$ also raises the effective $\alpha$ at constant body orientation (because the descent velocity component $w$ is unchanged but the horizontal component shrinks). Quantitatively, if the disc descends at $\dot{z} < 0$ in steady glide:

$$
\alpha_{\text{tail}} - \alpha_{\text{still}} \approx \frac{w_{\text{glide}}\,(U_0)}{V_{\text{still}}^2}. \tag{48}
$$

For $w_{\text{glide}} = V_{\text{still}} \sin\gamma \approx 4$ m/s, $V_{\text{still}} = 12$ m/s, $U_0 = 5$ m/s tailwind: $\Delta\alpha \approx 4 \cdot 5 / 144 \approx 0.14$ rad $\approx 8°$. A 5 m/s tailwind therefore pushes $\alpha$ from $\sim 7°$ trim to $\sim 15°$ — well above the trim. By (27), this drives a substantial anhyzer turn rate. **Tailwinds make backhands turn more**, a well-known practitioner observation now given an analytical form.

Headwinds do the opposite. A 5 m/s headwind reduces $\alpha$ by $\sim 8°$, pushing the disc below trim and accelerating the fade.

### 8.3 Crosswind drift

For a crosswind $\mathbf{U} = U_c\,\hat{\mathbf{n}}_2$ perpendicular to the intended throw direction, the disc's airspeed direction is rotated by $\arctan(U_c/V_e) \sim 5°/12° = 23°$ for $U_c = 5$ m/s. The trim angle of attack is unchanged (it depends only on $\alpha$ relative to the disc plane), but the disc's groundtrack rotates by the same angle. Cumulatively over a throw, this is a substantial lateral drift — for $T = 4$ s flight: $\Delta y \approx U_c T = 20$ m. **The leading lateral-error term for an outdoor throw is simply $U_c T$.** Any wind-correction scheme starts here.

### 8.4 Wind-gradient effects (the lakefront IBL)

When $\mathbf{U} = \mathbf{U}(z)$ varies with height — as in the Sean Earl Field internal boundary layer, where $|U|$ grows rapidly from the surface up to the IBL top at $z \sim 5$ m — a climbing disc experiences a *changing* relative wind. Differentiating (5),

$$
\frac{d\mathbf{V}_{\text{rel}}}{dt} = \dot{\mathbf{v}} - \frac{d\mathbf{U}}{dt} = \dot{\mathbf{v}} - (\mathbf{v}\cdot\nabla)\mathbf{U}. \tag{49}
$$

The second term acts as an effective extra acceleration. For an updraft-equivalent gradient $\partial U_x/\partial z > 0$ encountered as the disc climbs (vertical velocity $w_z > 0$), the airspeed shifts toward more tailwind — *exactly* the same effect as a sudden tailwind gust at constant altitude, with magnitude proportional to $w_z \partial U / \partial z$. For the lakefront IBL, $\partial U/\partial z \sim 1$ s$^{-1}$ near the surface; a climbing disc at $w_z = 2$ m/s sees a 2 m/s tailwind gust over a one-second climb, which by (48) shifts $\alpha$ by $\sim 3°$. This is enough to bias the trajectory of low backhand throws, and is the analytical justification for caring about the CFD-resolved $U(z)$ field rather than a height-averaged surface wind.

---

## 9. Optimal release conditions

Given the trim analysis, the optimal release for maximum still-air range is:
1. Release at the angle of attack that puts the disc *near but above* the trim $\alpha_M^\star = 6.7°$, so that the disc decelerates through its high-L/D range.
2. Release speed $V_0$ large enough that the disc enters the trim glide *after* the desired range, not before.
3. Spin rate large enough to suppress wobble — practical $\Omega \gtrsim 50$ rad/s, well within typical throw range.

In the presence of a steady wind $\mathbf{U}$:
- Tailwind: release at *lower* $\alpha_0$ than still-air, to counter the wind-induced $\alpha$ increase (48).
- Headwind: release at higher $\alpha_0$ (and lower $V_0$, since induced drag is amplified by the higher airspeed against headwind).
- Crosswind: aim heading at $\arctan(U_c/V_e)$ upwind of the target, plus an iterative correction for the integrated drift $U_c T$.

The two-page summary of practical implications for the paper, written in the throwability-metric language of Project 1, will quantify these statements into a $Q(\mathbf{x}, \theta_{\text{wind}}, \text{throw type})$ surface evaluated over the field.

In the presence of a spatially varying wind field $\mathbf{U}(\mathbf{x})$ from CFD (the central case of the manuscript), the optimal release is found by **trajectory optimization** in the 6-DOF model: a Pontryagin / shooting-method search over $(V_0, \alpha_0, \Omega_0, \theta_{\text{throw}})$ minimizing landing-point error for each target. This is straightforward but numerical; the analytical content is the 6-DOF model itself, and the closed-form approximations of §6–§8 act as the initial conditions for the search.

---

## 10. Summary of analytical results

| Result | Symbol | Closed form | Numerical (Hummel coeffs, $V = 15$ m/s) |
|---|---|---|---|
| Reynolds number | Re | $Vd/\nu$ | $2.7\times10^5$ |
| Spin moment of inertia | $I_s$ | $\tfrac{1}{2}mr^2$ | $1.59\times10^{-3}$ kg·m² |
| Transverse moment of inertia | $I_p$ | $\tfrac{1}{4}mr^2$ | $7.97\times10^{-4}$ kg·m² |
| Precession rate (roll) | $p$ | $-M/(I_p\Omega)$ | quasi-static |
| Precession rate (pitch) | $q$ | $+R/(I_p\Omega)$ | quasi-static |
| Trim angle of attack | $\alpha_{M^\star}$ | $\alpha_0 - C_{M0}/C_{M\alpha}$ | $6.7°$ |
| Trim glide L/D | $(L/D)_{\text{trim}}$ | (39) | $\approx 3.4$ |
| Max-L/D angle | $\alpha^\star$ | (36) | $\approx 17°$ (linearization) |
| Spin-down time | $\tau_{\text{spin}}$ | (17) | $\sim 80$–$400$ s |
| Phugoid period | $T_{\text{phug}}$ | $\pi V_e \sqrt{2}/g$ | $\approx 5.4$ s |
| Precession period | $T_{\text{prec}}$ | $2\pi/\Omega$ | $\approx 0.06$ s |
| Precession damping | $\tau_{\text{prec,decay}}$ | (45) | $\sim 30$–$100$ s |
| Tailwind $\alpha$ shift | $\Delta\alpha$ | $w_g U_0 / V_{\text{still}}^2$ | $+8°$ at $U_0 = 5$ m/s |
| Crosswind drift | $\Delta y$ | $U_c T$ | $20$ m at $U_c = 5$ m/s, $T=4$ s |

---

## 11. Open analytical extensions worth pursuing

1. **Beyond the linear-coefficient model.** The single most-impactful improvement is replacing (11) with the tabulated wind-tunnel $C_L(\alpha), C_D(\alpha), C_M(\alpha)$ from Potts & Crowther (2002). The trim and L/D analyses port over directly with the linearization replaced by table lookups; the trim equation $C_M(\alpha) = 0$ becomes a 1D root-find rather than (38). Worth doing.

2. **Advance-ratio dependence of $C_L, C_D$.** Potts & Crowther's data also shows weak AdvR dependence at AdvR $\gtrsim 0.5$. For low-spin throws (sloppy releases) this matters. A perturbative correction $C_{L} \to C_L (1 + \kappa\,\mathrm{AdvR})$ from the wind-tunnel fits would give a closed-form modifier.

3. **Roisman-style asymptotic analysis of the turn.** The integrated turn angle (46) under a model phugoid trajectory ($\alpha(t)$ given by linearized longitudinal dynamics) should give a closed-form prediction of net hyzer/anhyzer as a function of release $\alpha_0$ and $V_0$. To my knowledge nobody has written this down; it's a tractable 3–4-page calculation.

4. **Non-axisymmetric corrections for the rim profile.** The Ultrastar's rim is slightly non-flat (a "lip"). This introduces a small $C_S(\alpha)$ side force and breaks the body-frame $\hat{\mathbf{b}}_1$-along-airflow assumption used in §3. A perturbative treatment quantifies the residual "side slip" that all practitioners feel but no published analysis treats.

5. **Coupled stability with an inhomogeneous wind $\mathbf{U}(\mathbf{x})$.** The wind-gradient analysis (§8.4) is currently only a single-term WKB-type estimate. A full Floquet treatment of the linearized 6-DOF in a streamwise-periodic wind (representative of a building wake) would predict resonance windows where the trajectory becomes especially sensitive to launch parameters — a publishable analytical result that would naturally accompany the CFD.

Items 1, 3, and 5 are the strongest candidates for inclusion in the manuscript theory section; 2 and 4 are journal-club material.

---

## References (verified)

- Crowther, W.J. & Potts, J.R. (2007). Simulation of a spin-stabilised sports disc. *Sports Engineering* **10**(1), 3–21.
- Etkin, B. & Reid, L.D. (1996). *Dynamics of Flight: Stability and Control.* 3rd ed., Wiley. (Reference for the phugoid mode derivation; the disc analog (41)–(42) follows the same algebra.)
- Hubbard, M. & Hummel, S.A. (2000). Simulation of Frisbee Flight. 5th Math & Computers in Sport, Sydney.
- Hummel, S.A. (2003). *Frisbee Flight Simulation and Throw Biomechanics.* MS Thesis, UC Davis.
- Kamaruddin, N.M., Potts, J.R. & Crowther, W.J. (2018). Aerodynamic performance of flying discs. *Aircraft Eng. Aerosp. Tech.* **90**(2), 390–397.
- Lissaman, P.B.S. & Hubbard, M. (2010). Maximum range of flying discs. *Procedia Engineering* **2**(2), 2529–2535.
- Lorenz, R.D. (2005). Flight and attitude dynamics measurements of an instrumented Frisbee. *Meas. Sci. Technol.* **16**(3), 738.
- Lukes, R., Hart, J. et al. (2022). Disc golf trajectory modelling combining CFD and rigid body dynamics. *Sports Engineering* **25**, 24.
- Potts, J.R. & Crowther, W.J. (2002). Frisbee™ aerodynamics. AIAA 2002-3150.
- Stilley, G.D. & Carstens, D.L. (1972). *Adaptation of the Frisbee Flight Principle to Delivery of Special Ordnance.* Sandia Labs SC-DR-72-0427.

import pathlib

p = pathlib.Path(__file__).parent.parent / "STRESS_AND_DILATANCY_BRIEF.md"
t = p.read_text(encoding="utf-8")
EM = "—"

reps = [
 (r"— a self-similar field whose maximum sits", r"This is a self-similar field whose maximum sits"),
 (r"$6\times10^{-5}$ — roughly **220", r"$6\times10^{-5}$, roughly **220"),
 (r"0.35 at the same point — within a factor of 3 of unity instead of a factor of 100.",
  r"0.35 at the same point, within a factor of 3 of unity instead of a factor of 100."),
 (r"### Stage 0 — the air cushion (before contact)", r"### Stage 0: the air cushion (before contact)"),
 (r"### Stage 1 — the water hammer (nanoseconds)", r"### Stage 1: the water hammer (nanoseconds)"),
 (r"### Stage 2 — Wagner (the one that matters)", r"### Stage 2: Wagner (the one that matters)"),
 (r"### Stage 3 — the lamella", r"### Stage 3: the lamella"),
 (r"**5.4 MPa** — the largest stress anywhere in this problem, 700× the dynamic",
  r"**5.4 MPa**, the largest stress anywhere in this problem and 700× the dynamic"),
 (r"**Step 1 — how big is the contact patch?**", r"**Step 1. How big is the contact patch?**"),
 (r"**Step 2 — replace the drop by a disc.**", r"**Step 2. Replace the drop by a disc.**"),
 (r"**Step 3 — unsteady Bernoulli.**", r"**Step 3. Unsteady Bernoulli.**"),
 (r"**Step 4 — the trick.**", r"**Step 4. The trick.**"),
 (r"prefactor included** — verified", r"prefactor included.** This is verified"),
 (r"- $\rho U_0^2$ — the natural stress scale.", r"- $\rho U_0^2$ sets the stress scale."),
 (r"- $t^{*-1/2}$ — the amplitude.", r"- $t^{*-1/2}$ is the amplitude."),
 (r"- $(1-\xi^2)^{-1/2}$ — **the shape, and this is the stress localization.**",
  r"- $(1-\xi^2)^{-1/2}$ is **the shape, and this is the stress localization.**"),
 (r"bolted onto $t=0$ — **it is what makes", r"bolted onto $t=0$. **It is what makes"),
 (r"piezo sensor — which is exactly what the 2026 paper did, and something your rig does not currently do.",
  r"piezo sensor, which is exactly what the 2026 paper did and something your rig does not currently do."),
 (r"**Treat 5.5 as a soft number** — Philippi", r"**Treat 5.5 as a soft number.** Philippi"),
 (r"inertia — acoustic impedance when waves matter, mass density when they do not.",
  r"inertia: acoustic impedance when waves matter, mass density when they do not."),
 (r"stiffness — it is **the switch", r"stiffness. It is **the switch"),
 (r"$J\lesssim4\times10^{-3}$ — deep in the liquid-like corner.",
  r"$J\lesssim4\times10^{-3}$, deep in the liquid-like corner."),
 (r"| Shah, Arora & Driscoll (2022) — **your own 830 nm particles** |",
  r"| Shah, Arora & Driscoll (2022), **your own 830 nm particles** |"),
 (r"a factor of 11 — and once", r"a factor of 11, and once"),
 (r"divergence — sharper than the", r"divergence, sharper than the"),
 (r"pressure's inverse square root — and **exactly zero on the axis**, by symmetry.",
  r"pressure's inverse square root, and **exactly zero on the axis** by symmetry."),
 (r"**But be honest about what that is worth** — a prediction spanning",
  r"**But be honest about what that is worth.** A prediction spanning"),
 (r"are free of every uncertainty above — viscosity, drop size, criterion, and prefactor all cancel:",
  r"are free of every uncertainty above, because viscosity, drop size, criterion and prefactor all cancel:"),
 (r"That is the right direction — your nominal", r"That is the right direction, because your nominal"),
 (r"not a disc at the centre** — see Figure 4.", r"not a disc at the centre** (Figure 4)."),
 (r"If $\phi_c\approx0.47$ — which is what", r"If $\phi_c\approx0.47$, which is what"),
 (r"on these same 830 nm particles — then at $\phi = 0.44$ your",
  r"on these same 830 nm particles, then at $\phi = 0.44$ your"),
 (r"at $51°$ — a factor of 3.", r"at $51°$, a factor of 3."),
 (r"inflate $F_d$ badly — carried through,", r"inflate $F_d$ badly. Carried through,"),
 (r"- **interface pull-in** — the free surface is dragged",
  r"- **interface pull-in.** The free surface is dragged"),
 (r"$-10$ to $-100$ kPa — **saturating at exactly one atmosphere.**",
  r"$-10$ to $-100$ kPa, **saturating at exactly one atmosphere.**"),
 (r"0.35 there — order unity, off by a factor of 3. The uncapped number is 0.0093 — off by a factor of 100.",
  r"0.35 there, order unity and off by a factor of 3. The uncapped number is 0.0093, off by a factor of 100."),
 (r"somewhere around 30 m/s** — reachable with a pressurised jet or",
  r"somewhere around 30 m/s**, reachable with a pressurised jet or"),
 (r"| (A) $U_0/D_0$ | Shah, Arora & Driscoll 2022 — **your fitted axis** |",
  r"| (A) $U_0/D_0$ | Shah, Arora & Driscoll 2022, **your fitted axis** |"),
 (r"   drops' — the §8 landmine, again.", r"   drops'. This is the §8 landmine again."),
 (r"your system — and if $\phi_c$ really is 0.47,", r"your system. And if $\phi_c$ really is 0.47,"),
 (r"construct — permeability, water viscosity, a friction angle, $\Delta\phi$.",
  r"construct built from permeability, water viscosity, a friction angle and $\Delta\phi$."),
 (r"rather than one driving the other — which is in direct",
  r"rather than one driving the other, which is in direct"),
 (r"**The water hammer is a red herring** — 5.4 MPa, 0.5 ns, one part in a billion of the drop.",
  r"**The water hammer is a red herring.** 5.4 MPa, 0.5 ns, one part in a billion of the drop."),
 (r"> $t^{-1/2}$ — so 'the impact stress' isn't a number,",
  r"> $t^{-1/2}$, so 'the impact stress' isn't a number,"),
 (r"permeability goes as $a^2$ — but the Darcy pressure behind it",
  r"permeability goes as $a^2$. But the Darcy pressure behind it"),
 (r'**If asked "is this published?"** — the stress field yes, thoroughly.',
  r'**If asked "is this published?"** The stress field, yes, thoroughly.'),
 (r'**If asked "how sure are you?"** — Table in §10, and the phrase to use is',
  r'**If asked "how sure are you?"** Point at the table in §10. The phrase to use is'),
 (r"## Appendix A — symbols", r"## Appendix A: symbols"),
 (r"| $\Lambda_F$ | $\langle\tau\rangle/\sigma^\ast$ — proposed shear gauge |",
  r"| $\Lambda_F$ | $\langle\tau\rangle/\sigma^\ast$, the proposed shear gauge |"),
 (r"## Appendix B — sources actually read", r"## Appendix B: sources actually read"),
 (r'"inertio-dilatant" — the latter appears nowhere in the literature.',
  r'"inertio-dilatant", which appears nowhere in the literature.'),
 (r"## Appendix C — how to rerun", r"## Appendix C: how to rerun"),
]

missed = [a[:70] for a, b in reps if a not in t]
for a, b in reps:
    t = t.replace(a, b)
p.write_text(t, encoding="utf-8")
print("MISSED:", missed if missed else "none")

lines = t.split("\n")
prose = [(i + 1, l) for i, l in enumerate(lines) if EM in l and not l.lstrip().startswith("|")]
print("prose lines still holding an em dash:", len(prose))
for i, l in prose:
    print(f"  {i}: {l}")

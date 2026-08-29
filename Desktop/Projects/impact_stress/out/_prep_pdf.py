"""Insert figures into the brief and fix inline-math spacing that pandoc rejects."""
import pathlib
import re

root = pathlib.Path(__file__).parent.parent
p = root / "STRESS_AND_DILATANCY_BRIEF.md"
t = p.read_text(encoding="utf-8")

# --- 1. report malformed inline math (space adjacent to a delimiter) ----------
bad = []
for i, line in enumerate(t.split("\n"), 1):
    if line.startswith("$$"):
        continue
    for m in re.finditer(r"(?<!\$)\$([^$\n]+)\$(?!\$)", line):
        s = m.group(1)
        if s.startswith(" ") or s.endswith(" "):
            bad.append((i, m.group(0)))
print("malformed inline math:", bad)

# --- 2. fix the known ones ---------------------------------------------------
fixes = [
    (r"Peak force $F_{\max} = 0.79\,\rho D_0^2U_0^2 = $ **37 mN**",
     r"Peak force $F_{\max} = 0.79\,\rho D_0^2U_0^2 =$ **37 mN**"),
]
for a, b in fixes:
    t = t.replace(a, b)

# --- 3. insert figures -------------------------------------------------------
figs = [
    ("### Where the compressible stage plugs back in",
     r"![The self-similar wall stress. Left: pressure, shape "
     r"$(1-\xi^2)^{-1/2}$, flat across the patch and divergent at the rim. Right: "
     r"shear stress, shape $\xi/(1-\xi^2)$, exactly zero on the axis and more sharply "
     r"peaked at the rim; the dashed line is $\sigma^*$ for the silica. Each curve is "
     r"one instant. The shapes are identical and only the amplitude changes."
     "](figures/fig1_selfsimilar_fields.png)\n\n"
     "### Where the compressible stage plugs back in"),
    ("**How firm is the rim amplification?**",
     r"![Every stress scale of one impact on one axis, with $\sigma^*$ "
     r"in red for comparison. The water hammer is the largest and the least relevant."
     "](figures/fig2_stress_ladder.png)\n\n"
     "**How firm is the rim amplification?**"),
    ("### The direction that is actually useful: invert it",
     r"![The inertial gauge $\rho U_0^2$ (grey) never approaches "
     r"$\sigma^*$ over the accessible velocity range. The shear gauge "
     r"$\langle\tau\rangle$ (blue, three viscosity choices) crosses it right where the "
     r"jamming onset is measured (green lines)."
     "](figures/fig3_gauge_comparison.png)\n\n"
     "### The direction that is actually useful: invert it"),
    ("---\n\n## 6. The inertial–dilatant number",
     r"![Where in the contact patch the suspension is above $\sigma^*$ "
     r"(shaded). The frictional region is an annulus at the rim that fills inward as the "
     r"impact proceeds, and it appears earlier and reaches further in at higher speed."
     "](figures/fig4_frictional_annulus.png)\n\n"
     "---\n\n## 6. The inertial–dilatant number"),
    ("### Sensitivities",
     r"![ID against particle size at fixed $\Delta\phi$. Solid: Eq. 6 as "
     r"published, $\mathrm{ID}\propto a^2$. Dashed: capped at the pressure the pore fluid "
     r"can sustain, which removes the size dependence entirely. The two vertical lines are "
     r"your silica and the paper's cornstarch."
     "](figures/fig5_ID_vs_particle_size.png)\n\n"
     "### Sensitivities"),
]
for anchor, block in figs:
    if "](figures/" in t and anchor in block and block.split("](figures/")[1][:8] in t:
        continue  # already inserted
    assert anchor in t, anchor[:60]
    t = t.replace(anchor, block, 1)

p.write_text(t, encoding="utf-8")
print("done; figures referenced:", t.count("](figures/"))

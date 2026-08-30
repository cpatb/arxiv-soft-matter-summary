# calibration_transfer.png

Regenerate: `python calibration/calibrate.py && python calibration/report.py`.
Full argument: `../../CALIBRATION_2026-08-30.md`.

## Left panel — a flow curve does not report the packing of the sample

Vertical axis is the **residual**, φ inferred minus φ made, not φ inferred
against φ made. Plotted that way because the house rule keeps zero on a linear
axis, and on the 0-to-0.65 axis the raw comparison needs, a 0.04 offset is
invisible. Zero, the dashed line, is where every point would sit if a
Wyart-Cates fit returned the volume fraction the sample was weighed out at.

- **Filled circles**: the free three-parameter fit, on curves the *source lab*
  calls steady.
- **Open circles**: the free fit on that lab's own unsteady curves (Guy et al.,
  φ ≥ 0.56, temporal averages with wall slip; Richards et al., the three their
  deposited plotting script colours red).
- **Open triangles**: the fit-free plateau inversion,
  φ = φ₀(1 − η_r^(−1/n)) from the lowest-stress viscosity. No optimiser, no
  starting guess, no branch cut. Absent for the 404 nm PMMA because that series
  is Brownian at every stress it samples (Pe = 2×10⁻⁴ to 3×10⁻³), so the
  frictionless branch it inverts is not what was measured.
- **Lines**: the free-fit regression over all usable curves of that material.

The residual rises with φ in all three materials and crosses zero near φ ≈ 0.5
to 0.58. The error is largest where the suspension is most dilute and vanishes
as the packing approaches jamming — which is the regime where the curve carries
information about the distance to jamming and nothing else.

## Right panel — φ_m from one material's own curves

Each bar spans the smallest to largest φ_m the free fit returns across the
curves of a single material in a single paper. **Black tick** is that paper's
own published φ_m where it quotes one (Guy 0.558(5), Richards 0.457).

Two bars per material: **all usable curves** (solid) and **the source lab's own
steady subset** (pale). The restriction takes Guy's 3770 nm from 0.058 wide to
0.013 and corn starch from 0.060 to 0.012 — the dramatic version of the spread
is largely an artefact of curves those labs already disowned, which is what
`dst_fit`'s gate N8 found and how it must be quoted.

What survives is still load-bearing: 0.012 to 0.013 is six or seven times the
±0.002 tolerance a suspension's volume fraction is actually set to, and jamming
is decided by sign(φ − φ_m), so the verdict flips right across the remaining
band. It flips for all three materials on their steady subsets.

## Caveats a reader should carry

- Corn starch is in **weight** fractions, as its authors use, because the
  granules are porous. Its φ axis is not comparable with the PMMA ones.
- The 404 nm series contributes no plateau points and its bootstrap intervals
  contain the identity. It is drawn because its exclusion is itself the finding
  of §3 of the brief, not because it supports the map.
- The regression lines are drawn over all usable curves; the steady-subset
  regressions have roughly triple the interval width and are in the brief's
  table rather than on the plot.

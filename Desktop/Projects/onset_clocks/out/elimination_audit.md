# elimination_audit.png

Regenerate: `python run_audit.py`. Full argument: `../README.md`.

## Left panel — the contrast that carries the result

Vertical axis is the measured arrested footprint `beta_p`, the thing every
candidate has to explain. It spans 0.309 to 1.962, a factor of 6.35.

Two candidates, each on its own x-axis because they live on incomparable
scales — the point is the shape, not a common abscissa.

- **Green circles, bottom axis: the jamming deficit** `(phi0-phi)/(phi0-phiM)`,
  the strongest survivor of the audit (Spearman rho = +0.708, p = 2.5e-12). A
  pure packing group: no stress, no rate, no strain, no time in it.
- **Red squares, top axis: the accumulated strain over one impact time on the
  field's standard shear rate** `gammadot = U0/D0`. Its dynamic range over the
  dataset is 1.00x — the quantity is *identically 1* for every drop by
  construction, since the strain is `(U0/D0)(D0/U0)`. The scatter visible on
  that axis is float rounding, nothing else.

The red column is the finding. A criterion that does not move cannot explain a
response that moves by a factor of six. The stress clocks are only slightly
better: `f` at either gauge spans 1 to 4%, which on this axis would be a column
of similar width.

## Right panel — the audit

One bar per candidate, spanning 1 to its own dynamic range over the 73 drops.
Log axis, because the candidates span more than a decade between them.

- **Dashed line: 6.35x**, the response's own range. It is a reference, not a
  threshold — a steep enough function could in principle map a smaller range
  onto a larger one, which is why the verdict comes from the degeneracy test
  rather than from this line.
- **Red: FALSIFIED.** Drops whose group agrees to 5% retain more than half the
  typical spread in `beta_p`, so knowing the group has told you nothing.
- **Green: survives.**

Note that dynamic range and verdict are *not* the same ordering: `Sigma/sigma*`
and `We` both have large ranges and are still falsified, because their ties are
degenerate. That is the point of using the degeneracy test rather than the range.

## Caveats a reader should carry

- The strongest survivor is a monotone relabelling of `phi`, one of the two
  knobs the campaign turned. It orders the packing axis and carries no velocity
  or salt dependence at all, so surviving is not the same as explaining.
- Every material constant comes from one calibration
  (`Github/Drop Impact/closure/params_fitted.yaml`) that failed its own
  acceptance gate.
- The Mobaseri lamella clocks are undefined for the lowest-velocity,
  highest-`phi` drops, where Re falls to order 1; those drops are excluded from
  those two bars only, and the counts differ between bars for that reason.

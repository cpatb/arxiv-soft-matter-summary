# rehab_smartfluid progress ledger

Newest first. Records what failed and why, not only what worked.

---

## 2026-09-01: layer 3, and closed forms behind the layer-2 map

**Done**

- `derivations/shortfall_asymptotics.py` + `SHORTFALL_ASYMPTOTICS.md`. The
  quasi-steady shortfall now has closed forms, verified symbolically and
  numerically at four packings: the `A -> 0` plateau is exactly
  `eta(0)/eta(f_s) = sigma_lub/sigma_steady`, matching the integrator to 1e-6;
  the initial slope matches to 1e-7; the crossover `A*` inverts in closed form.
- `layer3/semiactive.py` + 6 tests. Semi-active plant with the field acting on
  `sigma*`. Three properties established: no instantaneous force authority
  (relative degree one), bandwidth `omega_c = gamma_dot/gamma_0` that rises with
  speed, and **zero authority at rest**.
- `CONTROL_2026-09-01.md` + `scripts/run_control_study.py` + figure. The bandwidth
  result, the closed-loop tracking study, and what both mean for the product.
- 19 tests pass from the repo root.

**Changed a conclusion**

- **`SCALING` §4's "quasi-steady is safe above A ~ 10" is superseded.** The closed
  form gives `A* = 4.59, 6.91, 7.49, 7.73` at `phi = 0.40, 0.46, 0.50, 0.52`. It
  rises with packing fraction, so there is no single threshold, and "10" was a
  guess at the top of a range read off a table.
- **`DESIGN` §11's "STF or MR first" is answered, against the framing of the
  question.** Actuator bandwidth spans 0.99 to 25.5 Hz over the whole published
  design space, against 10 Hz for weight acceptance and 33 Hz for a stumble. Over
  most of the space the controller is slower than the event. **Go STF-passive
  first**; semi-active earns its cost only in the high-rate short-clock corner,
  which is a material-selection target rather than a control-design one.
- **`A` is the same number twice.** It counts the material's strain budget and the
  actuator's time constants in an event. Layers 2 and 3 fail together.
- **The closed loop needs about 10x more strain than the material does.** Passive
  shortfall is gone by `A* = 4.6-7.7`; closed-loop tracking error is still 8.6% at
  `phi = 0.40` and 22% at `phi = 0.50` when `A = 13`, and needs `A ~ 44` for 5%.
  Not expected. The passive crossover asks whether a step reaches its peak; the
  controller lags a *moving* demand, so its error survives the material settling.

**Failed, and why**

- **arXiv fetches time out.** Gillissen et al. (1911.00280), which
  `FIELD_SURVEY` names as the model to port in phase 1, could not be retrieved.
  Two attempts. Still queued.
- **A heredoc patch silently did not apply and I relaunched the slow script.**
  `\\n` inside a `<<'PY'` heredoc becomes a real newline in the Python source, so
  the match against a file containing a literal backslash-n failed, the assertion
  fired, and because the heredoc and the run command were separate lines rather
  than chained, the *unpatched* script ran again. This is the exact trap already
  in `lessons.md`: the quoting protects the heredoc from bash, not from Python.
  Fixed by using the Edit tool.
- **The control study was too slow to run interactively**, twice. The demand was
  re-solving a fixed point at every ODE right-hand side. Now tabulated on a
  velocity grid, with a grid-doubling adequacy check printed beside the result.

**Open**

- Gillissen port still blocked on fetching the paper.
- The 30 ms stumble duration is still unsourced, and the stumble bandwidth verdict
  (1.31x short at the best corner) is marginal enough that this number decides it.
- `gamma_0` for silica/PEG specifically is still the single measurement that would
  most change the conclusions.

---

## 2026-08-31: phase 0 opens

**Done**

- `layer2/activation.py` + 12 tests, reusing `dst_fit` for Wyart-Cates rather than
  restating it, and carrying its gauge lesson (`phi_0` is not a design knob).
- `SCALING_2026-08-31.md`: `A = gamma_dot tau / gamma_0` is order one in gait
  (0.63-16 walking, 0.19-4.8 stumble). Shortfall amplifies with packing: at
  `A = 0.3`, 79% of design stress at `phi = 0.40` and 5.9% at `phi = 0.52`.
  Reachable shear rate collapses through the gait band at `phi_m`.
- `FIELD_SURVEY_2026-08-31.md`: two literatures that do not cite each other.

**Failed, and why**

- **A fixture in a dead regime nearly hid everything.** The first test material
  used a water-like carrier, putting `sigma` at 1.4 Pa against `sigma* = 50`, so
  `f_ss = 6e-16` and nothing thickened. Every test passed vacuously except the
  discriminating negative. Fixed by *solving* for the carrier viscosity that puts
  the onset at the reference rate; the answer, 0.07 Pa.s, is PEG-like, which is
  what the published devices use. A fixture guard now prevents the drift.
- **The S-shape test probed rates the material cannot reach**, 1e2 to 1e6 1/s on a
  packing whose steady curve turns over at 4.14 1/s.

**Changed the plan**

- The survey found Gillissen, Ness, Peterson, Wilson and Cates already published
  the tensorial time-dependent model. **The constitutive model is not the novel
  part**; the design map and the application are.

# A real high-concentration monoclonal antibody — the venture's wedge fluid, at last

Downloaded 2026-09-08. **CC BY 4.0.** This is the third licensable dataset in the project and
the first one that is the fluid the venture is actually about. `OFFER.md` §6 names it as
capability gap 3: *"No protein or high-concentration formulation curve, so nothing here has
been demonstrated on the wedge's own fluid. The examples are human blood and a microgel."*

| field | value |
|---|---|
| Zenodo record | 14409284 |
| DOI | `10.5281/zenodo.14409284` |
| licence | **CC BY 4.0** — commercial use permitted with attribution |
| deposit creators | Prašnikar, Bjelošević Žiberna, Križišnik, Roškar, Grabnar, Žula, Ahlin Grabnar |
| source paper | Prašnikar et al., *"Additive effects of the new viscosity-reducing and stabilizing excipients for monoclonal antibody formulation"*, Int. J. Pharm. **674** (2025) 125451, `10.1016/j.ijpharm.2025.125451` |
| paper access | **open access**; the accepted PDF is in the University of Ljubljana repository (`repozitorij.uni-lj.si`, record 167794) |
| material | one model monoclonal antibody, 58 formulations with viscosity-reducing excipients |
| instrument | **RheoSense VROC / m-VROC (Initium One Plus)**, microfluidic MEMS viscometer-rheometer |
| conditions | 25 °C |

| file | bytes | md5 |
|---|---|---|
| `raw/Table_S1.xlsx` | 89195 | `78a8fb89ed151b1c4ad5612b482854ff` |

The `Measurements` sheet carries **58 formulations**: label, excipient and its concentration,
**mAb concentration 159.3–212.2 mg/mL**, **average viscosity 13.57–47.07 mPa·s** with a
standard deviation, a viscosity-reduction percentage with a confidence interval, and SEC
stability columns. The `Test compound properties` sheet carries the excipients' chemistry.

## The two sentences from the Methods that decide what may be done

Both are quoted verbatim from §2.2.4 of the open-access PDF, not paraphrased:

> Viscosity was measured using a viscometer-rheometer (RheoSense VROC (Viscometer/
> Rheometer-On-a-Chip) or m–VROC, Initium One Plus, USA) based on microfluidic and MEMS
> (Micro-Electro-Mechanical Systems) technologies. **First, we measured the viscosity of
> highly concentrated protein solutions at shear rates of 960–3329 s⁻¹ … categorized them as
> Newtonian. Consequently, samples were measured at 25 °C and shear rates within the
> described range to ensure a constant flow rate.**

> The viscosity reduction was calculated as a relative decrease in viscosity (%) of the
> sample containing the test compound compared to the starting formulation without the test
> compound at the same mAb concentration, which was calculated using the exponential
> calibration curve.

**1. The fluid is Newtonian by the source's own test over a stated range, not by
assumption.** That licenses a one-parameter model and it names the range the model was
verified over. Everything downstream that says "Newtonian" points here.

**2. That range is needle-adjacent, and it is why this dataset matters more than the other
two licensable ones.** 960–3329 1/s, against a wall rate of 1.8e4 1/s in the case study's
needle: **0.74 decades of extrapolation**. The DMEM records stop at 1 1/s (four to five
decades) and the Carbopol-in-propylene-glycol records at 1e3 1/s. The blood donors
extrapolate 1.42 decades. **This is the least-extrapolated fluid in the project.**

## Three things about this deposit that decide what is NOT done

**1. The `viscosity reduction [%]` column does not reproduce, and is not used.** The Methods
define it as `100·(1 − η/η_base(c))` with `η_base` an exponential (Mooney) calibration curve.
Back-solving `η_base = η/(1 − red/100)` and fitting one exponential gives **R² = 0.4575** on
the 49 rows of the main series — a single calibration curve would collapse that to numerical
precision. Either more than one curve was used or the column carries processing the Methods
do not describe. The case study uses the **measured viscosities** and never this column.
What does survive is the **slope**: 0.03131 mL/mg on the main series and 0.03233 on the
higher-concentration one, agreeing to 3.2 %, so viscosity doubles every **21.8 mg/mL**.

**2. There are two series, and the labels are the only thing that says so.** Nine rows carry
a `200 mAb` suffix and sit at 185.5–212.2 mg/mL; the other 49 sit at 159.3–184.2. Pooling
them is what takes the Mooney back-solve from R² = 0.945 (series B alone) to 0.458.

**3. No record from this deposit can enter the registry, and it is a property of the schema
rather than of the data.** A registry record needs a fit and an identifiability sidecar
written by `biorheo.identify`. A Newtonian model has one parameter and each formulation here
has **one** viscosity: `fit_rate_curve` raises `ValueError: 1 points cannot determine 1
parameters`. That is demonstrated in the case study's G25.7 rather than asserted.

**This is worth stating plainly rather than treated as a shortcoming, because it is also
what a client's data usually looks like** — a table of formulation viscosities from a
plate-based or microfluidic screen, not a flow curve. `case_studies/mab_injectability/`
therefore reads the deposited table directly, which is the arriving-data path a paid study
would actually take.

## What is deliberately not touched

The SEC monomer / aggregate / accelerated-stability columns and the DLS `kD` and `B22`
columns are real data and are **out of scope**: stability is not rheology and this package
has nothing to say about it. An excipient that fixes viscosity and destabilises the protein
is not a win, and that trade-off is in the source for a reader who needs it.

## Two sibling deposits from the same group, examined and not used

- **Zenodo 10371967** (CC BY 4.0), *"The Search for Novel Proline Analogs…"* — the earlier
  screen. Relative viscosity of a model mAb with 25 mM of each test compound, one viscosity
  per compound. Same schema boundary, less concentration range.
- **Zenodo 15790919** (CC BY 4.0) — a **literature review** table of viscosity-reduction
  percentages from other people's papers, with molecular descriptors. Second-hand numbers
  without their curves; this repo's standing rule excludes it.

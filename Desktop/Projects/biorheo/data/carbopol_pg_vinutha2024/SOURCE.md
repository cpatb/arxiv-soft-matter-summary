# Carbopol 974p in propylene glycol (Vinutha et al. 2024) — a LICENSABLE yield-stress fluid at three temperatures

Downloaded 2026-09-08. **CC BY 4.0.** This is the second commercially licensable dataset in
the project, after `dmem_sints2024/` (NEXT_STEPS item 22), and it is the one that matters
for the wedge: the DMEM records are a shear-thinning medium whose rates stop at 1 1/s, while
this is a **yield-stress microgel measured to 1e3 1/s** — the same family as the 15
Carbopol records the injectability case study ships, all of which are CC BY-NC-SA and
therefore cannot appear in a public demo.

| field | value |
|---|---|
| Zenodo record | 11389467 |
| DOI | `10.5281/zenodo.11389467` |
| licence | **CC BY 4.0** — commercial use permitted with attribution |
| deposit creators | Del Gado, Emanuela; Trappe, Veronique; Vinutha, H. A. |
| published | 2024-05-29 |
| source paper | Vinutha, Marchand, Caggioni, Vasisht, Del Gado, Trappe, *"Memory of shear flow in soft jammed materials"*, arXiv:2312.00251 |
| material | 2 wt% Carbopol 974p (Lubrizol) in propylene glycol; polyacrylic acid microgel, mean radius ~10 um, estimated volume fraction ~74 % |
| instrument | Anton Paar MCR 300, **stress-controlled**, cone-plate, cone angle 1 deg, **both surfaces roughened to minimise wall slip**, Peltier plate + Peltier hood, T to 0.1 deg C |
| temperatures | 10, 20, 50 deg C |

## Why this dataset and not another

Four CC BY 4.0 candidates were listed in NEXT_STEPS item 20. This one was taken first
because of a property none of the others has: **the temperature enters the physics through
one channel the authors measured separately.** They chose propylene glycol precisely
because its viscosity moves by about a factor of ten over 10-50 deg C, which lets them vary
the viscous stress independently of the shear rate. So the material is the same at all three
temperatures and only the suspending medium changes.

That makes it the first real test bed for `constitutive/shift.py` (item 21). Until now the
one-parameter time-temperature shift has been exercised on synthetic pairs and on the DMEM
records, where the collapse was **pre-registered to fail and did** (G22.6). Here a collapse
is expected to succeed, and — the part that makes it an anchor rather than a demo — **the
shift factor is published, not fitted**: see the two numbers below.

## The numbers this deposit publishes that biorheo does not compute

Both are carried in `raw/Fig 1(c) and Fig S3 residual stress_exp.txt`, in named columns, and
both are corroborated by the paper's Methods section.

| quantity | 10 deg C | 20 deg C | 50 deg C |
|---|---|---|---|
| `medium viscosity [Pa s]` (propylene glycol) | 0.102 | 0.055 | 0.012 |
| `yield stress [Pa]` | 8.9 | 8.9 | 8.9 |

Two things follow.

1. **The yield stress is quoted as ONE number for all three temperatures.** The paper states
   it directly: *"The yield stress is independent of the solvent viscosity, which denotes
   that sigma_y is exclusively set by the intrinsically solid characteristics of the dense
   ensemble of soft spheres."* Item 21's derivation says a one-parameter shift leaves the
   yield stress invariant. Those are the same statement, reached independently.

2. **The shift factor is therefore predicted, not fitted.** If the temperature acts only
   through the medium viscosity, `a = eta_s(T) / eta_s(T_ref)`, which is exactly
   `shift.shift_from_two_viscosities`. `shift.collapse_shift` fits `a` from the flow curves
   alone and never sees `eta_s`. Comparing the two is a check of a FITTED quantity against a
   number from outside this repo.

**Read the caveat before using the 8.9 Pa.** It is the authors' fit to the *three-component
model* of ref. [18] (`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`), which biorheo does not
implement. biorheo fits Herschel-Bulkley. A yield stress is a fit-objective convention
before it is a material constant — that is a recorded lesson in this project, from Horner's
Casson tau_y — so a Herschel-Bulkley `tau_y` is NOT expected to equal 8.9 Pa, and a gate
that demanded it would be comparing two model forms and calling the difference an error.

## Files

`raw/` holds the deposit verbatim. All five experimental files are committed; the two
simulation-only families of `.dat` files were not downloaded, because this repo has no use
for a molecular-dynamics trajectory.

| file | bytes | md5 |
|---|---|---|
| `Fig 1(a) flow curves_exp.txt` | 3489 | `cac7d83d62a6bf5d20fb90b956199a3e` |
| `Fig 1(c) and Fig S3 residual stress_exp.txt` | 2677 | `b5d2d8dba6fde52599964f015144aafc` |
| `Fig S1 shear modulus_exp.txt` | 1471 | `2150a32ce075204092a03b55dd94669f` |
| `Fig 1(b) and Fig S3 selected flow cessation tests_exp.txt` | 518956 | `8bb9a7f68c8992ddf8bdea3e24b368d0` |
| `All flow cessation tests_exp.txt` | 1778408 | `a6dee4ad7296df4385081825f6fe9278` |

Encoding is **cp1252, CRLF** — the degree sign in every `temperature [°C]` header is a
single 0xB0 byte and a UTF-8 read raises. Every file is a **wide** layout: several
independent blocks laid side by side across the columns, separated by an empty column, with
the header repeated in each block. An importer that reads the first three columns gets one
temperature and silently discards the other two, which is the same shape as the trap
recorded for `dmem_sints2024/` and for the Horner workbook's concatenated ramps.

## The four traps in this deposit

Numbers here are printed by `scripts/carbopol_pg_candidate_scan.py`; see
`candidate_scan.md`. They are repeated in prose only because they decide what may become a
record.

1. **`Fig 1(a)` is three blocks, not one curve**, at 10 / 20 / 50 deg C with 48 / 67 / 60
   points, and the blocks have different lengths, so pairing the rate and stress columns by
   filtering each for numeric cells independently misaligns them. Pair row by row.

2. **The 20 deg C block is not instrument output and must not become a record.** Four
   independent signals agree, and no one of them would be enough alone:
   - its rates do not sit on the programmed logarithmic grid (3 of 67, against 37 of 48 at
     10 deg C and 50 of 60 at 50 deg C);
   - its rate **goes backwards twice** inside a monotonically decreasing ramp
     (32.4 -> 33.6 and 24.4 -> 26.4 1/s), which a programmed sweep cannot do;
   - its stresses are rounded to **3 significant figures on 62 of 67 points**, against 7-8
     at 10 deg C and 4-5 at 50 deg C;
   - it is the only block whose independent replicate (trap 3) is offset **one-sidedly**,
     reading 6-14 % high on 4 of 4 matched rates.

   Whether it was digitised from the figure or measured under a different protocol cannot
   be settled from the deposit, and it does not need to be: it cannot be certified as a raw
   curve at the grade of the other two, and this repo does not add a record it cannot
   certify. **The exclusion costs the reference temperature the medium viscosity is quoted
   at (0.055 Pa s), and it does not cost the shift test**, which runs 10 -> 50 deg C, the
   pair with the largest viscosity contrast (8.5x) and the two clean curves.

3. **There is a SECOND flow curve in the deposit, and nothing says so.** The two flow
   cessation files carry a `preshear-stress [Pa]` column at 16 preshear rates
   (1e-3 to 100 1/s) at each of the three temperatures. That is a steady flow curve of the
   same material on the same instrument under a different protocol, i.e. an independent
   replicate — the only one this project has for any registry record. At matched rates it
   agrees with `Fig 1(a)` to a median of **2.5 %** at 10 deg C and **2.9 %** at 50 deg C
   (worst 6.6 % and 11.6 %). That spread is a measurement of the material and the protocol,
   and it is the scale a fitted parameter's confidence interval should be compared against.

4. **`Fig S1` is not this sample.** It is shear modulus against Carbopol concentration over
   1.7-2.6 wt%, used to estimate the volume fraction, and it carries a second block of
   reference data from Mason et al. PRL 1995. It is not a flow curve and no record comes
   from it.

## What this dataset cannot do

It does not reach needle rates. The flow curves stop at 1e3 1/s against 1e4-1e5 at the wall
of a fine needle — one to two decades of extrapolation. That is far better than the DMEM
records' four to five, and far worse than the 1.42 decades measured for the Horner blood
donors, and any injectability number taken from these records has to say so.

The suspending medium is propylene glycol, not water, and the material is a model system for
a rheology study, not a formulation anyone injects. It is a **licensable stand-in for the
Carbopol records**, not a biofluid.

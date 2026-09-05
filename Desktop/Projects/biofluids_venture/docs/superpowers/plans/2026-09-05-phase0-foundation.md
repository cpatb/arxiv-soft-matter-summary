# Phase 0 Foundation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `biorheo` monorepo in which a cited fluid record loads from a registry, a rate-space constitutive model with declared regularisation is identified from data by the existing `dst_fit` machinery with a per-parameter verdict, and one closed-form injectability prediction carries uncertainty from that fit and matches a hand-computed anchor.

**Architecture:** Three packages under one `src/biorheo/` namespace: `registry` (YAML fluid records with provenance validation), `constitutive` (rate-space GNF and structural models expressed as `dst_fit.models.Model` instances so the existing fitter and profiler work unchanged), `identify` (thin layer over `dst_fit` that turns profile-likelihood intervals into verdicts), plus `apps/injectability.py`. `dst_fit` stays where it is and is consumed as an editable sibling dependency; it is not moved.

**Tech Stack:** Python 3.12, numpy 2, scipy 1.17, sympy, pint, PyYAML, pytest. `dst-fit` 0.1.0 editable from `Desktop/Projects/dst_fit`.

## Global Constraints

- Every number in a record or a test has a source (DOI, dataset DOI, or a hand computation shown in the test).
- No `si_fit` anywhere; no "illustrative" parameters in `registry/records/`.
- Regularisation knobs are mandatory keyword arguments with no defaults on anything a solver export could consume.
- Out-of-envelope evaluation raises; missing provenance raises at load.
- A fit that hits a bound returns an interval and the verdict `bounded_only`, never a point estimate.
- Gates are pytest tests that exit non-zero; each unit's gate is listed in `05_ROADMAP.md`.
- Files written by scripts use `newline="\n"` and UTF-8.
- Windows paths are never handled with sed/perl.

## Deviation from the spec, recorded

Spec unit 0.1 said "move `dst_fit` into `packages/identify`". This plan keeps `dst_fit` at `Desktop/Projects/dst_fit` and depends on it. Reasons: (1) the package name (D1) is unresolved, so a physical move now would be moved again; (2) the Obsidian dst-fit card and `_validate.py` cite the current path; (3) the gate ("tests green, byte-identical figures") is satisfied by an unchanged package more cheaply than by a moved one. The move is one `git mv` plus a path edit once D1 is decided. `fit_model`, `profile_likelihood` and `Model` are generic in the independent variable (they only require `func(x, *params, **fixed)`), so rate-space models are `Model` instances with `x = gammadot`; nothing in `dst_fit` changes.

## File structure

```
Desktop/Projects/biorheo/
  pyproject.toml
  README.md
  src/biorheo/__init__.py
  src/biorheo/units.py                 pint registry, one instance
  src/biorheo/registry/__init__.py
  src/biorheo/registry/schema.py       FluidRecord, Envelope, Source, Parameter dataclasses
  src/biorheo/registry/provenance.py   DOI syntax check + resolution with on-disk cache
  src/biorheo/registry/loader.py       load_record(path) -> FluidRecord, raises on defects
  src/biorheo/registry/records/*.yaml
  src/biorheo/registry/doi_cache.json
  src/biorheo/constitutive/__init__.py
  src/biorheo/constitutive/base.py     RateModel (dst_fit Model + control tag + guess + limits)
  src/biorheo/constitutive/gnf.py      power_law, carreau_yasuda, cross, herschel_bulkley
  src/biorheo/constitutive/regularise.py  papanastasiou, biviscosity, clamp; induced_error
  src/biorheo/constitutive/structural.py  Moore-Coussot lambda kinetics (unit 0.5)
  src/biorheo/identify/__init__.py
  src/biorheo/identify/fit.py          fit_rate_curve -> IdentifiedModel with verdicts
  src/biorheo/identify/synthetic.py    synthetic truth + recovery harness
  src/biorheo/apps/__init__.py
  src/biorheo/apps/injectability.py    Buckingham-Reiner, power-law pipe flow, force with UQ
  derivations/                          derive-and-check outputs
  tests/test_units.py
  tests/test_registry_schema.py
  tests/test_provenance.py
  tests/test_loader.py
  tests/test_gnf.py
  tests/test_regularise.py
  tests/test_identify.py
  tests/test_structural.py
  tests/test_injectability.py
  tests/test_gate_0_1_dstfit_unchanged.py
```

---

### Task 1: Monorepo skeleton and the 0.1 gate

**Files:**
- Create: `pyproject.toml`, `README.md`, `src/biorheo/__init__.py`, `src/biorheo/units.py`
- Test: `tests/test_units.py`, `tests/test_gate_0_1_dstfit_unchanged.py`

**Interfaces:**
- Produces: `biorheo.units.ureg` (pint `UnitRegistry`), `biorheo.units.Q_`; `biorheo.__version__`.

- [ ] **Step 1: Write the gate test**

```python
# tests/test_gate_0_1_dstfit_unchanged.py
"""Unit 0.1 gate: dst_fit is consumed unchanged. Its suite is green from here and
its one figure regenerates byte-identically (hash recorded 2026-09-05)."""
import hashlib, subprocess, sys
from pathlib import Path
import pytest

DST_FIT = Path(r"C:\Users\compl\Desktop\Projects\dst_fit")
FIG = DST_FIT / "figures" / "fig_phim_drift.png"
FIG_MD5 = "<fill from Step 2>"

def test_dst_fit_importable_at_pinned_version():
    import dst_fit
    assert dst_fit.__version__ == "0.1.0"

@pytest.mark.slow
def test_dst_fit_suite_green():
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x"], cwd=DST_FIT,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:]

def test_figure_hash_unchanged():
    assert hashlib.md5(FIG.read_bytes()).hexdigest() == FIG_MD5
```

- [ ] **Step 2: Record the figure hash**

Run: `python -c "import hashlib;print(hashlib.md5(open(r'C:\Users\compl\Desktop\Projects\dst_fit\figures\fig_phim_drift.png','rb').read()).hexdigest())"` and paste into `FIG_MD5`.

- [ ] **Step 3: pyproject and units**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"
[project]
name = "biorheo"
version = "0.0.1"
requires-python = ">=3.12"
dependencies = ["numpy>=2.0", "scipy>=1.14", "sympy>=1.12", "pint>=0.24", "pyyaml>=6.0", "dst-fit>=0.1.0"]
[project.optional-dependencies]
dev = ["pytest>=8"]
[tool.setuptools.packages.find]
where = ["src"]
[tool.setuptools.package-data]
"biorheo.registry" = ["records/*.yaml", "doi_cache.json"]
[tool.pytest.ini_options]
testpaths = ["tests"]
pythonpath = ["src"]
markers = ["slow: runs another package's suite", "network: needs internet"]
```

```python
# src/biorheo/units.py
"""One pint registry for the whole package. Boundaries convert once; internals are SI floats."""
import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity
def to_si(value, unit: str) -> float:
    """Convert a (value, unit string) pair to the SI base magnitude as a float."""
    return float(Q_(value, unit).to_base_units().magnitude)
```

- [ ] **Step 4: Run and commit** — `pytest tests/test_units.py tests/test_gate_0_1_dstfit_unchanged.py -v` → PASS.

---

### Task 2: Registry schema and provenance (unit 0.2)

**Files:**
- Create: `src/biorheo/registry/schema.py`, `provenance.py`, `loader.py`, `doi_cache.json`
- Test: `tests/test_registry_schema.py`, `tests/test_provenance.py`, `tests/test_loader.py`

**Interfaces:**
- Produces: `FluidRecord(id, fluid, model, parameters: dict[str, Parameter], envelope: Envelope, sources: list[Source], licence, identifiability: dict|None, notes)`; `Parameter(value, unit, si, stderr|None)`; `Envelope(gammadot_min, gammadot_max, temperature_C, geometry)`; `Source(doi, role, title)`; `load_record(path) -> FluidRecord`; `resolve_doi(doi, cache) -> ResolvedDOI|None`; `ProvenanceError`, `EnvelopeError`.

Rules enforced by `load_record`:
1. `sources` non-empty and every `doi` matches `^10\.\d{4,9}/\S+$` and resolves (cache first, then Crossref `https://api.crossref.org/works/{doi}`; a record whose DOI is not in the cache and cannot be resolved raises `ProvenanceError`).
2. Every parameter carries `unit`; `si` is computed by `to_si` at load, never stored.
3. `envelope.gammadot_min < gammadot_max`, both > 0.
4. `identifiability` is either absent or was written by `biorheo.identify` (carries `written_by: biorheo.identify` and a `fit_hash`); a hand-written block raises.
5. `licence` present; `commercial_use: false` records are loadable but `FluidRecord.commercial_ok` is False.

Test cases: synthetic record with `doi: 10.0000/not-a-real-doi` → `ProvenanceError`; record with cached DOI loads; parameter without unit → `ValueError`; envelope inverted → `EnvelopeError`; hand-written identifiability → `ProvenanceError`; `FluidRecord.eta_at(gammadot)` outside envelope → `EnvelopeError`.

---

### Task 3: GNF family with regularisations (unit 0.4)

**Files:**
- Create: `src/biorheo/constitutive/base.py`, `gnf.py`, `regularise.py`, `derivations/gnf_limits.md`
- Test: `tests/test_gnf.py`, `tests/test_regularise.py`

**Interfaces:**
- Produces: `RateModel` = `dst_fit.models.Model` subclass fields + `control="rate"`, `guess(x, eta) -> np.ndarray`, `limits() -> dict[str, str]`; instances `POWER_LAW(k, n)`, `CARREAU_YASUDA(eta_0, eta_inf, lam, a, n)`, `CROSS(eta_0, eta_inf, lam, m)`, `HERSCHEL_BULKLEY(tau_y, k, n)` (rate space: `eta = tau_y/gd + k gd^(n-1)`); `regularise(model, scheme, **knobs) -> RateModel` with `scheme in {"papanastasiou", "biviscosity", "clamp"}`, knobs `m` (Papanastasiou, 1/s), `eta_max` (biviscosity/clamp, Pa s); `induced_error(model, reg, gammadot) -> np.ndarray` = `(eta_reg - eta)/eta`.

Gate (from roadmap 0.4): where regularisation is inactive (`gd*m >> 1`, `eta << eta_max`) the three agree with unregularised HB to 1e-6 relative; where active the disagreement equals the closed form: Papanastasiou `-tau_y exp(-m gd)/gd / eta`; clamp `(eta_max - eta)/eta` for `eta > eta_max`. SymPy limits in `derivations/gnf_limits.md`: Carreau-Yasuda → `eta_0` as `gd→0`, → `eta_inf` as `gd→∞`, → power law with `k = eta_0 lam^(n-1)` for `lam gd >> 1` and `eta_inf=0`; HB → power law at `tau_y=0`, → Bingham at `n=1`. Three-point numeric spot checks per model with hand-computed values shown in the test.

---

### Task 4: Identification layer with verdicts (unit 0.5a)

**Files:**
- Create: `src/biorheo/identify/fit.py`, `synthetic.py`
- Test: `tests/test_identify.py`

**Interfaces:**
- Produces: `fit_rate_curve(model: RateModel, gammadot, eta, *, level=0.95, seed=0) -> IdentifiedModel`; `IdentifiedModel(model, params: dict[str, ParamEstimate], fit: dst_fit.FitResult, envelope: Envelope)`; `ParamEstimate(value, lo, hi, verdict)`, `verdict in {"measurable", "bounded_only", "degenerate", "not_identifiable"}`; `synthetic_rate_curve(model, truth, gammadot, noise_cv, seed) -> eta`; `recovery_study(model, truth, gammadot, n_seeds=20) -> RecoveryReport(within_interval: dict[str,int])`.

Verdict rules (written before the run): `measurable` if profile interval finite on both sides and neither end within `1e-6` of a bound; `bounded_only` if one end is at a bound; `not_identifiable` if `profile_likelihood` reports the interval open on both sides; `degenerate` if `|corr| > 0.99` with another parameter AND interval open on at least one side. Gate: on 20 seeds at 3 % lognormal noise, truth inside the 95 % interval for ≥ 17/20 per parameter for `HERSCHEL_BULKLEY` on `gd in logspace(-2, 2, 25)`; on `gd in logspace(1, 2, 8)` (no yield plateau) `tau_y` returns `bounded_only` or `not_identifiable` in ≥ 18/20 seeds.

---

### Task 5: Structural family (unit 0.5b)

**Files:** `src/biorheo/constitutive/structural.py`, `tests/test_structural.py`
- Port Moore–Coussot lambda kinetics from `Desktop/Projects/bio_fluids_bridge` (read its source; do not re-derive from memory), as `dlambda_dt = (1 - lam)/t_r - alpha*gd*lam` with `eta = eta_inf (1 + beta lam)`; steady state `lam_ss = 1/(1 + alpha t_r gd)`, which is a `RateModel` with `param_names=("eta_inf","beta","alpha_tr")` — note `alpha` and `t_r` enter the steady curve only as the product, so a flow curve identifies the product only (the `parameters-that-enter-as-a-ratio` lesson made a gate: `recovery_study` on the 3-parameter steady model returns `degenerate`/`bounded_only` for `alpha` and `t_r` fitted separately, `measurable` for the product).

---

### Task 6: Injectability case study (unit 0.6)

**Files:** `src/biorheo/apps/injectability.py`, `tests/test_injectability.py`, `derivations/buckingham_reiner.md`, `case_studies/injectability_surrogates/` (script + md)

**Interfaces:**
- Produces: `flow_rate_hb(dp, R, L, tau_y, k, n) -> Q` (Herschel–Bulkley pipe flow, closed form for n=1 Buckingham–Reiner and numerical quadrature of the general HB velocity profile otherwise; plug radius `r_p = 2 tau_y L / dp`, returns 0 if `r_p >= R`); `pressure_for_flow_rate(Q, ...)` by `scipy.optimize.brentq`; `injection_force(Q, needle: Needle, syringe: Syringe, fluid: dict) -> Force` = `dp * pi R_s^2` plus friction term supplied by the caller; `force_with_uq(identified: IdentifiedModel, ...) -> (force, lo, hi)` by pushing each parameter through its profile interval and by sampling the fit covariance (N=2000, seed 0), reporting both.

Anchors (hand-computed, shown in test): Newtonian `Q = pi R^4 dp / (8 mu L)` at three parameter sets; power law `Q = pi n/(3n+1) R^3 (dp R/(2 k L))^(1/n)`; Bingham Buckingham–Reiner `Q = pi R^4 dp/(8 mu_p L) [1 - 4/3 (r_p/R) + 1/3 (r_p/R)^4]`. The report states which parameters the force is insensitive to (sensitivity `d ln F / d ln p` at the fit).

---

## Self-review

- Spec coverage: 0.1 → Task 1; 0.2 → Task 2; 0.3 (records) → written during Task 2 as the loader's fixtures plus a `records/` commit after Task 4 supplies the identifiability block for Carbopol; 0.4 → Task 3; 0.5 → Tasks 4–5; 0.6 → Task 6; 0.7 founder-only.
- Type consistency: `RateModel` is used by Tasks 3–6 with the same fields; `IdentifiedModel` by 4 and 6; `Envelope` by 2, 4, 6.
- Placeholder: `FIG_MD5` is filled in Step 2 of Task 1, not left.

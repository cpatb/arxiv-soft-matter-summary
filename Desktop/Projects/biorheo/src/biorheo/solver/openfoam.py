"""Export a registry record to an OpenFOAM ``transportProperties`` viscosity model and
read it back (the L3 solver bridge, first rung: Herschel-Bulkley and Casson).

Ground truth is the OpenFOAM v2606 (ESI) source, read in WSL at
``/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels/``
on 2026-09-05. The two ``calcNu`` bodies, verbatim:

HerschelBulkley/HerschelBulkley.C (lines 56-67)::

    min(nu0_, (tau0_ + k_*rtone*pow(tone*sr(), n_)) / max(sr(), VSMALL))

Casson/Casson.C (lines 60-70)::

    max(nuMin_, min(nuMax_, sqr(sqrt(tau0_/max(strainRate(), VSMALL)) + sqrt(m_))))

with ``VSMALL = 1.0e-300`` (OpenFOAM/primitives/Scalar/doubleScalar/doubleScalar.H:64),
``strainRate() = sqrt(2.0)*mag(symm(fvc::grad(U_)))`` (viscosityModel.C:62), which is
the engineering shear rate in simple shear, and coefficients read from the sub-dictionary
``<Model>Coeffs`` as ``k`` (dimViscosity), ``n`` (dimless), ``tau0``
(dimViscosity/dimTime), ``nu0`` (dimViscosity) for HB and ``m``, ``tau0``, ``nuMin``,
``nuMax`` for Casson (HerschelBulkley.C:88-91, Casson.C:90-93). Every one of them is
KINEMATIC, so a record's dynamic parameters are divided by a density the caller must
supply; records carry no density and this module will not guess one.

Parameter map (biorheo -> OpenFOAM), exact, no fitting involved::

    herschel_bulkley  eta = tau_y/gd + k gd^(n-1) = (tau_y + k gd^n)/gd
        tau0 = tau_y/rho   k = k/rho   n = n   nu0 = eta_max/rho
    casson            eta = (sqrt(tau_y/gd) + sqrt(eta_c))^2
        tau0 = tau_y/rho   m = eta_c/rho   nuMin = eta_c/rho   nuMax = eta_max/rho

``eta_max`` is required: OpenFOAM has no unregularised yield-stress law, and its scheme
is exactly biorheo's ``regularise(model, "clamp", eta_max=...)`` (``min(eta, eta_max)``),
so the exported law equals the clamped record, not the record. ``nuMin`` defaults to
the Casson high-rate limit ``eta_c``, the largest floor that never changes the law
(``eta >= eta_c`` for every rate), so the only knob is ``eta_max`` for both models.
``export_report`` says at which of the requested rates the cap binds; the caller
chooses ``eta_max`` above ``eta_at(gammadot_min)`` if the cap must stay outside the
measured envelope. ``coefficients`` refuses an ``eta_max`` at or below the viscosity at
the TOP of the envelope, where the exported law would be Newtonian everywhere.

``k`` is written into a slot OpenFOAM declares as ``dimViscosity``, and for n != 1 that
is dimensionally false: ``k_si`` is Pa s^n, so ``k/rho`` is m^2 s^(n-2). It is
numerically right because OpenFOAM carries ``tone = 1 s`` and ``rtone = 1 /s`` around
``pow(tone*sr(), n_)`` and everything here is SI. Do NOT "improve" the file by writing
an explicit ``[0 2 -1 0 0 0 0]`` bracket on ``k``: that would assert something untrue.

**THE ENVELOPE GUARD DOES NOT SURVIVE THE EXPORT.** ``FluidRecord.eta_at`` raises
``EnvelopeError`` outside the measured range; a dictionary cannot. OpenFOAM evaluates
``calcNu`` at whatever shear rate each cell has, and a run started from rest has every
cell at zero, decades below ``gammadot_min``. The only bound on the extrapolation is the
clamp, and the clamp is deliberately set outside the envelope. The exported file says so
in a comment, ``export_report`` returns ``envelope`` and ``residual_structure``, and
that is the whole of what this layer can enforce.

Gate (tests/test_openfoam_export.py): the round trip record -> text -> parsed
coefficients -> ``nu_openfoam`` reproduces ``rho * record.eta_at`` on 20 rates inside
the envelope to a relative 1e-10 wherever the cap does not bind, and the export
refuses a record with a ``degenerate`` or ``not_identifiable`` parameter unless
``allow_unidentified=True``. ``nu_openfoam`` is this module's transcription of the C++
above, so the round trip certifies the dictionary and the parameter map, not OpenFOAM.
Two anchors outside it: the test that runs ``foamDictionary -precision 17`` on the
written file (skipped without WSL) shows OpenFOAM's own parser recovers every
coefficient bit-exactly, and ``scripts/openfoam_couette_anchor.py`` runs ``simpleFoam``
on a 1-D plane-Couette case built from the exported dictionary and compares OpenFOAM's
own ``nu`` field and ``wallShearStress`` against ``record.eta_at`` without touching
``nu_openfoam`` at all (``data/openfoam/couette_anchor.md``).

That anchor also found the thing this module cannot fix. The export is exact; the SOLVE
is not. On a 20-cell 1-D problem ``simpleFoam``'s outer iteration converges at a rate set
by the local Bingham number ``B = tau_y/(k gd^n)``, and at the bottom of the Carbopol
envelope (B ~ 29) it returns a viscosity over 100 % wrong with ``rc = 0`` and no
warning. ``cap_binds`` does not catch it -- the clamp was ten times clear of the envelope
throughout. A record exported here is a correct dictionary, not a converged simulation.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..registry.schema import FluidRecord

OPENFOAM_VERSION = "v2606"
OPENFOAM_SRC = "/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels"
VSMALL = 1.0e-300  # doubleScalar.H:64

# biorheo model name -> (OpenFOAM transportModel name, coefficient keys in file order)
_MODELS = {
    "herschel_bulkley": ("HerschelBulkley", ("k", "n", "tau0", "nu0")),
    "casson": ("Casson", ("m", "tau0", "nuMin", "nuMax")),
}
_BAD_VERDICTS = ("degenerate", "not_identifiable")

_HEADER = """FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


class ExportRefused(ValueError):
    """The record's own identifiability sidecar says a parameter is not a measurement."""


def coefficients(record: FluidRecord, rho: float, eta_max: float) -> dict[str, float]:
    """Kinematic OpenFOAM coefficients for ``record`` at density ``rho`` (kg/m^3), with
    the clamp ``eta_max`` (Pa s). Pure arithmetic; no verdict check."""
    if not (rho > 0 and np.isfinite(rho)):
        raise ValueError(f"rho must be a positive density in kg/m^3, got {rho!r}")
    if not (eta_max > 0 and np.isfinite(eta_max)):
        raise ValueError(f"eta_max must be a positive viscosity in Pa s, got {eta_max!r}")
    if record.model.name not in _MODELS:
        raise NotImplementedError(
            f"no OpenFOAM {OPENFOAM_VERSION} export for model '{record.model.name}'; "
            f"available: {sorted(_MODELS)}")
    # A clamp at or below the viscosity at the TOP of the envelope makes the exported law
    # a constant-viscosity Newtonian fluid over the whole measured range, silently. Both
    # families are shear-thinning, so eta_at(gammadot_max) is that minimum. This also
    # SUBSUMES the Casson-specific check that used to live here (eta_max <= eta_c would
    # put nuMax at or below the nuMin floor and defeat OpenFOAM's max(nuMin, min(nuMax,.))
    # ordering): eta_c is Casson's infinite-rate limit and eta > eta_c at every finite
    # rate, so eta_c < eta_at(gammadot_max) always and this branch fires first.
    eta_hi = float(np.min(record.eta_at(np.array([record.envelope.gammadot_max]))))
    if eta_max <= eta_hi:
        raise ValueError(
            f"eta_max = {eta_max:g} Pa s is at or below the record's viscosity at the top "
            f"of its envelope ({eta_hi:g} Pa s at {record.envelope.gammadot_max:g} 1/s); "
            "the clamp would bind everywhere and the exported law would be Newtonian. "
            "For a Casson record this also covers the nuMin floor, which is lower still.")
    p = {name: float(v) for name, v in zip(record.model.param_names, record.params_si)}
    if record.model.name == "herschel_bulkley":
        return {"k": p["k"] / rho, "n": p["n"], "tau0": p["tau_y"] / rho, "nu0": eta_max / rho}
    assert eta_max > p["eta_c"], "unreachable: eta_c < eta_at(gammadot_max) <= eta_max"
    return {"m": p["eta_c"] / rho, "tau0": p["tau_y"] / rho,
            "nuMin": p["eta_c"] / rho, "nuMax": eta_max / rho}


def structured_flag(record: FluidRecord) -> dict | None:
    """The sidecar's ``residual_structure`` block, or None. A record can have every
    parameter graded ``measurable`` and still be a systematically wrong FORM -- 13 of the
    15 Carbopol records are flagged ``structured``, and Casson misfits 19 of 20 blood
    donors. Identifiability grades the parameters; this grades the model."""
    ident = record.identifiability
    rs = None if ident is None else ident.get("residual_structure")
    return rs if isinstance(rs, dict) else None


def check_exportable(record: FluidRecord, allow_unidentified: bool = False,
                     strict_form: bool = False) -> None:
    """Raise ``ExportRefused`` for a parameter the sidecar calls degenerate or not
    identifiable, or for a record without a sidecar, unless overridden. With
    ``strict_form`` also refuse a record whose residuals are flagged ``structured``,
    i.e. whose constitutive FORM is a systematic misfit even where every parameter is
    determined. That is off by default because it would refuse most of the registry --
    but it is written into the exported file either way, so a case that quotes an
    exported number can be traced back to a refuted form."""
    verdicts = record.verdicts
    if verdicts is None:
        if allow_unidentified:
            return
        raise ExportRefused(f"{record.id}: no identifiability sidecar; the parameters have "
                            "no verdicts. Pass allow_unidentified=True to export anyway.")
    bad = {k: v for k, v in verdicts.items() if v in _BAD_VERDICTS}
    if bad and not allow_unidentified:
        raise ExportRefused(f"{record.id}: parameter verdict(s) {bad} are not measurements; "
                            "an export would put an unidentified number into a solver. "
                            "Pass allow_unidentified=True to override, and say so in the case.")
    rs = structured_flag(record)
    if strict_form and rs and rs.get("structured"):
        raise ExportRefused(
            f"{record.id}: residual_structure is `structured` (z = {rs.get('z')}, "
            f"{rs.get('n_runs')} runs against {rs.get('expected_runs')} expected), so the "
            f"{record.model.name} FORM misfits this data systematically even though every "
            "parameter is measurable. Pass strict_form=False to export it anyway.")


def export_transport_properties(record: FluidRecord, rho: float, eta_max: float, *,
                                allow_unidentified: bool = False,
                                strict_form: bool = False) -> str:
    """The full ``constant/transportProperties`` text for ``record``."""
    check_exportable(record, allow_unidentified, strict_form)
    of_name, keys = _MODELS[record.model.name]
    c = coefficients(record, rho, eta_max)
    verd = record.verdicts or {}
    rs = structured_flag(record)
    lines = [_HEADER,
             f"// biorheo export of record {record.id} ({record.fluid})",
             f"// model {record.model.name}; kinematic coefficients = dynamic / rho, rho = {rho:.10g} kg/m^3",
             f"// clamp eta_max = {eta_max:.10g} Pa s (OpenFOAM {OPENFOAM_VERSION} form, {OPENFOAM_SRC}/{of_name})",
             f"// envelope {record.envelope.gammadot_min:g} to {record.envelope.gammadot_max:g} 1/s; "
             f"verdicts {verd if verd else 'none (allow_unidentified)'}",
             "// WARNING: the envelope above is a COMMENT. FluidRecord.eta_at refuses to",
             "// evaluate outside it; OpenFOAM will not. Every cell whose shear rate falls",
             "// below the lower bound -- which includes every cell of a run started from",
             "// rest -- is an extrapolation of this fit, bounded only by the clamp.",
             ("// WARNING: residual_structure is `structured` "
              f"(z = {rs.get('z'):.4g}, {rs.get('n_runs')} runs vs "
              f"{rs.get('expected_runs'):.4g} expected): the {record.model.name} FORM "
              "misfits this data systematically, even where every parameter is measurable."
              ) if rs and rs.get("structured") else
             "// residual_structure: not flagged (or no sidecar).",
             "",
             f"transportModel  {of_name};",
             "",
             f"{of_name}Coeffs",
             "{"]
    lines += [f"    {k:<8}{c[k]:.17g};" for k in keys]
    lines += ["}", "",
              "// ************************************************************************* //", ""]
    return "\n".join(lines)


def write_transport_properties(record: FluidRecord, rho: float, eta_max: float, path: str | Path, *,
                               allow_unidentified: bool = False,
                               strict_form: bool = False) -> Path:
    path = Path(path)
    path.write_text(export_transport_properties(record, rho, eta_max,
                                                allow_unidentified=allow_unidentified,
                                                strict_form=strict_form),
                    encoding="ascii", newline="\n")
    return path


# A dictionary entry: `key value;`, `key [dims] value;`, or OpenFOAM's full
# dimensionedScalar form `key key [dims] value;` (what its own tutorials write). The
# optional repeated name is why the second group exists; without it the tutorial form
# parsed to an EMPTY dict and `parse_transport_properties` returned success.
# An entry starts at the beginning of the text or just after `{`, `;` or a newline --
# NOT at the start of a LINE, which is what it used to require: a one-line block
# `{ k 0.0046; n 0.5; }` (legal OpenFOAM) then yielded only its first entry. The
# delimiter is a lookbehind so it is not consumed and consecutive entries still match.
_ENTRY = re.compile(r"(?:(?<=[{;\n])|^)\s*(\w+)\s+(?:\1\s+)?(?:\[[^\]]*\]\s*)?"
                    r"([-+0-9.eE]+)\s*;")
# `(?<![A-Za-z0-9_])` matters: without it `oldHerschelBulkleyCoeffs` (a stale block a
# user might keep beside the live one) matched first and its coefficients were returned
# with no error.
_MODEL_ENTRY = re.compile(r"(?<![A-Za-z0-9_])transportModel\s+(\w+)\s*;")


def parse_transport_properties(text: str, *, expect: str | None = None
                               ) -> tuple[str, dict[str, float]]:
    """``(transportModel, coefficients)`` from dictionary text. Accepts plain,
    dimensioned (``[0 2 -1 0 0 0 0] value``) and full dimensionedScalar
    (``k k [..] value``) entries; strips ``//`` comments first.

    ``expect`` names a biorheo model whose full coefficient set must be present; the
    parse then raises rather than returning a partial dict. Leave it None only when
    reading a dictionary whose contents are not known in advance.
    """
    text = re.sub(r"//[^\n]*", "", text)
    m = _MODEL_ENTRY.search(text)
    if m is None:
        raise ValueError("no transportModel entry")
    of_name = m.group(1)
    block = re.search(r"(?<![A-Za-z0-9_])" + re.escape(of_name) + r"Coeffs\s*\{([^}]*)\}", text)
    if block is None:
        raise ValueError(f"no {of_name}Coeffs block")
    coeffs = {k: float(v) for k, v in _ENTRY.findall(block.group(1))}
    if expect is not None:
        want = set(_MODELS[expect][1])
        if not want.issubset(coeffs):
            raise ValueError(f"{of_name}Coeffs is missing {sorted(want - set(coeffs))}; "
                             f"parsed {sorted(coeffs)}")
    return of_name, coeffs


def nu_openfoam(of_name: str, c: dict[str, float], sr) -> np.ndarray:
    """Kinematic viscosity OpenFOAM v2606 computes from ``c`` at strain rate ``sr``:
    the two ``calcNu`` bodies in the module docstring, transcribed term for term."""
    sr = np.asarray(sr, dtype=float)
    if of_name == "HerschelBulkley":
        return np.minimum(c["nu0"], (c["tau0"] + c["k"] * sr ** c["n"]) / np.maximum(sr, VSMALL))
    if of_name == "Casson":
        return np.maximum(c["nuMin"], np.minimum(
            c["nuMax"], (np.sqrt(c["tau0"] / np.maximum(sr, VSMALL)) + np.sqrt(c["m"])) ** 2))
    raise NotImplementedError(of_name)


def export_report(record: FluidRecord, rho: float, eta_max: float, gammadot, *,
                  allow_unidentified: bool = False, strict_form: bool = False) -> dict:
    """Round trip at the rates ``gammadot`` (all inside the envelope): export, parse,
    evaluate the OpenFOAM form, compare with ``rho * record.eta_at``. Returns the
    text, the parsed coefficients, the relative error where the cap does not bind, and
    the rates where it does (there the OpenFOAM law is ``eta_max`` by construction)."""
    text = export_transport_properties(record, rho, eta_max,
                                       allow_unidentified=allow_unidentified,
                                       strict_form=strict_form)
    of_name, c = parse_transport_properties(text, expect=record.model.name)
    gd = np.asarray(gammadot, dtype=float)
    nu_of = nu_openfoam(of_name, c, gd)
    nu_rec = record.eta_at(gd) / rho
    binds = nu_rec > eta_max / rho
    rel = np.abs(nu_of - nu_rec) / nu_rec
    return {
        "transportModel": of_name, "coefficients": c, "text": text, "gammadot": gd,
        "nu_openfoam": nu_of, "nu_record": nu_rec, "cap_binds": binds,
        "max_rel_err_unbound": float(rel[~binds].max()) if (~binds).any() else float("nan"),
        "n_bound": int(binds.sum()),
        "rate_where_cap_binds": float(gd[binds].max()) if binds.any() else None,
        # the two things the dictionary cannot enforce, carried so a caller can
        "envelope": (record.envelope.gammadot_min, record.envelope.gammadot_max),
        "residual_structure": structured_flag(record),
    }

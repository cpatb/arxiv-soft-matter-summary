"""Export a registry record to an OpenFOAM ``transportProperties`` viscosity model and
read it back (the L3 solver bridge: Herschel-Bulkley, Casson, Cross and Carreau-Yasuda).

Ground truth is the OpenFOAM v2606 (ESI) source, read in WSL at
``/usr/lib/openfoam/openfoam2606/src/transportModels/incompressible/viscosityModels/``
on 2026-09-05 (the two yield-stress families) and 2026-09-05 (the other two). The four
``calcNu`` bodies, verbatim:

HerschelBulkley/HerschelBulkley.C (lines 56-67)::

    min(nu0_, (tau0_ + k_*rtone*pow(tone*sr(), n_)) / max(sr(), VSMALL))

Casson/Casson.C (lines 60-70)::

    max(nuMin_, min(nuMax_, sqr(sqrt(tau0_/max(strainRate(), VSMALL)) + sqrt(m_))))

CrossPowerLaw/CrossPowerLaw.C (lines 53-57)::

    (nu0_ - nuInf_)/(scalar(1) + pow(m_*strainRate(), n_)) + nuInf_

BirdCarreau/BirdCarreau.C (lines 52-59)::

    nuInf_ + (nu0_ - nuInf_)*pow(scalar(1) + pow(k_*strainRate(), a_), (n_ - 1.0)/a_)

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
    cross             eta = eta_inf + (eta_0 - eta_inf)/(1 + (lam gd)^m)
        nu0 = eta_0/rho   nuInf = eta_inf/rho   m = lam   n = m        <- SWAPPED NAMES
    carreau_yasuda    eta = eta_inf + (eta_0 - eta_inf)[1 + (lam gd)^a]^((n-1)/a)
        nu0 = eta_0/rho   nuInf = eta_inf/rho   k = lam   n = n   a = a

**THE CROSS NAMES COLLIDE WITH THEIR MEANINGS SWAPPED.** biorheo's ``lam`` is a time
constant in seconds and OpenFOAM calls it ``m`` (declared ``dimTime``,
CrossPowerLaw.C:77); biorheo's ``m`` is a dimensionless exponent and OpenFOAM calls it
``n`` (``dimless``, line 78). A map written by matching names would write an exponent
into a time slot and a time into an exponent slot, and it would not raise anywhere --
both are plain scalars in the dictionary. ``tests/test_openfoam_export.py`` carries that
wrong map explicitly as a discriminator: it must FAIL the round trip by more than 1e-3,
or the round trip passing at 1e-10 has proved nothing.

**THE CLAMP IS NOT UNIVERSAL.** For Herschel-Bulkley and Casson, ``nu0``/``nuMax`` is a
regularisation the caller picks and ``eta_max`` is required. ``cross`` and
``carreau_yasuda`` have no yield stress: their ``eta_0`` is a FITTED PARAMETER and the
law is already bounded by it at every rate, so there is nothing to clamp, ``nu0`` is the
record itself, and ``coefficients`` REFUSES an ``eta_max`` rather than accepting and
ignoring one.

**``a`` IS OPTIONAL IN OPENFOAM, NOT ABSENT.** NEXT_STEPS item 6 was written believing
v2606's ``BirdCarreau`` has no ``a`` exponent and that the export would have to refuse
``a != 2``. It has one: ``a_`` is ``getOrDefault("a", dimensionedScalar("a", dimless, 2))``
(BirdCarreau.C:81-88, 120-124) and ``calcNu`` uses it, so the map is total and no refusal
is warranted. What IS true is the mirror image -- because ``a`` is optional with a silent
default of 2, a dictionary that omits it is a legal ``BirdCarreau`` running as a Carreau
fluid. So this module always WRITES ``a``, and ``parse_transport_properties`` fills a
missing one with 2 rather than reporting it missing, which is what OpenFOAM means by the
file. That ``a`` is really honoured was decided by running ``simpleFoam``, not by reading
the C++: see ``data/openfoam/couette_anchor.md``, gate G4.

``eta_max`` is required for the two yield-stress families: OpenFOAM has no unregularised yield-stress law, and its scheme
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
    "cross": ("CrossPowerLaw", ("nu0", "nuInf", "m", "n")),
    "carreau_yasuda": ("BirdCarreau", ("nu0", "nuInf", "k", "n", "a")),
}
# The families whose zero-shear plateau is a FITTED PARAMETER rather than a clamp, so
# they take no `eta_max` at all. See the module docstring, "The clamp is not universal".
_NO_CLAMP = ("cross", "carreau_yasuda")
# `a` is the one coefficient OpenFOAM reads with `getOrDefault` rather than `readEntry`
# (BirdCarreau.C:81-88, 120-124): a dictionary that omits it is a legal BirdCarreau that
# runs as a Carreau fluid. The parser mirrors that default; it does not raise.
_BIRDCARREAU_A_DEFAULT = 2.0
_OPTIONAL_COEFFS = {"BirdCarreau": {"a": _BIRDCARREAU_A_DEFAULT}}
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


def coefficients(record: FluidRecord, rho: float, eta_max: float | None = None
                 ) -> dict[str, float]:
    """Kinematic OpenFOAM coefficients for ``record`` at density ``rho`` (kg/m^3), with
    the clamp ``eta_max`` (Pa s) for the yield-stress families and ``None`` for the two
    that do not take one. Pure arithmetic; no verdict check."""
    if not (rho > 0 and np.isfinite(rho)):
        raise ValueError(f"rho must be a positive density in kg/m^3, got {rho!r}")
    if record.model.name not in _MODELS:
        raise NotImplementedError(
            f"no OpenFOAM {OPENFOAM_VERSION} export for model '{record.model.name}'; "
            f"available: {sorted(_MODELS)}")
    p = {name: float(v) for name, v in zip(record.model.param_names, record.params_si)}
    if record.model.name in _NO_CLAMP:
        # Refusing rather than ignoring is the point: eta_0 is a fitted parameter here,
        # so accepting a clamp and dropping it would put the caller's regularisation
        # request into the same silent-drop class this bridge exists to prevent.
        if eta_max is not None:
            raise ValueError(
                f"{record.id}: model '{record.model.name}' takes no eta_max. Its "
                f"zero-shear plateau eta_0 = {p['eta_0']:g} Pa s is a FITTED PARAMETER of "
                "the record, not a regularisation, and the law is bounded by it at every "
                "rate; there is nothing to clamp. Pass eta_max=None. (Herschel-Bulkley "
                "and Casson do require one: OpenFOAM has no unregularised yield-stress "
                "law.)")
        if record.model.name == "cross":
            # NAME COLLISION, both ways: biorheo's `lam` (a time) is OpenFOAM's `m`, and
            # biorheo's `m` (an exponent) is OpenFOAM's `n`. A map written by matching
            # names would put an exponent into a dimTime slot. CrossPowerLaw.C:56,75-78.
            return {"nu0": p["eta_0"] / rho, "nuInf": p["eta_inf"] / rho,
                    "m": p["lam"], "n": p["m"]}
        # carreau_yasuda -> BirdCarreau: a straight rename, `lam` -> `k` (dimTime), and
        # `a` is written ALWAYS even when it equals 2, because omitting it is legal and
        # silently means 2. BirdCarreau.C:55-59,77-88.
        return {"nu0": p["eta_0"] / rho, "nuInf": p["eta_inf"] / rho,
                "k": p["lam"], "n": p["n"], "a": p["a"]}
    if eta_max is None:
        raise ValueError(
            f"{record.id}: model '{record.model.name}' requires eta_max. OpenFOAM has no "
            "unregularised yield-stress law, so the exported dictionary is the CLAMPED "
            "record and the clamp has to be chosen here.")
    if not (eta_max > 0 and np.isfinite(eta_max)):
        raise ValueError(f"eta_max must be a positive viscosity in Pa s, got {eta_max!r}")
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
                     strict_form: bool = False,
                     strict_identifiability: bool = False) -> None:
    """Raise ``ExportRefused`` for a parameter the sidecar calls degenerate or not
    identifiable, or for a record without a sidecar, unless overridden. With
    ``strict_form`` also refuse a record whose residuals are flagged ``structured``,
    i.e. whose constitutive FORM is a systematic misfit even where every parameter is
    determined. That is off by default because it would refuse most of the registry --
    but it is written into the exported file either way, so a case that quotes an
    exported number can be traced back to a refuted form.

    ``strict_identifiability`` additionally refuses ``bounded_only``. It exists because a
    pre-registered gate FAILED on exactly that verdict (NEXT_STEPS item 9,
    ``data/horner2020/multistart_scatter.md``): the Herschel-Bulkley fit to donor G moves
    4.149 -- 415 % -- in ``tau_y`` under nothing but a change of the optimiser's seed,
    while ``tau_y`` there is graded ``bounded_only`` and so passes the default check.
    ``bounded_only`` means the profile ran into a bound on one side; it does NOT mean the
    value is pinned, and a flat direction can put it anywhere.

    It is OFF by default, and deliberately so: over the 43 export-accepted BLOOD records
    in the three families the registry ships (casson, cross, carreau_yasuda) the worst
    seed spread is 5.997e-05, four orders below the bar, and the one demonstrated
    counterexample is in blood Herschel-Bulkley, which is fitted for the model comparison
    and never shipped as a record.

    **Scope, now closed on both sets.** The first seed study
    (``scripts/multistart_seed_scatter.py``) fitted the 20 BLOOD curves only, leaving the
    15 Carbopol records unmeasured; ``scripts/multistart_seed_scatter_carbopol.py`` covers
    them (NEXT_STEPS item 10). All 15 are all-``measurable``, their worst seed spread is
    5.409e-06, their verdict sets do not change with the seed, and every one of the 225
    (seed, parameter, record) intervals still contains the value the record ships from
    ``carbopol/fit_hb.py`` -- so the loader's containment rule is not resting on seed 0
    being lucky.

    Turning this flag on refuses 7 of the 58 records the default accepts. All 7 are
    Carreau-Yasuda blood records; no Carbopol record carries a ``bounded_only`` at all.

    A caller about to QUOTE an exported parameter, rather than run a simulation with it,
    should pass ``strict_identifiability=True``."""
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
    if strict_identifiability and not allow_unidentified:
        weak = {k: v for k, v in verdicts.items() if v == "bounded_only"}
        if weak:
            raise ExportRefused(
                f"{record.id}: parameter verdict(s) {weak} are bounded on one side only. "
                "That is not the same as determined -- a `bounded_only` parameter has "
                "been measured moving 415 % under nothing but a change of the fitter's "
                "seed (see data/horner2020/multistart_scatter.md). Pass "
                "strict_identifiability=False to export it anyway, or "
                "allow_unidentified=True to override every verdict check.")
    rs = structured_flag(record)
    if strict_form and rs and rs.get("structured"):
        raise ExportRefused(
            f"{record.id}: residual_structure is `structured` (z = {rs.get('z')}, "
            f"{rs.get('n_runs')} runs against {rs.get('expected_runs')} expected), so the "
            f"{record.model.name} FORM misfits this data systematically even though every "
            "parameter is measurable. Pass strict_form=False to export it anyway.")


def export_transport_properties(record: FluidRecord, rho: float, eta_max: float | None = None, *,
                                allow_unidentified: bool = False,
                                strict_form: bool = False,
                                strict_identifiability: bool = False) -> str:
    """The full ``constant/transportProperties`` text for ``record``. ``eta_max`` is
    required for the yield-stress families and refused for the other two."""
    check_exportable(record, allow_unidentified, strict_form, strict_identifiability)
    of_name, keys = _MODELS[record.model.name]
    c = coefficients(record, rho, eta_max)
    verd = record.verdicts or {}
    rs = structured_flag(record)
    clamp_line = (
        f"// clamp eta_max = {eta_max:.10g} Pa s (OpenFOAM {OPENFOAM_VERSION} form, "
        f"{OPENFOAM_SRC}/{of_name})" if eta_max is not None else
        f"// no clamp: {record.model.name} has no yield stress and its zero-shear plateau "
        f"eta_0 is a fitted parameter, so nu0 below is the RECORD, not a regularisation "
        f"(OpenFOAM {OPENFOAM_VERSION} form, {OPENFOAM_SRC}/{of_name})")
    lines = [_HEADER,
             f"// biorheo export of record {record.id} ({record.fluid})",
             f"// model {record.model.name}; kinematic coefficients = dynamic / rho, rho = {rho:.10g} kg/m^3",
             clamp_line,
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


def write_transport_properties(record: FluidRecord, rho: float, eta_max: float | None,
                               path: str | Path, *,
                               allow_unidentified: bool = False,
                               strict_form: bool = False,
                               strict_identifiability: bool = False) -> Path:
    path = Path(path)
    path.write_text(export_transport_properties(record, rho, eta_max,
                                                allow_unidentified=allow_unidentified,
                                                strict_form=strict_form,
                                                strict_identifiability=strict_identifiability),
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
# MALFORMED and ABSENT are distinguished by a RESIDUE check: whatever `_ENTRY` does not
# consume from the coefficients block must be whitespace. It matters only because one
# coefficient is optional -- `_ENTRY` simply fails to match a broken entry, so before
# 2026-09-05 dropping the semicolon after `a` made it read as absent and `getOrDefault`'s
# 2 was filled in silently, giving a law 11 % wrong with no error anywhere. OpenFOAM's
# own dictionary reader raises on a present-but-unparsable entry; so does this now.
# A per-chunk `split(";")` test was tried first and does NOT work: the text after the
# last semicolon has no semicolon of its own, so exactly the dropped-semicolon case it
# was meant to catch looked well-formed.


def parse_transport_properties(text: str, *, expect: str | None = None,
                               defaults_applied: set[str] | None = None
                               ) -> tuple[str, dict[str, float]]:
    """``(transportModel, coefficients)`` from dictionary text. Accepts plain,
    dimensioned (``[0 2 -1 0 0 0 0] value``) and full dimensionedScalar
    (``k k [..] value``) entries; strips ``//`` comments first.

    ``expect`` names a biorheo model whose full coefficient set must be present; the
    parse then raises rather than returning a partial dict. Leave it None only when
    reading a dictionary whose contents are not known in advance.

    One coefficient is OPTIONAL in OpenFOAM itself: ``BirdCarreau``'s ``a`` is read with
    ``getOrDefault("a", 2)``, so a dictionary that omits it is legal and runs as a Carreau
    fluid. This function mirrors that -- a missing ``a`` is FILLED with 2, not reported
    missing -- because reporting it missing would disagree with the solver about what the
    file means. biorheo's own export always writes it. Pass ``defaults_applied`` a set to
    learn which keys were filled; a caller that quotes a number from a dictionary it did
    not write should say so.

    ABSENT and MALFORMED are different, and only the first gets a default. Every
    ``;``-terminated chunk of the coefficients block must be a well-formed entry, so a
    dropped semicolon raises instead of silently reading as an omitted key.
    """
    # Both OpenFOAM comment forms, block first so a `//` inside one cannot truncate it.
    # The block form was added 2026-09-06: the residue check below rejects anything
    # `_ENTRY` cannot consume, which turned a perfectly legal `/* ... */` inside a
    # coefficients block into a hard parse error. That only ever bit a dictionary
    # biorheo did not write -- which is exactly the case `expect=None` exists for.
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//[^\n]*", "", text)
    m = _MODEL_ENTRY.search(text)
    if m is None:
        raise ValueError("no transportModel entry")
    of_name = m.group(1)
    block = re.search(r"(?<![A-Za-z0-9_])" + re.escape(of_name) + r"Coeffs\s*\{([^}]*)\}", text)
    if block is None:
        raise ValueError(f"no {of_name}Coeffs block")
    body = block.group(1)
    coeffs = {k: float(v) for k, v in _ENTRY.findall(body)}
    residue = _ENTRY.sub("", body).strip()
    if residue:
        # Two different failures share this check and the caller needs to tell them
        # apart: a MALFORMED entry (a dropped semicolon, a comma decimal) means the file
        # is wrong, while an UNSUPPORTED CONSTRUCT (`$var` expansion, `#include`,
        # `#calc`) means the file is fine and this parser is not OpenFOAM. Both must
        # refuse -- silently skipping either would let an optional key fall back to its
        # default and change the law -- but the message says which.
        unsupported = sorted({t for t in ("$", "#") if t in residue})
        hint = (f" This looks like an OpenFOAM construct this parser does not implement "
                f"({', '.join(unsupported)}...): expansion, includes and computed "
                "entries are not evaluated here. Resolve them with "
                "`foamDictionary -expand` first." if unsupported else
                " A malformed entry is NOT treated as an absent one: OpenFOAM's "
                "dictionary reader raises here, and for an optional key like "
                "BirdCarreau's `a` a silent default would change the law with no error.")
        raise ValueError(
            f"{of_name}Coeffs contains text this parser cannot read as an entry: "
            f"{residue!r}.{hint}")
    for k, default in _OPTIONAL_COEFFS.get(of_name, {}).items():
        if k not in coeffs:
            coeffs[k] = default
            if defaults_applied is not None:
                defaults_applied.add(k)
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
    if of_name == "CrossPowerLaw":
        # (nu0_ - nuInf_)/(scalar(1) + pow(m_*strainRate(), n_)) + nuInf_
        # `m` is the TIME here and `n` the exponent -- the reverse of biorheo's names.
        return (c["nu0"] - c["nuInf"]) / (1.0 + (c["m"] * sr) ** c["n"]) + c["nuInf"]
    if of_name == "BirdCarreau":
        # nuInf_ + (nu0_ - nuInf_)*pow(1 + pow(k_*strainRate(), a_), (n_ - 1.0)/a_)
        # `a` absent means 2, exactly as OpenFOAM's getOrDefault does.
        a = c.get("a", _BIRDCARREAU_A_DEFAULT)
        return c["nuInf"] + (c["nu0"] - c["nuInf"]) * (1.0 + (c["k"] * sr) ** a) ** ((c["n"] - 1.0) / a)
    raise NotImplementedError(of_name)


def export_report(record: FluidRecord, rho: float, eta_max: float | None, gammadot, *,
                  allow_unidentified: bool = False, strict_form: bool = False,
                  strict_identifiability: bool = False) -> dict:
    """Round trip at the rates ``gammadot`` (all inside the envelope): export, parse,
    evaluate the OpenFOAM form, compare with ``rho * record.eta_at``. Returns the
    text, the parsed coefficients, the relative error where the cap does not bind, and
    the rates where it does (there the OpenFOAM law is ``eta_max`` by construction).
    With ``eta_max=None`` -- the only option for ``cross`` and ``carreau_yasuda`` -- no
    cap exists and ``cap_binds`` is False everywhere."""
    text = export_transport_properties(record, rho, eta_max,
                                       allow_unidentified=allow_unidentified,
                                       strict_form=strict_form,
                                       strict_identifiability=strict_identifiability)
    of_name, c = parse_transport_properties(text, expect=record.model.name)
    gd = np.asarray(gammadot, dtype=float)
    nu_of = nu_openfoam(of_name, c, gd)
    nu_rec = record.eta_at(gd) / rho
    binds = (nu_rec > eta_max / rho) if eta_max is not None else np.zeros(gd.shape, dtype=bool)
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

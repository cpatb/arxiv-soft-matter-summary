"""The deliverable: one page a client receives, generated from a registry record.

This module computes nothing new. Every number it prints comes from a function that is
anchored elsewhere in this package -- :mod:`biorheo.apps.injectability` for the force,
:mod:`biorheo.apps.design_space` for the bore exponent and the derating, the record's own
identifiability sidecar for the intervals and verdicts, and
:mod:`biorheo.solver.openfoam` for the export round trip. Its job is to put them on one
page in an order a buyer can act on, **and to refuse to print what the data does not
support**.

Why a module and not a script
-----------------------------
The venture's first product is a fixed-price study whose deliverable is a verdict, a
confidence interval and an identifiability statement (see
``../../biofluids_venture/docs/market_research_2026-09-06/00_SYNTHESIS.md``, rank 1). The
market research found no seller anywhere that states an uncertainty on a fitted
constitutive law, checks that a solver reproduces it, or tells the buyer which parameters
the answer does not depend on. That is the differentiator, so the page that carries it is
product code and is tested, not a one-off script.

What the report refuses
-----------------------
A refusal is the feature. In order of severity:

``no_sidecar``
    The record carries no identifiability block, so no interval and no verdict exist. The
    report prints the nominal force **labelled as a point estimate with no uncertainty**
    and states that plainly. It does not silently omit the interval.
``open_interval``
    A profile end sits on a box bound. That end is a BOUND, not a limit, so the affected
    force endpoint is printed as an inequality and the parameter is named. This repo's
    standing rule (``NEXT_STEPS.md``): do not read a sidecar ``hi`` as a limit without
    checking ``hi_at_bound``.
``non_commercial``
    The record's licence does not permit commercial use. The page is stamped and the
    reason is given. **Every record in the registry as of 2026-09-07 is in this state**,
    so a paid study needs the client's own curve.
``extrapolated``
    The device imposes a shear rate outside the fitted envelope. The distance is reported
    in decades and the force is labelled an extrapolation. It is not suppressed: an
    extrapolated number with its distance stated is more useful than no number, and the
    distance is the thing a reviewer asks for.

None of these is a warning appended at the end. Each changes how the headline is worded.

What is deliberately NOT here
-----------------------------
No new physics, no fitting, no plotting, and no aggregate over records -- a deliverable is
about one fluid. The omissions of the force model itself (barrel pressure drop, entrance
loss, kinetic term, wall slip, plunger stiction) are listed on the page verbatim from
:mod:`biorheo.apps.injectability`, not paraphrased, so they cannot drift apart.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

import numpy as np

from ..registry.schema import FluidRecord
from .design_space import bore_exponent, derating, newtonian_rule_error
from .injectability import (
    INSENSITIVE, Needle, Syringe, injection_force, sensitivities,
)

__all__ = ["DeviceSpec", "ReportData", "build", "render", "OMISSIONS"]

# Quoted, not paraphrased, from injectability.py's module docstring so the two cannot
# drift apart. `tests/test_report.py` asserts each phrase still appears there.
OMISSIONS = (
    "the barrel's own pressure drop",
    "entrance and exit losses",
    "the kinetic term",
    "the transient at plunger start",
    "wall slip",
)


@dataclass(frozen=True)
class DeviceSpec:
    """Everything about the device and the operating point. All SI.

    ``needle_radius`` is a radius, not a gauge: a gauge number is a range, not a
    dimension, and this package will not turn one into the other.
    """
    needle_radius: float
    needle_length: float
    barrel_radius: float
    flow_rate: float                     # m^3/s
    friction: float = 0.0                # N, a property of the syringe not the fluid
    density: float | None = None         # kg/m^3, only to report a wall Reynolds number
    force_budget: float | None = None    # N, for the derating section
    bore_ratio: float = 0.7              # the "finer needle" question
    bore_tolerance: float = 0.05         # declared manufacturing error on the bore
    label: str = "unnamed device"

    @property
    def needle(self) -> Needle:
        return Needle(radius=self.needle_radius, length=self.needle_length)

    @property
    def syringe(self) -> Syringe:
        return Syringe(barrel_radius=self.barrel_radius)


@dataclass
class ReportData:
    record_id: str
    fluid: str
    model: str
    spec: DeviceSpec
    flags: tuple[str, ...]
    nominal: dict
    parameters: dict
    hull: dict
    sensitivity: dict
    envelope: dict
    bore: dict
    derating: dict | None
    provenance: dict
    notes: tuple[str, ...] = field(default_factory=tuple)


def _intervals(record: FluidRecord) -> tuple[dict[str, tuple[float, float]], tuple[str, ...]]:
    block = record.identifiability
    if block is None:
        return {}, ()
    out, open_names = {}, []
    for name, p in block["parameters"].items():
        out[name] = (float(p["lo"]), float(p["hi"]))
        if p.get("lo_at_bound") or p.get("hi_at_bound"):
            open_names.append(name)
    return out, tuple(open_names)


def build(record: FluidRecord, spec: DeviceSpec) -> ReportData:
    """Assemble every number the page needs. Raises nothing a caller can act on: a
    record the report cannot fully describe produces FLAGS, not an exception, because a
    partial answer with its limits stated is the deliverable."""
    model, params = record.model, record.params_si
    needle, syringe = spec.needle, spec.syringe
    Q = spec.flow_rate

    flags: list[str] = []
    if record.identifiability is None:
        flags.append("no_sidecar")
    if not record.commercial_ok:
        flags.append("non_commercial")

    force = injection_force(Q, needle, syringe, model, params,
                            friction=spec.friction, rho=spec.density)

    iv, open_names = _intervals(record)
    if open_names:
        flags.append("open_interval")

    # one-at-a-time profile push-through hull. NOT a joint 95 % region: it is the box
    # that region is inscribed in, so it is wider. The page says which every time.
    # An interval END can lie just OUTSIDE the model's own box: the sidecar stores the
    # value the walk reported at the bound, and `carreau_yasuda`'s `a` comes back as
    # 10.001 against a bound of 10. Pushing that through the force solver raises. The end
    # IS the bound and the excess is how it was written down, so it is clamped -- and the
    # clamp is NOTED on the page rather than done silently.
    lo_b, hi_b = (np.asarray(b, dtype=float) for b in model.bounds)
    per_param, forces, clamped, unreachable = {}, [force.force], [], []
    for j, name in enumerate(model.param_names):
        if name not in iv:
            continue
        pair: list[float] = []
        for end, val in zip(("lo", "hi"), iv[name]):
            v = float(val)
            v_c = float(np.clip(v, lo_b[j], hi_b[j]))
            if v_c != v:
                clamped.append(f"`{name}` = {v:.6g} clamped to the model bound {v_c:.6g}")
            p = params.copy()
            p[j] = v_c
            try:
                pair.append(injection_force(Q, needle, syringe, model, p,
                                            friction=spec.friction).force)
            except ValueError as e:
                # The pressure inversion ran out of bracket: at this end of the interval
                # the fluid cannot be delivered through this needle at this rate at ANY
                # attainable pressure. That is a real answer and it is not a number, so
                # the hull end is reported OPEN rather than silently dropped.
                unreachable.append(f"`{name}` at its {end} end ({v_c:.6g}): {e}")
                pair.append(float("nan"))
        per_param[name] = tuple(pair)
        forces += [f for f in pair if np.isfinite(f)]
    if clamped:
        flags.append("interval_end_outside_model_bound")
    if unreachable:
        flags.append("hull_unbounded_above")
    hull = {"lo": float(min(forces)), "hi": float(max(forces)),
            "relative_width": float((max(forces) - min(forces)) / force.force),
            "per_parameter": per_param,
            "kind": "one-at-a-time profile push-through hull",
            "wider_than": "the joint 95 % region, which is an ellipse inscribed in this box",
            "clamped": tuple(clamped), "unreachable": tuple(unreachable)}

    sens = sensitivities(Q, needle, syringe, model, params, friction=spec.friction)
    sensitivity = {"d_lnF_d_lnp": sens,
                   "insensitive_to": tuple(k for k, v in sens.items() if abs(v) < INSENSITIVE),
                   "threshold": INSENSITIVE,
                   "dominant": max(sens, key=lambda k: abs(sens[k])) if sens else None}

    lo_env = float(record.envelope.gammadot_min)
    hi_env = float(record.envelope.gammadot_max)
    gd_w = float(force.wall_shear_rate)
    outside = not (lo_env <= gd_w <= hi_env)
    if outside:
        flags.append("extrapolated")
    decades = (float(np.log10(gd_w / hi_env)) if gd_w > hi_env
               else float(np.log10(lo_env / gd_w)) if gd_w < lo_env else 0.0)
    envelope = {"gammadot_min": lo_env, "gammadot_max": hi_env,
                "wall_shear_rate": gd_w, "outside": outside, "decades_past": decades}

    b = bore_exponent(Q, needle, syringe, model, params)
    nre = newtonian_rule_error(Q, needle, syringe, model, params, bore_ratio=spec.bore_ratio)
    bore = {"exponent": b, "newtonian_exponent": -4.0,
            "newtonian_rule_safe_within_10pct": bool(abs(nre["over_prediction"]) < 0.10),
            "tolerance_term": abs(b) * spec.bore_tolerance, **nre}

    der = None
    if spec.force_budget is not None:
        try:
            if iv:
                d = derating(spec.force_budget, needle, syringe, model, params, iv,
                             friction=spec.friction, open_parameters=open_names)
                der = {"budget_N": spec.force_budget,
                       "q_nominal_m3_s": d.q_nominal, "q_upper_m3_s": d.q_upper_bound,
                       "factor": d.factor, "binding_parameter": d.binding_parameter,
                       "binding_end": d.binding_end,
                       "open_parameters": list(d.open_parameters)}
            else:
                der = {"budget_N": spec.force_budget,
                       "error": "no interval: a derating factor needs one"}
        except ValueError as e:
            der = {"budget_N": spec.force_budget, "error": str(e)}

    parameters = {}
    for name in model.param_names:
        entry: dict[str, object] = {"value": float(record.parameters[name].si),
                                    "unit": record.parameters[name].unit}
        if record.identifiability is not None:
            p = record.identifiability["parameters"][name]
            entry.update({"lo": float(p["lo"]), "hi": float(p["hi"]),
                          "verdict": p["verdict"],
                          "lo_at_bound": bool(p.get("lo_at_bound")),
                          "hi_at_bound": bool(p.get("hi_at_bound"))})
        parameters[name] = entry

    provenance = {
        "sources": [{"doi": s.doi, "role": s.role} for s in record.sources],
        "licence": record.licence.name,
        "commercial_use": bool(record.commercial_ok),
        "generated_on": date.today().isoformat(),
        "generated_by": "biorheo.apps.report",
    }

    return ReportData(record_id=record.id, fluid=record.fluid, model=model.name,
                      spec=spec, flags=tuple(flags), nominal={
                          "force_N": force.force, "dp_Pa": force.dp,
                          "wall_stress_Pa": force.wall_stress,
                          "wall_shear_rate_1_s": force.wall_shear_rate,
                          "plug_fraction": force.plug_fraction,
                          "reynolds": force.reynolds, "friction_N": force.friction},
                      parameters=parameters, hull=hull, sensitivity=sensitivity,
                      envelope=envelope, bore=bore, derating=der, provenance=provenance)


def _headline(d: ReportData) -> str:
    """The one sentence, worded by the flags rather than followed by them."""
    F = d.nominal["force_N"]
    ml_min = d.spec.flow_rate * 6e7
    if "no_sidecar" in d.flags:
        core = (f"**{F:.3g} N** at {ml_min:.3g} mL/min — a POINT ESTIMATE with **no "
                "uncertainty**, because this record carries no identifiability block")
    elif "open_interval" in d.flags:
        names = ", ".join(sorted(
            n for n, p in d.parameters.items()
            if p.get("lo_at_bound") or p.get("hi_at_bound")))
        core = (f"**{F:.3g} N** at {ml_min:.3g} mL/min, with a one-at-a-time interval of "
                f"[{d.hull['lo']:.3g}, {d.hull['hi']:.3g}] N that is a **BOUND, not a limit** "
                f"— the profile for {names} ran into a box bound")
    else:
        core = (f"**{F:.3g} N** at {ml_min:.3g} mL/min, "
                f"[{d.hull['lo']:.3g}, {d.hull['hi']:.3g}] N over the 95 % "
                f"one-at-a-time parameter hull ({d.hull['relative_width'] * 100:.3g} % wide)")
    if "hull_unbounded_above" in d.flags:
        core += (". **The upper end of that interval is OPEN**: at one end of one "
                 "parameter's interval this fluid cannot be delivered through this needle "
                 "at this rate at any attainable pressure, so the interval quoted is a "
                 "lower bound on the uncertainty and not the uncertainty")
    if "extrapolated" in d.flags:
        core += (f". **The device runs {d.envelope['decades_past']:.2f} decades outside the "
                 "measured range**, so this is an extrapolation of the fitted law")
    return core


def render(d: ReportData) -> str:
    """The client-facing page. Every number in it is in :class:`ReportData`."""
    s = d.spec
    L = [f"# Injectability assessment — {d.fluid}", "",
         f"Device: **{s.label}**. Fluid record `{d.record_id}`, model `{d.model}`. "
         f"Generated {d.provenance['generated_on']} by `{d.provenance['generated_by']}`.", ""]

    if "non_commercial" in d.flags:
        L += ["> **NOT LICENSABLE AS A DELIVERABLE.** The source data for this record is "
              f"licensed *{d.provenance['licence']}* and does not permit commercial use. "
              "This page is a demonstration of the method. A paid study must run on the "
              "client's own curve or on a curve whose licence permits it.", ""]

    L += ["## The answer", "", _headline(d) + ".", "",
          "| quantity | value |", "|---|---|",
          f"| plunger force | {d.nominal['force_N']:.4g} N |",
          f"| pressure drop across the needle | {d.nominal['dp_Pa'] / 1e5:.4g} bar |",
          f"| wall shear stress | {d.nominal['wall_stress_Pa']:.4g} Pa |",
          f"| wall shear rate | {d.nominal['wall_shear_rate_1_s']:.4g} 1/s |"]
    if np.isfinite(d.nominal["reynolds"]):
        L.append(f"| wall Reynolds number | {d.nominal['reynolds']:.3g} |")
    if d.nominal["plug_fraction"] > 0:
        L.append(f"| unyielded plug, fraction of the bore | {d.nominal['plug_fraction']:.3g} |")
    if d.nominal["friction_N"]:
        L.append(f"| of which syringe friction (supplied, not modelled) | "
                 f"{d.nominal['friction_N']:.4g} N |")

    L += ["", "## Where your data stops and the device starts", "",
          f"The fitted envelope runs **{d.envelope['gammadot_min']:.4g} to "
          f"{d.envelope['gammadot_max']:.4g} 1/s**. This device imposes a wall shear rate of "
          f"**{d.envelope['wall_shear_rate']:.4g} 1/s**."]
    if d.envelope["outside"]:
        L += ["", f"That is **{d.envelope['decades_past']:.2f} decades outside the measured "
              "range.** Every number above is the fitted law evaluated where nothing was "
              "measured. The law was chosen to fit the data you have; nothing in that data "
              "constrains its shape here. Measuring one point at the device's own rate is "
              "the cheapest thing that would change this page."]
    else:
        L += ["", "The device operates **inside** the measured range. No extrapolation is "
              "involved in the force above."]

    L += ["", "## What the answer depends on — and what it does not", "",
          "`d ln F / d ln p` at the fit: a value near zero means the force does not care "
          "about that parameter at this operating point.", "",
          "| parameter | value | 95 % interval | verdict | sensitivity |",
          "|---|---|---|---|---|"]
    for name, p in d.parameters.items():
        iv = ("—" if "lo" not in p else
              f"[{p['lo']:.4g}, {p['hi']:.4g}]"
              + ("*" if p.get("lo_at_bound") or p.get("hi_at_bound") else ""))
        L.append(f"| `{name}` | {p['value']:.5g} {p['unit']} | {iv} | "
                 f"{p.get('verdict', '—')} | "
                 f"{d.sensitivity['d_lnF_d_lnp'].get(name, float('nan')):+.3f} |")
    if any(p.get("lo_at_bound") or p.get("hi_at_bound") for p in d.parameters.values()):
        L += ["", "\\* that interval ran into a box bound on at least one side. It is a "
              "**bound, not a confidence limit**."]
    if d.hull["unreachable"]:
        L += ["", "**At one or more interval ends this fluid cannot be delivered through "
              "this needle at this rate at any attainable pressure**, so the corresponding "
              "force endpoint does not exist and the interval above is a LOWER BOUND on the "
              "uncertainty: " + "; ".join(d.hull["unreachable"])]
    if d.hull["clamped"]:
        L += ["", "One or more interval ends were stored just outside the model's own "
              "parameter box and were clamped to it before the force was evaluated: "
              + "; ".join(d.hull["clamped"])
              + ". The end IS the bound; the excess is how it was written down."]
    insens = d.sensitivity["insensitive_to"]
    L += ["", (f"**The force is insensitive to {', '.join('`' + n + '`' for n in insens)}** "
               f"(|d ln F / d ln p| < {d.sensitivity['threshold']}). Effort spent pinning "
               "those parameters down will not move this answer."
               if insens else
               "**No parameter is insensitive here** — every one of them moves the force "
               f"by more than {d.sensitivity['threshold']} per decade of itself."),
          "",
          f"The dominant parameter is `{d.sensitivity['dominant']}` at "
          f"{d.sensitivity['d_lnF_d_lnp'][d.sensitivity['dominant']]:+.3f}."
          if d.sensitivity["dominant"] else ""]

    b = d.bore
    L += ["", "## A finer needle", "",
          f"`d ln F / d ln R` for this fluid at this operating point is **{b['exponent']:.4g}**, "
          f"against **{b['newtonian_exponent']:.0f}** for a Newtonian fluid. At a bore "
          f"{(1 - b['bore_ratio']) * 100:.0f} % narrower the force goes to "
          f"**{b['F_true']:.4g} N**; a Newtonian `R^-4` sizing rule would have predicted "
          f"**{b['F_newtonian_rule']:.4g} N**, "
          + (f"an over-prediction of **{b['over_prediction'] * 100:.3g} %**."
             if b["over_prediction"] > 0 else
             f"a difference of {b['over_prediction'] * 100:.3g} %."),
          "",
          ("The Newtonian rule is safe for this fluid here, to within 10 %."
           if b["newtonian_rule_safe_within_10pct"] else
           "**The Newtonian rule is not safe for this fluid here.** Sizing a needle with "
           "`R^-4` will mis-state the force by more than 10 %."),
          "",
          f"A **{s.bore_tolerance * 100:.0f} %** error in the bore moves the force by "
          f"**{b['tolerance_term'] * 100:.3g} %** "
          f"(|b| x the declared tolerance). Compare that with the "
          f"{d.hull['relative_width'] * 100:.3g} % parameter hull above: "
          + ("**the bore tolerance is the larger term**, and tightening the fit will not "
             "help until the geometry is pinned down."
             if b["tolerance_term"] > d.hull["relative_width"] else
             "**the parameter uncertainty is the larger term**, so a better fit — or "
             "more data — buys more than a tighter bore.")]

    if d.derating is not None:
        L += ["", "## Against a force budget", ""]
        if "error" in d.derating:
            L += [f"No derating factor: {d.derating['error']}"]
        else:
            dd = d.derating
            L += [f"Budget **{dd['budget_N']:.4g} N**. At the fitted parameters the budget is "
                  f"reached at **{dd['q_nominal_m3_s'] * 6e7:.4g} mL/min**; against the upper "
                  f"end of the parameter interval it is reached at "
                  f"**{dd['q_upper_m3_s'] * 6e7:.4g} mL/min**.", "",
                  f"**Derating factor {dd['factor']:.4g}** — a factor of "
                  f"{1 / dd['factor']:.3g} in deliverable flow rate, set by `"
                  f"{dd['binding_parameter']}` at its {dd['binding_end']} end. "
                  "The **ratio** is the usable number; the absolute rate inherits every "
                  "omission listed below."]
            if dd["open_parameters"]:
                L += ["", "This factor uses at least one interval that runs into a bound "
                      f"(`{'`, `'.join(dd['open_parameters'])}`), so it is itself a bound."]

    L += ["", "## What this calculation leaves out", "",
          "Laminar, fully developed, no-slip flow of a generalised-Newtonian fluid in a "
          "circular needle. Not included, and they do not all point the same way:", ""]
    L += [f"- {o}" for o in OMISSIONS]
    L += ["", "Syringe friction is a property of the device, not the fluid; it is added only "
          "if the caller supplies it. Wall slip in particular is a boundary condition, not a "
          "constitutive extrapolation, and no flag on this page covers it — detecting it "
          "needs two needle radii at the same wall stress.", "",
          "## Provenance", "",
          "| field | value |", "|---|---|",
          f"| record | `{d.record_id}` |",
          f"| constitutive model | `{d.model}` |",
          f"| licence | {d.provenance['licence']} |",
          f"| commercial use permitted | {'yes' if d.provenance['commercial_use'] else '**NO**'} |"]
    for src in d.provenance["sources"]:
        L.append(f"| source ({src['role']}) | https://doi.org/{src['doi']} |")
    L += [f"| needle | radius {s.needle_radius * 1e3:.4g} mm, length {s.needle_length * 1e3:.4g} mm |",
          f"| syringe barrel radius | {s.barrel_radius * 1e3:.4g} mm |",
          f"| flow rate | {s.flow_rate * 6e7:.4g} mL/min |",
          f"| flags | {', '.join(d.flags) if d.flags else 'none'} |", "",
          "A gauge number is a range, not a dimension; the radius above is the one supplied "
          "with the request and no standard was consulted.", ""]
    return "\n".join(x for x in L if x is not None)

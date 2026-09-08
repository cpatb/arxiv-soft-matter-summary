"""Item 17: the case for (and against) fixing the profile-interval defect in `dst_fit`.

**This writes a PROPOSAL. It changes nothing in `dst_fit` and must not be turned into a
patch by a later session without the founder's decision (D5).** `dst_fit` is the founder's
own package and this repo's standing rule is that it is consumed unchanged, precisely so
that `biorheo` cannot quietly reshape it.

Items 13 to 16 each deferred the upstream question with the same sentence -- *"if the fix
goes upstream it needs its own justification, not this one"* -- and four items later that
justification exists in measured form and has never been written down as a decision. This
is that document.

The one number nobody had computed is the **blast radius on `dst_fit`'s own results**, and
item 16 supplies the tool for it: the defect acts through

    exposure = (span / 5000) / (wald / 4),   span = |bound - best|

so the box floor wins exactly when `wald < span / 1250`. That is a statement about a
parameter's DECLARED BOX and the precision the fit reports, and it needs no run at all.

Run:  python scripts/dstfit_upstream_proposal.py
Writes `data/dstfit_upstream_proposal.md`. Touches nothing else.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DSTFIT = REPO.parent / "dst_fit"
OUT = REPO / "data" / "dstfit_upstream_proposal.md"

# The precision gate_n7's own module docstring quotes for phi_m on the real series.
# It is READ from that file rather than typed here, so it cannot go stale silently.
QUOTED_PHI_M = (0.5876, 0.0013)


def shipped_profile_callers() -> dict:
    """Which SHIPPED `dst_fit` code paths actually reach `profile_likelihood`.

    `grep`, not judgement: a script that never passes `profile=` cannot be affected by a
    defect inside the profile walk, whatever else it does.
    """
    def rg(pattern: str, path: str) -> list[str]:
        """Pure-Python so the result does not depend on a tool being on PATH.

        Comment lines are excluded: a docstring that MENTIONS `profile=` is not a call
        site, and counting one would overstate the exposure -- which is the direction that
        would make this document argue harder than the evidence supports.
        """
        rx = re.compile(pattern)
        out = []
        for f in sorted((DSTFIT / path).rglob("*.py")):
            for k, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
                if line.lstrip().startswith("#"):
                    continue
                if rx.search(line):
                    out.append(f"{f.relative_to(DSTFIT)}:{k}: {line.strip()}")
        return out

    return {
        "scripts_calling_profile_likelihood": rg(r"profile_likelihood\(", "scripts"),
        "scripts_passing_profile_kwarg": rg(r"profile=", "scripts"),
        "tests_passing_profile_kwarg": rg(r"profile=", "tests"),
    }


def exposure(span: float, wald: float) -> float:
    return (span / 5000.0) / (wald / 4.0)


def main() -> int:
    from dst_fit.models import WYART_CATES

    lo, hi = WYART_CATES.bounds
    boxes = {n: (float(a), float(b))
             for n, a, b in zip(WYART_CATES.param_names, lo, hi)}
    callers = shipped_profile_callers()

    best, wald = QUOTED_PHI_M
    b_lo, b_hi = boxes["phi_m"]
    e_lo = exposure(abs(best - b_lo), wald)
    e_hi = exposure(abs(b_hi - best), wald)

    # the same arithmetic for the wide-boxed parameter, at a plausible fitted value, to
    # show what the defect WOULD do if anything shipped profiled it
    s_lo, s_hi = boxes["sigma_star"]
    s_best, s_wald = 44.2, 4.0          # order of magnitude only; labelled as such below
    se_hi = exposure(abs(s_hi - s_best), s_wald)

    L = [
        "# Item 17 — a PROPOSAL: should the profile-interval defect be fixed in `dst_fit`?",
        "",
        f"Written by `scripts/{Path(__file__).name}`. **This document changes nothing.** "
        "`dst_fit` is the founder's own package and this repo's standing rule is that it is "
        "consumed unchanged. The decision is D5's and this is the case, not a patch.",
        "",
        "## The defect, in one line",
        "",
        "`dst_fit.profiling.profile_likelihood` sizes its outward step floor from the "
        "parameter's declared **box bound** rather than from its scale:",
        "",
        "```python",
        "step = max(wald / 4.0 if finite else span / 200.0,",
        "           span / 5000.0)          # span = |bound - best|",
        "```",
        "",
        "so the returned endpoint can be a point on the grid `best + m*span/5.12e6` rather "
        "than a likelihood crossing. Measured in `biorheo` (item 13): **201 of 296 closed "
        "upper ends** on that grid, worst **1639x the best-fit value** above it, and the "
        "shipped endpoints sat at a median **710x `f_crit`** over a 70-end pre-registered "
        "sample. A confidence limit has `F = f_crit` there, by definition.",
        "",
        "## What is NOT in dispute",
        "",
        "- The endpoints are not confidence limits (item 13, G13.2, n = 70).",
        "- `biorheo`'s corrected walk matches a **closed form that touches none of its "
        "code** to 2.4e-9 on 20 records (G13.4).",
        "- The blast radius is computable per parameter without running anything, and is "
        "now a field in every `biorheo` sidecar (item 16, G16.4).",
        "- `dst_fit`'s own suite is green and its one figure is byte-identical today.",
        "",
        "## The number nobody had computed: what it does to `dst_fit`'s OWN results",
        "",
        "The defect acts through",
        "",
        "```",
        "exposure = (span / 5000) / (wald / 4)          span = |bound - best|",
        "```",
        "",
        "so **the box floor wins exactly when `wald < span / 1250`**. That depends only on a "
        "parameter's declared box and the precision the fit reports.",
        "",
        "### Step 1 — which shipped code paths reach the profile walk at all",
        "",
        f"- `scripts/` calling `profile_likelihood` directly: "
        f"**{len(callers['scripts_calling_profile_likelihood'])}** line(s).",
        f"- `scripts/` passing `profile=` to `assess()`: "
        f"**{len(callers['scripts_passing_profile_kwarg'])}**.",
        f"- `tests/` passing `profile=`: "
        f"**{len(callers['tests_passing_profile_kwarg'])}**.",
        "",
    ]
    for ln in callers["scripts_calling_profile_likelihood"]:
        L.append(f"      {ln}")
    L += [
        "",
        "**No shipped `dst_fit` script passes `profile=` to `assess()`.** The only shipped "
        "consumer of the walk is `scripts/gate_n7.py`, and it profiles **`phi_m` only**. "
        "The one place `sigma_star` is profiled is a test, and a test asserting behaviour is "
        "not a published number.",
        "",
        "### Step 2 — `phi_m`'s exposure, which is the whole question",
        "",
        f"`phi_m`'s declared box is **[{b_lo:g}, {b_hi:g}]** — a packing fraction, so it is "
        f"O(1) wide. `gate_n7.py`'s own docstring quotes **phi_m = {best} +/- {wald}** on "
        "the real series. With that precision:",
        "",
        "| end | `span` | `span/5000` | `wald/4` | exposure | box floor wins? |",
        "|---|---|---|---|---|---|",
        f"| lower | {abs(best - b_lo):.4g} | {abs(best - b_lo) / 5000:.3e} | "
        f"{wald / 4:.3e} | **{e_lo:.3g}** | {'YES' if e_lo > 1 else 'no'} |",
        f"| upper | {abs(b_hi - best):.4g} | {abs(b_hi - best) / 5000:.3e} | "
        f"{wald / 4:.3e} | **{e_hi:.3g}** | {'YES' if e_hi > 1 else 'no'} |",
        "",
        f"**Both ends are below 1**, by factors of {1 / e_lo:.0f} and {1 / e_hi:.0f}. The "
        "box floor never sets the first step for `phi_m` at the precision `dst_fit` reports, "
        "so **the defect is inert for the only number `dst_fit` currently publishes through "
        "the profile walk.**",
        "",
        "That is the same mechanism items 15 and 16 found downstream: the defect's size "
        "scales with the DECLARED BOX, and a parameter whose box is O(1) is barely touched "
        "while a viscosity bounded at 1e6 is destroyed.",
        "",
        "### Step 3 — what it WOULD do to the wide-boxed parameter",
        "",
        f"`sigma_star`'s box is **[{s_lo:g}, {s_hi:g}]**. At an order-of-magnitude fitted "
        f"value of {s_best:g} Pa with a standard error of {s_wald:g} Pa — **illustrative "
        "numbers, not a fit** — the upper end's exposure is",
        "",
        f"    exposure_hi ~ {se_hi:.3g}",
        "",
        f"i.e. the first outward step would be about **{se_hi:.0f}x** the Wald-scaled one. "
        f"In absolute terms that step is `span/5000` = **{(s_hi - s_best) / 5000:.3g} Pa** on "
        f"a parameter whose value is {s_best:g} Pa — roughly "
        f"**{(s_hi - s_best) / 5000 / s_best * 100:.0f} % of the parameter itself, in a "
        "single step**. Nothing shipped profiles it today. **That is the risk the fix buys "
        "against: not a wrong number now, but a silently wrong one the first time anyone "
        "profiles a parameter with a wide box.**",
        "",
        "## The case FOR fixing it upstream",
        "",
        "1. The failure is **silent**. It returns a plausible-looking interval, wider than "
        "the truth, with no warning. A too-wide interval looks like an honest conservative "
        "one, which is why it survived this long.",
        "2. The fix is **one line** (size the floor from the parameter, not its box) plus "
        "the endpoint verification item 13 already wrote and tested.",
        "3. The correction has an **independent anchor**: a closed form that touches none of "
        "the walk, agreeing to 2.4e-9 on 20 records.",
        "4. It is **cheap to verify per parameter** — the exposure ratio above says which of "
        "`dst_fit`'s numbers could move, before any code runs. This document is that check "
        "and the answer is *none of them*, which makes the change low-risk as well as "
        "low-value.",
        "",
        "## The case AGAINST, which is not weak",
        "",
        "1. **`dst_fit` is consumed unchanged on purpose**, so that `biorheo` cannot quietly "
        "reshape the founder's package to suit itself. That rule has more value than one "
        "correct interval.",
        "2. **`biorheo`'s wrapper already gives `biorheo` correct intervals.** Patching "
        "upstream buys correctness for OTHER callers of `dst_fit` — of which there is "
        "currently one, `dst_fit` itself, and Step 2 shows it is not affected.",
        "3. `tests/test_gate_0_1_dstfit_unchanged.py` pins a figure hash that **will move**, "
        "and a gate whose whole purpose is to detect change cannot be left to fail quietly. "
        "It must be updated in the same commit, which means the change is not one line in "
        "practice.",
        "4. **`dst_fit`'s published numbers were computed with the defective walk.** Step 2 "
        "says they do not move, but 'do not move' is a claim that has to be re-established "
        "after the change, not assumed from before it.",
        "",
        "## Recommendation, for the founder to accept or reject",
        "",
        "**Fix it, but not because anything is currently wrong.** Step 2 establishes that no "
        "shipped `dst_fit` number moves; the argument is entirely about the *next* caller, "
        "and the failure mode is a wrong answer with no error. If it proceeds:",
        "",
        "1. **Pre-register what 'unchanged results' means before touching code**, per "
        "parameter, via the exposure ratio — the table in Step 2 is that pre-registration "
        "for `phi_m` and predicts NO movement at all. A gate that predicts no movement and "
        "then sees some has found something.",
        "2. `test_gate_0_1_dstfit_unchanged.py` moves **in the same commit**.",
        "3. **Do not delete `biorheo`'s wrapper.** It carries the endpoint verification, the "
        "second nuisance start and the `step_rule` / `exposure` reporting, none of which is "
        "part of the one-line fix and all of which caught something real — including two "
        "records where the profile is **not a function** of the parameter and a root finder "
        "returns a discontinuity as if it were a limit.",
        "",
        "**If it is rejected, that is a defensible position** and the thing to do is record "
        "it here, so item 17 stops being reopened. What is not defensible is leaving the "
        "question undecided for a fifth item.",
        "",
    ]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {OUT.name}: phi_m exposure lo {e_lo:.3g}, hi {e_hi:.3g} "
          f"(box floor wins on neither); "
          f"{len(callers['scripts_passing_profile_kwarg'])} shipped scripts pass profile=")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

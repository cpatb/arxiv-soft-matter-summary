"""test_shape_campaign.py -- regression tests for the shape CAMPAIGN.

`test_shape.py` tests the GEOMETRY (shape.py, 19 anchors). This file tests the runner,
the grader and the gates: the things that decide whether a number reaches a claims page.

Three rules this file is built around, all from lessons.md.

1. **A gate must be able to fail.** Every gate below is exercised twice -- once on data
   that must pass, once on data perturbed in the specific way the gate exists to catch.
   A gate that only ever sees passing data has not been tested (2026-09-07).
2. **A LOADER IS A LOADER WHEREVER IT LIVES** (2026-09-13). The helpers here build
   records through the SAME `grade_shape` functions the grader calls, never through a
   private copy that can drift from it.
3. **An ast-based scope guard on the import lines.** The one job of this session is the
   shape campaign; a test file that quietly grew an import of the closure or the ternary
   modules would be scope creep that nobody notices. The guard names what may be
   imported and fails on anything else.
"""

from __future__ import annotations

import ast
import json
import math
import os
import re
import runpy
import sys

import numpy as np
import pytest

import grade_shape
import shape

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------------------
# 0. THE SCOPE GUARD
# --------------------------------------------------------------------------------------

# What the shape campaign is allowed to depend on. `mixture_closure`, `kladder`,
# `ternary`, `ugap` and the additivity block are OUT OF SCOPE for this session by the
# author's instruction, and an import of one of them from the shape runner or grader
# would mean the campaign had quietly acquired a dependency on an unrelated open
# question.
ALLOWED_PROJECT_IMPORTS = {
    "packing", "psd", "measure", "shape", "phase2_phiJ", "phase2_shape",
    "grade_shape",
}
OUT_OF_SCOPE = {
    "mixture_closure", "mixture_model", "kladder", "phase12_kladder",
    "phase12_ternary", "phase12_mixture", "ugap", "predict_ternary",
    "predict_kladder", "predict_ugap", "transfer_kladder", "transfer_law",
    "shape_control", "firstpassage", "dst_window",
}


def _module_imports(path: str) -> set[str]:
    """Top-level module names imported by `path`, by AST rather than by regex.

    A regex over the source counts a module named in a docstring or a comment; the AST
    counts only imports the interpreter will execute.
    """
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                names.add(a.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:
                names.add(node.module.split(".")[0])
    return names


@pytest.mark.parametrize("fname", ["phase2_shape.py", "grade_shape.py",
                                   "test_shape_campaign.py"])
def test_scope_guard_no_out_of_scope_imports(fname):
    """The shape campaign imports nothing from the closure / ternary / u-gap blocks."""
    names = _module_imports(os.path.join(HERE, fname))
    leaked = names & OUT_OF_SCOPE
    assert not leaked, (
        f"{fname} imports out-of-scope modules {sorted(leaked)}. The shape campaign is "
        f"item 1; the closure (item 4) and the off-locus PSD (item 5) are explicitly "
        f"not this session's work.")


@pytest.mark.parametrize("fname", ["phase2_shape.py", "grade_shape.py"])
def test_scope_guard_project_imports_are_declared(fname):
    """Every first-party import is on the allowed list -- so a new dependency has to be
    added here deliberately rather than appearing."""
    names = _module_imports(os.path.join(HERE, fname))
    local = {n for n in names
             if os.path.exists(os.path.join(HERE, n + ".py"))}
    undeclared = local - ALLOWED_PROJECT_IMPORTS
    assert not undeclared, (
        f"{fname} imports first-party modules {sorted(undeclared)} that are not on "
        f"ALLOWED_PROJECT_IMPORTS. Add them here deliberately or remove them.")


def test_scope_guard_can_fail():
    """The guard itself is exercised on a file that DOES violate it, so a guard that
    silently matched nothing could not pass (lessons.md 2026-09-07: a check that must
    pass has decided nothing)."""
    src = "import kladder\nimport numpy as np\n"
    tree = ast.parse(src)
    names = {a.name.split(".")[0] for node in ast.walk(tree)
             if isinstance(node, ast.Import) for a in node.names}
    assert names & OUT_OF_SCOPE == {"kladder"}


# --------------------------------------------------------------------------------------
# 1. S0a -- the bit-identity gate, pinned
# --------------------------------------------------------------------------------------

def test_s0a_reference_values_are_the_committed_fabric_ones():
    """The S0a reference values are the F1 numbers in PREREG_FABRIC_2026-09-18.md, and
    they are asserted against the PREREG TEXT rather than against a second copy of the
    same constants -- so a silent edit to either side is caught."""
    with open(os.path.join(HERE, "PREREG_FABRIC_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    for (fam, cv), val in grade_shape.S0A_REFERENCE.items():
        assert repr(val) in text, (
            f"{fam} CV {cv}: grade_shape's reference {val!r} does not appear in "
            f"PREREG_FABRIC_2026-09-18.md. One of the two has drifted.")


def _gate_record(fam, cv, phi_por, phi_union=None, **kw):
    phi_union = phi_por if phi_union is None else phi_union
    r = dict(key=f"k_{fam}_{cv}", tag="shpgate", family=fam, cv=cv, alpha=1.0,
             phi=phi_union, phi_union=phi_union, phi_porosity=phi_por,
             z_all=4.0, z_nr=5.0, rattler_frac=0.2, n_contacts=2000,
             n_contacts_raw=2000, n_contacts_internal=0, dup_contact_frac=0.0,
             converged=True, seed=1, p_hat=1e-3, n_placed=1500, n_requested=1500,
             mass_over_vol_solid=1.0, eps_corrugation=0.0, span_over_deq=1.0)
    r.update(kw)
    return r


def _full_gate_set(perturb=None):
    recs = []
    for (fam, cv), val in grade_shape.S0A_REFERENCE.items():
        kw = dict((perturb or {}).get((fam, cv), {}))
        phi_por = kw.pop("phi_porosity", val)
        recs.append(_gate_record(fam, cv, phi_por, **kw))
    return recs


def test_s0a_passes_on_the_committed_values():
    g = grade_shape.g_s0a(_full_gate_set(), lambda *a, **k: None)
    assert g.passed is True


def test_s0a_fails_on_a_one_ulp_perturbation():
    """The gate is BIT identity, so one ulp must fail it. If this test passes with a
    tolerance somewhere, the gate is not the gate it says it is."""
    ref = grade_shape.S0A_REFERENCE[("uniform", 0.3)]
    nudged = np.nextafter(ref, 1.0)
    assert nudged != ref
    g = grade_shape.g_s0a(_full_gate_set({("uniform", 0.3):
                                            {"phi_porosity": nudged}}),
                          lambda *a, **k: None)
    assert g.passed is False


def test_s0a_fails_when_a_gate_cell_is_missing():
    """A gate that did not run cannot pass -- the failure mode where two of three cells
    are green and the third was never launched."""
    recs = [r for r in _full_gate_set() if r["family"] != "lognormal"]
    g = grade_shape.g_s0a(recs, lambda *a, **k: None)
    assert g.passed is False


def test_s0a_fails_if_the_cells_were_not_run_at_alpha_1():
    g = grade_shape.g_s0a(_full_gate_set({("uniform", 0.3): {"alpha": 1.5}}),
                          lambda *a, **k: None)
    assert g.passed is False


# --------------------------------------------------------------------------------------
# 2. S0b -- and the anti-vacuity assertion
# --------------------------------------------------------------------------------------

def test_s0b_passes_at_the_measured_separation():
    recs = _full_gate_set()
    for r in recs:
        r["phi_union"] = r["phi_porosity"] + 1e-14
    g = grade_shape.g_s0b(recs)
    assert g.passed is True


def test_s0b_fails_if_the_two_routes_are_EXACTLY_equal():
    """A deviation of exactly zero means `phi_union` and `phi_porosity` are the same
    calculation, so the gate would be vacuous rather than passing (lessons.md
    2026-09-06). It must fail, not pass."""
    g = grade_shape.g_s0b(_full_gate_set())     # phi_union == phi_porosity by default
    assert g.passed is False


def test_s0b_fails_when_the_routes_disagree_physically():
    recs = _full_gate_set()
    for r in recs:
        r["phi_union"] = r["phi_porosity"] * 2.0    # the un-corrected double count
    g = grade_shape.g_s0b(recs)
    assert g.passed is False


# --------------------------------------------------------------------------------------
# 3. S0c -- the clump construction, which S0a can no longer test
# --------------------------------------------------------------------------------------

def _prod_record(alpha, fam="lognormal", cv=0.3, seed=1, p_hat=1e-3, phi=0.62, **kw):
    ds = shape.sub_sphere_diameter(1.0, alpha)
    ratio = 2.0 * ds ** 3
    r = dict(key=f"p_{fam}_{cv}_{alpha}_{seed}_{p_hat}", tag="shp", family=fam, cv=cv,
             alpha=alpha, seed=seed, p_hat=p_hat, p_hat_target=p_hat, phi=phi,
             phi_union=phi, phi_porosity=phi * ratio, converged=True,
             n_placed=1500, n_requested=1500, mass_over_vol_solid=1.0,
             eps_corrugation=shape.corrugation(1.0, alpha),
             span_over_deq=shape.span_factor(alpha),
             cv_realised=cv, d32_realised=1.3, d43_realised=1.5,
             z_all=4.0, z_nr=5.0, rattler_frac=0.2, n_contacts=2000,
             n_contacts_raw=2200, n_contacts_internal=0, dup_contact_frac=0.09)
    r.update(kw)
    return r


def test_s0c_passes_on_a_correctly_built_clump_set():
    recs = [_prod_record(a) for a in (1.25, 1.5, 1.75, 2.0)]
    g = grade_shape.g_s0c(recs, lambda *a, **k: None)
    assert g.passed is True


def test_s0c_fails_on_a_wrong_sub_sphere_diameter():
    """The predicted ratio is alpha-DEPENDENT (2.0000 -> 1.0000), which is what makes
    this gate able to fail: no single wrong constant satisfies five different values."""
    recs = [_prod_record(a) for a in (1.25, 1.5, 1.75, 2.0)]
    recs[1]["phi_porosity"] = recs[1]["phi_union"] * 2.0 * \
        shape.sub_sphere_diameter(1.0, 1.75) ** 3      # the ratio for the WRONG alpha
    g = grade_shape.g_s0c(recs, lambda *a, **k: None)
    assert g.passed is False


def test_s0c_fails_on_a_wrong_density():
    recs = [_prod_record(a) for a in (1.25, 1.5, 1.75, 2.0)]
    recs[2]["mass_over_vol_solid"] = 2.0 * shape.sub_sphere_diameter(1.0, 1.75) ** 3
    g = grade_shape.g_s0c(recs, lambda *a, **k: None)
    assert g.passed is False


def test_s0c_fails_when_no_clump_cell_exists():
    """If alpha = 1 is a sphere and no alpha > 1 cell ran, the clump construction is
    UNTESTED, and that must not read as a pass."""
    g = grade_shape.g_s0c([_prod_record(1.0)], lambda *a, **k: None)
    assert g.passed is False


def _union_volume_by_quadrature(ds: float, s: float) -> float:
    """Union volume of two spheres of diameter `ds` whose centres are `s` apart,
    computed by QUADRATURE and touching none of shape.py.

    For two spheres on the z axis the cross-section at height z is a disk whose radius
    is the LARGER of the two spheres' local radii, so

        V = integral pi * max(r^2 - (z-c1)^2, r^2 - (z-c2)^2, 0) dz

    with the max clipped at zero outside both. This is independent of `g_of_u`, of the
    lens formula, and of every constant in the module -- which is the point: gate S0c
    calls `shape.sub_sphere_diameter` on BOTH sides of its comparison, so it cannot
    detect an error in the closed form (referee, 2026-09-18). This is the outside
    anchor that can.
    """
    from scipy.integrate import quad
    r = 0.5 * ds
    c1, c2 = -0.5 * s, +0.5 * s

    def area(z):
        a1 = r * r - (z - c1) ** 2
        a2 = r * r - (z - c2) ** 2
        return math.pi * max(a1, a2, 0.0)

    lo, hi = c1 - r, c2 + r
    # split at the waist and at each centre so the integrand's kinks are on panel edges
    pts = sorted({lo, c1, 0.0, c2, hi})
    total = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        v, _e = quad(area, a, b, limit=200, epsabs=1e-13, epsrel=1e-13)
        total += v
    return total


@pytest.mark.parametrize("alpha", [1.0, 1.1, 1.25, 1.5, 1.75, 1.9, 2.0])
def test_clump_volume_against_independent_quadrature(alpha):
    """THE OUTSIDE ANCHOR. `shape.clump_volume` must equal (pi/6)*d_eq^3 by design; that
    identity is checked elsewhere but runs entirely through `g_of_u`. Here the SAME
    construction is measured by a quadrature that never imports the closed form, so a
    wrong coefficient in `g_of_u` moves one side and not the other."""
    d_eq = 1.0
    ds = shape.sub_sphere_diameter(d_eq, alpha)
    s = (alpha - 1.0) * ds
    v_quad = _union_volume_by_quadrature(ds, s)
    v_target = math.pi / 6.0 * d_eq ** 3
    assert abs(v_quad - v_target) < 1e-10, (alpha, v_quad, v_target)


def test_the_quadrature_anchor_can_fail():
    """A discriminating negative: perturb g_of_u's coefficient and the anchor must
    reject the construction it produces. Without this the anchor could be passing
    because both sides are wrong in the same way."""
    d_eq = 1.0
    alpha = 1.5
    u = alpha - 1.0
    bad_g = math.pi / 3.0 - (math.pi / 12.0) * (2.0 + u) * (1.0 - u) ** 2 * 1.01
    ds_bad = d_eq * ((math.pi / 6.0) / bad_g) ** (1.0 / 3.0)
    v = _union_volume_by_quadrature(ds_bad, u * ds_bad)
    assert abs(v - math.pi / 6.0) > 1e-4, v


def test_s0c_docstring_does_not_claim_to_test_the_closed_form():
    """S0c calls `shape.sub_sphere_diameter` on both sides, so `phi_porosity/phi_union`
    equals `2*(ds/d_eq)^3` for ANY value that function returns. The gate tests the
    CONSTRUCTION and the density, not the geometry, and its text must say so."""
    src = open(os.path.join(HERE, "grade_shape.py"), encoding="utf-8").read()
    i = src.index("def g_s0c")
    body = src[i:src.index("def g_s1")]
    # Join adjacent string literals so a phrase split across two source lines is found
    # as the reader will see it. Without this the test greps the SOURCE and the prose
    # is the CONCATENATION, which is how a correct message can fail a text assertion.
    body = re.sub(r'"\s*\n\s*(f?)"', "", body)
    assert "wrong g(u)" not in body, (
        "S0c cannot detect a wrong g(u): it shares `sub_sphere_diameter` with the "
        "construction. The independent anchor is "
        "test_clump_volume_against_independent_quadrature.")
    assert "cannot detect an error in the closed form" in body


def test_s0c_predicted_ratios_match_the_prereg_table():
    """The five predicted double counts named in the prereg amendment, recomputed."""
    expect = [2.0000, 1.4629, 1.1852, 1.0449, 1.0000]
    got = [2.0 * shape.sub_sphere_diameter(1.0, a) ** 3 for a in shape.ALPHA_LADDER]
    assert np.allclose(got, expect, atol=5e-5), got


# --------------------------------------------------------------------------------------
# 4. S1 -- the PSD must not move
# --------------------------------------------------------------------------------------

def test_s1_passes_when_the_psd_is_invariant():
    recs = [_prod_record(a) for a in (1.0, 1.5, 2.0)]
    g = grade_shape.g_s1(recs)
    assert g.passed is True


def test_s1_fails_when_a_descriptor_moves():
    recs = [_prod_record(a) for a in (1.0, 1.5, 2.0)]
    recs[2]["d32_realised"] = 1.3 + 1e-9
    g = grade_shape.g_s1(recs)
    assert g.passed is False


def test_s1_fails_when_there_is_no_alpha_1_rung_to_compare_against():
    """Zero comparisons is not a pass: an invariance gate with nothing to compare is the
    exact-zero signature (lessons.md 2026-09-09)."""
    recs = [_prod_record(a) for a in (1.5, 2.0)]
    g = grade_shape.g_s1(recs)
    assert g.passed is False


# --------------------------------------------------------------------------------------
# 5. The structural tag guard and the duplicate-key refusal
# --------------------------------------------------------------------------------------

def test_dedup_refuses_duplicates_that_disagree():
    """lessons.md 2026-09-09: 15 duplicated keys in this project returned a different
    phi from the same cell, and which one survived was set by the filename sort."""
    a = _prod_record(1.5)
    b = dict(a)
    b["phi"] = a["phi"] + 1e-4
    with pytest.raises(SystemExit):
        grade_shape.dedup([a, b], lambda *x, **k: None)


def test_dedup_accepts_exact_duplicates():
    a = _prod_record(1.5)
    out = grade_shape.dedup([a, dict(a)], lambda *x, **k: None)
    assert len(out) == 1


def test_forbidden_prefixes_cover_every_sibling_campaign():
    """The tag guard names the campaigns whose shards must never be admitted. If a new
    campaign appears in out/ with a phase2_ prefix this test is where it gets noticed."""
    for p in ("phase2_mixr", "phase2_tern", "phase2_kl", "phase2_ug", "phase2_fab"):
        assert p in grade_shape.FORBIDDEN_PREFIXES


def test_tag_guard_rejects_a_foreign_shard(tmp_path, monkeypatch):
    bad = tmp_path / "phase2_fab_lognormal.jsonl"
    bad.write_text(json.dumps(_prod_record(1.5)) + "\n", encoding="utf-8")
    monkeypatch.setattr(grade_shape, "SHARD_GLOB", str(tmp_path / "phase2_*.jsonl"))
    with pytest.raises(SystemExit):
        grade_shape.load_records(lambda *a, **k: None)


def test_gate_and_probe_tags_are_held_out_of_production(tmp_path, monkeypatch):
    """The `shpgate` cells are alpha = 1 reproductions of records already on disk and
    the `shptime` cells ran at a different N. Pooling either into an S3 average would be
    averaging a different question."""
    f = tmp_path / "phase2_shp_mixed.jsonl"
    rows = [_gate_record("uniform", 0.3, 0.61),
            dict(_prod_record(1.5), tag="shptime", n_requested=5400),
            _prod_record(1.5)]
    f.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    monkeypatch.setattr(grade_shape, "SHARD_GLOB", str(tmp_path / "phase2_shp*.jsonl"))
    prod, gate, _files = grade_shape.load_records(lambda *a, **k: None)
    assert len(prod) == 1 and len(gate) == 2
    assert all(r["tag"] == "shp" for r in prod)


def test_records_without_an_alpha_field_are_refused(tmp_path, monkeypatch):
    """Every record written before 2026-09-18 lacks `alpha`, and those are exactly the
    sphere records this campaign's alpha = 1 rung REPRODUCES -- so admitting one would
    let the reproduction stand in for the measurement."""
    r = _prod_record(1.5)
    del r["alpha"]
    f = tmp_path / "phase2_shp_noalpha.jsonl"
    f.write_text(json.dumps(r) + "\n", encoding="utf-8")
    monkeypatch.setattr(grade_shape, "SHARD_GLOB", str(tmp_path / "phase2_shp*.jsonl"))
    with pytest.raises(SystemExit):
        grade_shape.load_records(lambda *a, **k: None)


# --------------------------------------------------------------------------------------
# 6. The runner's cell key
# --------------------------------------------------------------------------------------

def test_cell_key_carries_alpha_and_cannot_alias_a_sphere_key():
    import phase2_phiJ
    import phase2_shape
    sphere = phase2_phiJ.cell_key(0.3, 0.3, 1, 1e-3, 1500, "uniform")
    shp1 = phase2_shape.cell_key(0.3, 0.3, 1, 1e-3, 1500, "uniform", 1.0)
    shp2 = phase2_shape.cell_key(0.3, 0.3, 1, 1e-3, 1500, "uniform", 1.5)
    assert shp1 != sphere, (
        "the shape runner's alpha = 1 key must NOT collide with the sphere key it "
        "reproduces, or a resume would serve the sphere record back as the gate's own "
        "answer")
    assert shp1 != shp2 and "a1.000" in shp1 and "a1.500" in shp2


def test_runner_rejects_an_alpha_outside_the_two_sphere_family():
    with pytest.raises(ValueError):
        shape.sub_sphere_diameter(1.0, 2.5)


# --------------------------------------------------------------------------------------
# 7. The grain-level contact graph -- the rewrite S0a checks on real data
# --------------------------------------------------------------------------------------

class _FakeInteraction:
    def __init__(self, i, j):
        self.id1, self.id2, self.isReal = i, j, True
        self.phys = type("P", (), {"normalForce": type("V", (), {
            "norm": staticmethod(lambda: 1.0)})()})()


class _FakeBody:
    def __init__(self, bid, clump_id=None, is_clump=False):
        self.id = bid
        self.isClump = is_clump
        self.isClumpMember = clump_id is not None
        self.clumpId = clump_id if clump_id is not None else bid


class _FakeO:
    def __init__(self, bodies, inters):
        self.bodies, self.interactions = bodies, inters


def test_contact_graph_is_a_noop_for_spheres():
    """For a sphere pack the deduplicated grain graph must equal the old
    per-interaction one body for body -- asserted, not assumed."""
    import packing
    bodies = [_FakeBody(i) for i in range(4)]
    inters = [_FakeInteraction(0, 1), _FakeInteraction(1, 2), _FakeInteraction(2, 3)]
    stats = {}
    counts, neigh, forces = packing.contact_graph(_FakeO(bodies, inters), stats=stats)
    assert counts == {0: 1, 1: 2, 2: 2, 3: 1}
    assert neigh[1] == {0, 2}
    assert forces.size == 3
    assert stats["n_raw"] == 3 and stats["dup_frac"] == 0.0


def test_contact_graph_deduplicates_a_clump_pair():
    """Two clumps whose four member spheres all touch are ONE grain contact, and the
    duplicate fraction is 0.75 -- the number probe_clump_s0.sh measured at alpha = 1."""
    import packing
    # bodies 0,1 are members of clump 100; bodies 2,3 of clump 101
    bodies = [_FakeBody(0, 100), _FakeBody(1, 100), _FakeBody(2, 101),
              _FakeBody(3, 101), _FakeBody(100, is_clump=True),
              _FakeBody(101, is_clump=True)]
    inters = [_FakeInteraction(a, b) for a in (0, 1) for b in (2, 3)]
    grain_of = {0: 100, 1: 100, 2: 101, 3: 101}
    stats = {}
    counts, neigh, _f = packing.contact_graph(_FakeO(bodies, inters),
                                              grain_of=grain_of, stats=stats)
    assert counts == {100: 1, 101: 1}, counts
    assert neigh[100] == {101}
    assert stats["n_raw"] == 4 and stats["n_distinct"] == 1
    assert stats["dup_frac"] == 0.75


def test_contact_graph_drops_internal_member_contacts():
    import packing
    bodies = [_FakeBody(0, 100), _FakeBody(1, 100), _FakeBody(100, is_clump=True)]
    inters = [_FakeInteraction(0, 1)]
    stats = {}
    counts, _n, _f = packing.contact_graph(
        _FakeO(bodies, inters), grain_of={0: 100, 1: 100}, stats=stats)
    assert counts == {100: 0}
    assert stats["n_internal"] == 1 and stats["n_raw"] == 0


def test_z_all_denominator_is_grains_not_bodies():
    """The bug this rewrite exists to prevent: with clumps, `len(counts)` must be the
    number of GRAINS. Counting bodies would put the clump body and both members in the
    denominator and make z_all wrong by a factor of three."""
    import packing
    bodies = [_FakeBody(0, 100), _FakeBody(1, 100), _FakeBody(100, is_clump=True)]
    counts, _n, _f = packing.contact_graph(_FakeO(bodies, []),
                                           grain_of={0: 100, 1: 100})
    assert len(counts) == 1


# --------------------------------------------------------------------------------------
# 7b. S3 -- the paired interaction statistic
# --------------------------------------------------------------------------------------

def _tab(interaction, noise=0.0, seeds=(1, 2, 3), alphas=(1.0, 2.0),
         cvs=(0.0, 0.5), fams=("lognormal", "uniform", "bidisperse")):
    """A synthetic phi_J table with a KNOWN interaction, so the estimator can be
    checked against a truth the code under test did not produce."""
    rng = np.random.default_rng(0)
    tab = {}
    for fam in fams:
        for a in alphas:
            for c in cvs:
                d = {}
                for s in seeds:
                    base = 0.62 - 0.05 * c
                    # the interaction: an extra CV effect that appears only at alpha_hi
                    if a == alphas[-1] and c == cvs[-1]:
                        base += interaction
                    # (phi_J, se, b_free, n_points, slope) -- the SAME tuple shape
                    # `grade_shape.phi_j_table` emits. The slope carries the
                    # contact-stiffness prediction (`d_eq/ds` = 1.2599 from alpha 1 to
                    # alpha 2), so it is built to satisfy it here.
                    slope = -5.0 * (1.2599 if a == alphas[-1] else 1.0)
                    d[s] = (base + noise * rng.normal(), 1e-4, float("nan"), 3, slope)
                tab[(fam, c, a)] = d
    return tab


def test_s3_recovers_a_known_interaction():
    g, s3 = grade_shape.g_s3(_tab(0.010), lambda *a, **k: None)
    for fam, (m, se, tcrit, res, n) in s3["inter"].items():
        assert abs(m - 0.010) < 1e-9, (fam, m)
        assert n == 3
        assert res is True or se == 0.0


def test_s3_multiplicity_corrected_stiffness_prediction_is_the_product_it_names():
    """The second stiffness prediction is (d_eq/ds) * m_lo/m_hi, hand-computed here; the
    synthetic slopes satisfy the FIRST prediction exactly, so with m_hi = 1.1719 the
    measured ratio must sit at 1.1719x the corrected one -- and without `mult` the column
    must say nan rather than silently assume m = 1."""
    lines = []
    g, _ = grade_shape.g_s3(_tab(0.0), lambda *a, **k: None, mult={1.0: 1.0, 2.0: 1.1719})
    lines = [ln for ln in g.lines if ln.startswith("| lognormal | ")]
    row = [c.strip() for c in lines[-1].strip("|").split("|")]
    pred = 1.0 / grade_shape.sub_sphere_diameter(1.0, 2.0)
    assert abs(float(row[4]) - round(pred, 4)) < 1e-4
    assert abs(float(row[6]) - round(pred / 1.1719, 4)) < 1e-4, row
    assert abs(float(row[7]) - round(1.1719, 3)) < 2e-3, row
    assert any("right sign and the wrong size" in ln or "SIGN of the miss" in ln for ln in g.lines)
    g0, _ = grade_shape.g_s3(_tab(0.0), lambda *a, **k: None)
    row0 = [c.strip() for c in [ln for ln in g0.lines if ln.startswith("| lognormal | ")][-1]
            .strip("|").split("|")]
    assert row0[6] == "nan" and row0[7] == "nan", row0
    assert not any("SIGN of the miss" in ln for ln in g0.lines)


def test_s3_reports_not_resolved_when_the_interaction_is_zero_plus_noise():
    """The discriminating negative: a design that called everything resolved would pass
    the test above and fail this one."""
    g, s3 = grade_shape.g_s3(_tab(0.0, noise=0.002), lambda *a, **k: None)
    assert s3["inter"], "no families were graded"
    assert not any(v[3] for v in s3["inter"].values()), s3["inter"]


def test_s3_critical_value_is_t_at_n_minus_1_not_a_flat_two():
    """At three seeds the two-sided 95 % value is 4.303. A flat `2*se` guard is an 18 %
    test and is exactly what this project had to retract on the contact-force CV
    channel."""
    _g, s3 = grade_shape.g_s3(_tab(0.010, noise=0.001), lambda *a, **k: None)
    tcrits = {round(v[2], 3) for v in s3["inter"].values()}
    assert tcrits == {4.303}, tcrits


def test_s3_is_paired_so_a_common_per_seed_offset_cancels_exactly():
    """The whole reason for pairing on the seed: a size-realisation offset shared by
    both alpha rungs must leave the interaction untouched."""
    a = _tab(0.010)
    b = {k: {s: (v[0] + 0.03 * s, *v[1:]) for s, v in d.items()} for k, d in a.items()}
    _g1, s1 = grade_shape.g_s3(a, lambda *x, **k: None)
    _g2, s2 = grade_shape.g_s3(b, lambda *x, **k: None)
    for fam in s1["inter"]:
        assert abs(s1["inter"][fam][0] - s2["inter"][fam][0]) < 1e-12


def test_s3_needs_two_complete_seed_quartets():
    """A seed present at one alpha and missing at the other cannot enter a paired
    difference, and one surviving seed is not a standard error."""
    tab = _tab(0.010)
    for c in (0.0, 0.5):
        tab[("uniform", c, 2.0)] = {1: tab[("uniform", c, 2.0)][1]}
    _g, s3 = grade_shape.g_s3(tab, lambda *a, **k: None)
    assert "uniform" not in s3["inter"]
    assert "lognormal" in s3["inter"]


def test_s3_fails_when_the_grid_is_not_two_dimensional():
    """One alpha rung is not a shape axis; the gate must not report a verdict on it."""
    tab = _tab(0.0, alphas=(1.0,))
    g, _s3 = grade_shape.g_s3(tab, lambda *a, **k: None)
    assert g.passed is False


# --------------------------------------------------------------------------------------
# 7c. The registered pressure ladder (amendment 5.6)
# --------------------------------------------------------------------------------------

def _ladder_records(fam="lognormal", cv=0.0, alpha=1.0, seed=1, phats=None):
    phats = phats if phats is not None else grade_shape.P_HAT_LADDER_SHAPE
    out = []
    for p in phats:
        r = _prod_record(alpha, fam=fam, cv=cv, seed=seed, p_hat=p,
                         phi=0.62 - 5.0 * p)
        r["p_hat_target"] = p
        r["key"] = f"{fam}_{cv}_{alpha}_{seed}_{p:.3e}"
        out.append(r)
    return out


def test_off_ladder_rungs_are_excluded_by_design():
    """A 3e-5 record exists on disk for some cells and not others. It must be excluded
    UNIFORMLY, or one cell's phi_J has a longer lever arm than its neighbour's."""
    recs = _ladder_records() + _ladder_records(phats=(3e-5,))
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    assert len(tab) == 1
    (fitted,) = tab.values()
    assert list(fitted) == [1]


def test_a_cell_missing_a_registered_rung_is_dropped_not_fitted_short():
    """EXACTLY the registered rungs, not 'at least three'."""
    recs = _ladder_records(phats=(1e-3, 3e-4))
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    assert tab == {} or all(not v for v in tab.values())


def test_the_registered_ladder_is_the_one_the_prereg_names():
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert "`p_hat in {1e-3, 3e-4, 1e-4}`" in text
    assert grade_shape.P_HAT_LADDER_SHAPE == (1e-3, 3e-4, 1e-4)
    assert "z_all` = 0.023" in text, (
        "the prereg must carry the measured evidence for dropping the rung, not just "
        "the decision")


def test_ladder_filter_keeps_a_matched_grid():
    """Two cells, one of which happens to have a stray extra rung: after filtering both
    must be fitted from the same three pressures."""
    recs = (_ladder_records(cv=0.0) + _ladder_records(cv=0.5)
            + _ladder_records(cv=0.0, phats=(3e-5,)))
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    assert len(tab) == 2
    for v in tab.values():
        assert list(v) == [1]
        assert v[1][3] == 3, "every fitted ladder must have exactly three points"


# --------------------------------------------------------------------------------------
# 7d. S5 -- the CV = 0 identity and the aliasing it licenses
# --------------------------------------------------------------------------------------

def _cv0_records(fams=("lognormal", "uniform"), phi=0.5977, differ_at_cv_hi=True):
    out = []
    for fam in fams:
        for a in (1.0, 2.0):
            for s in (1, 2, 3):
                for p in grade_shape.P_HAT_LADDER_SHAPE:
                    r = _prod_record(a, fam=fam, cv=0.0, seed=s, p_hat=p,
                                     phi=phi - 5.0 * p)
                    r["p_hat_target"] = p
                    r["key"] = f"{fam}_0_{a}_{s}_{p:.3e}"
                    out.append(r)
                    q = _prod_record(a, fam=fam, cv=0.5, seed=s, p_hat=p,
                                     phi=phi + 0.06 - 5.0 * p
                                     + (0.01 if (differ_at_cv_hi and fam == "uniform")
                                        else 0.0))
                    q["p_hat_target"] = p
                    q["key"] = f"{fam}_5_{a}_{s}_{p:.3e}"
                    out.append(q)
    return out


def test_s5_passes_when_cv0_is_identical_and_cv_hi_differs():
    recs = _cv0_records()
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    g, n = grade_shape.g_s5(recs, tab)
    assert g.passed is True
    assert n == 0, "nothing to alias when both families ran CV = 0"


def test_s5_fails_when_cv0_differs_between_families():
    recs = _cv0_records()
    for r in recs:
        if r["family"] == "uniform" and r["cv"] == 0.0:
            r["phi"] += 1e-12
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    g, _n = grade_shape.g_s5(recs, tab)
    assert g.passed is False


def test_s5_fails_when_the_families_do_NOT_differ_at_cv_hi():
    """The discriminating half. If CV = 0.5 is also identical, the CV = 0 identity only
    says the family flag does nothing -- which is a bug, not a licence to alias."""
    recs = _cv0_records(differ_at_cv_hi=False)
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    g, _n = grade_shape.g_s5(recs, tab)
    assert g.passed is False


def test_s5_aliases_the_cv0_rung_into_a_family_that_lacks_it():
    """The saving: bidisperse ran no CV = 0 cell, and must still get a D."""
    recs = _cv0_records()
    recs += [r for r in _cv0_records(fams=("bidisperse",))
             if r["cv"] != 0.0]
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    assert ("bidisperse", 0.0, 1.0) not in tab
    g, n = grade_shape.g_s5(recs, tab)
    assert g.passed is True
    assert n == 6, n
    assert ("bidisperse", 0.0, 1.0) in tab
    _g3, s3 = grade_shape.g_s3(tab, lambda *a, **k: None)
    assert "bidisperse" in s3["inter"]


def test_s5_alias_does_not_overwrite_a_family_that_ran_its_own():
    recs = _cv0_records()
    tab = grade_shape.phi_j_table(recs, lambda *a, **k: None)
    before = dict(tab[("uniform", 0.0, 1.0)])
    grade_shape.g_s5(recs, tab)
    assert tab[("uniform", 0.0, 1.0)] == before


def _clearance_failure(tag):
    return {"failed": True, "tag": tag, "alpha": 2.0, "cv": 0.5, "seed": 1,
            "p_hat_target": 1e-3, "n_requested": 1700, "family": "lognormal",
            "error": ("ValueError('cell too small for the periodic collider at cv=0.5, "
                      "N=1700, alpha=2.0 ... Use N >= 3017')")}


def test_s2_does_not_count_a_probe_cell_as_a_production_failure(tmp_path):
    """The clearance guard firing on a deliberately undersized PROBE cell is the guard
    working. Counting it against S2 reported a working guard as a failed gate, which is
    what the first version of this file did."""
    f = tmp_path / "phase2_shp_x.jsonl"
    rows = [_clearance_failure("shptime"), _prod_record(2.0)]
    f.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    g = grade_shape.g_s2([str(f)], [_prod_record(2.0)])
    assert g.passed is True


def test_s2_still_fails_on_a_real_production_clearance_error(tmp_path):
    f = tmp_path / "phase2_shp_y.jsonl"
    rows = [_clearance_failure("shp"), _prod_record(2.0)]
    f.write_text("\n".join(json.dumps(r) for r in rows) + "\n", encoding="utf-8")
    g = grade_shape.g_s2([str(f)], [_prod_record(2.0)])
    assert g.passed is False


def test_s2_fails_when_N_is_not_constant_across_the_grid():
    """phi_J drifts with N, so an N that moved with alpha would wear alpha's label."""
    recs = [_prod_record(1.0), dict(_prod_record(2.0), n_requested=1700)]
    g = grade_shape.g_s2([], recs)
    assert g.passed is False


def test_the_shared_reference_makes_a_cross_family_difference_CONSERVATIVE():
    """The sign of the CV = 0 correlation, DERIVED rather than asserted.

    The first version of prereg 5.8 said a cross-family comparison "understates its own
    uncertainty" if it treats the three families as independent. It OVERSTATES it: with
    D_A = X_A - C and D_B = X_B - C the shared C cancels exactly in the difference, so
    the naive independent variance adds 2*Var(C) that is not there. Pinned here because
    the claim is a sentence in a prereg and a sentence can drift back.
    """
    rng = np.random.default_rng(0)
    n = 200_000
    C = rng.normal(size=n)          # the shared CV = 0 reference
    X_A = rng.normal(size=n)
    X_B = rng.normal(size=n)
    D_A, D_B = X_A - C, X_B - C
    true_var = float(np.var(D_A - D_B, ddof=1))
    naive_var = float(np.var(D_A, ddof=1) + np.var(D_B, ddof=1))
    assert true_var < naive_var, (true_var, naive_var)
    assert abs(true_var - 2.0) < 0.05, true_var      # Var(X_A) + Var(X_B)
    assert abs(naive_var - 4.0) < 0.05, naive_var    # + 2*Var(C)


def test_prereg_records_the_corrected_sign_and_the_correction():
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert "HAD THE SIGN BACKWARDS" in text, (
        "the correction must stay visible; a silently replaced claim is how the wrong "
        "one comes back")
    assert "assuming independence is" in text or "It overstates it." in text


def test_cv0_psd_really_is_family_independent():
    """The claim S5 rests on, checked against `psd.make_psd` directly rather than
    against the records -- and with the discriminating negative at CV = 0.5."""
    from psd import make_psd
    ref = None
    for fam in ("lognormal", "uniform", "bidisperse"):
        rng = np.random.default_rng(1)
        r = make_psd(fam, 1.0, 0.0).radii_for_yade(500, rng, truncate=(0.2, 5.0))
        if ref is None:
            ref = r
        else:
            assert np.array_equal(r, ref)
    ref = None
    for fam in ("lognormal", "uniform", "bidisperse"):
        rng = np.random.default_rng(1)
        r = make_psd(fam, 1.0, 0.5).radii_for_yade(500, rng, truncate=(0.2, 5.0))
        if ref is None:
            ref = r
        else:
            assert not np.array_equal(r, ref), (
                f"{fam} at CV = 0.5 is identical to lognormal -- the family flag is "
                f"doing nothing, and the CV = 0 identity would be meaningless")


# --------------------------------------------------------------------------------------
# 8. The prereg amendment is on disk and says what the code does
# --------------------------------------------------------------------------------------

def test_prereg_amendment_records_the_alpha_1_seam():
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert "## 5. AMENDMENT" in text
    for token in ("SINGLE SPHERE", "0.7500", "0.7071", "S0a", "S0b", "S0c",
                  "contact multiplicity", "5379"):
        assert token in text, f"the amendment does not mention {token!r}"


def test_prereg_amendment_predates_any_production_cell():
    """The amendment's own provenance claim, checked rather than trusted: no record with
    a production tag may be older than the amendment. Since the amendment is committed
    with the runner, the checkable form is that the amendment names the gate tags it
    excludes."""
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert "shpgate" in text and "shptime" in text
    assert set(grade_shape.GATE_TAGS) == {"shpgate", "shptime"}


def test_span_factor_cubed_is_the_clearance_scaling_in_the_prereg():
    """The N-scaling table in amendment 5.4 is span_factor^3, recomputed here."""
    expect = [1.0000, 1.4285, 2.0000, 2.8005, 4.0000]
    got = [shape.span_factor(a) ** 3 for a in shape.ALPHA_LADDER]
    assert np.allclose(got, expect, atol=5e-4), got


def test_lognormal_cv05_alpha2_needs_the_N_the_prereg_reports():
    """5379 is the binding number in the whole design; recomputed from the PSD."""
    from psd import make_psd
    worst = 0
    for seed in range(1, 13):
        rng = np.random.default_rng(seed)
        r = make_psd("lognormal", 1.0, 0.5).radii_for_yade(1500, rng,
                                                           truncate=(0.2, 5.0))
        vs = float(np.sum(4.0 / 3.0 * math.pi * r ** 3))
        dmax = float(2.0 * r.max()) * shape.span_factor(2.0)
        worst = max(worst, math.ceil(0.72 * (2.5 * dmax) ** 3 / (vs / len(r))))
    assert worst == 5379, worst


def test_runner_and_grader_are_importable_without_yade():
    """Both must import on Windows: the grader runs there and the runner's argument
    parsing is tested there."""
    assert hasattr(grade_shape, "main")
    sys.modules.pop("phase2_shape", None)
    import phase2_shape
    assert phase2_shape.SHAPE_TAG == "shp"
    assert runpy is not None


# ---------------------------------------------------------------------------------
# Amendment 5.10: a ladder that mixes caps, and gate S6
# ---------------------------------------------------------------------------------

def _cap_record(cap, converged, iters, phi=0.5977, cv=0.0, seed=1, alpha=1.0,
                fam="lognormal", p_hat=1e-3):
    r = _prod_record(alpha, fam=fam, cv=cv, seed=seed, p_hat=p_hat, phi=phi)
    base = f"cv{cv:.3f}_mu0.300_s{seed}_p{p_hat:.3e}_N5400_a{alpha:.3f}"
    if fam != "lognormal":
        base += f"_{fam}"
    r["key"] = base if cap == grade_shape.DEFAULT_CAP else f"{base}_it{cap}"
    r["max_iter"] = cap
    r["converged"] = converged
    r["iters"] = iters
    r["tol_met"] = converged
    return r


def test_physical_key_strips_only_the_cap_suffix():
    d = _cap_record(200000, True, 19396)
    c = _cap_record(1000000, True, 19396)
    assert d["key"] != c["key"]
    assert grade_shape.physical_key(d) == grade_shape.physical_key(c)
    u = _cap_record(1000000, True, 19396, fam="uniform")
    assert grade_shape.physical_key(u).endswith("_uniform")


def test_prefer_raised_cap_uses_the_raised_record_and_prints_it():
    d = _cap_record(200000, False, 200000, phi=0.6212, cv=0.5)
    c = _cap_record(1000000, True, 237746, phi=0.6185, cv=0.5)
    lines = []
    out, multi = grade_shape.prefer_raised_cap([d, c], lines.append)
    assert len(out) == 1 and out[0]["max_iter"] == 1000000 and out[0]["phi"] == 0.6185
    assert len(multi) == 1
    assert any("using 1000000" in s for s in lines)


def test_prefer_raised_cap_leaves_a_single_cap_cell_alone():
    d = _cap_record(200000, True, 19396)
    out, multi = grade_shape.prefer_raised_cap([d], lambda *a, **k: None)
    assert out == [d] and multi == {}


def test_prefer_raised_cap_refuses_two_records_at_one_cap():
    with pytest.raises(SystemExit):
        grade_shape.prefer_raised_cap([_cap_record(200000, True, 1),
                                       _cap_record(200000, True, 1)],
                                      lambda *a, **k: None)


def _pair(phi_d=0.5977, phi_c=0.5977, iters_d=19396, iters_c=19396, z_c=4.0,
          conv_d=True):
    d = _cap_record(200000, conv_d, iters_d, phi=phi_d)
    c = _cap_record(1000000, True, iters_c, phi=phi_c)
    c["z_all"] = z_c
    return {grade_shape.physical_key(d): {200000: d, 1000000: c}}


def test_s6_passes_on_a_bit_identical_converged_pair():
    g = grade_shape.g_s6(_pair())
    assert g.passed is True and g.informative


def test_s6_fails_on_one_ulp_in_phi():
    g = grade_shape.g_s6(_pair(phi_c=np.nextafter(0.5977, 1.0)))
    assert g.passed is False and g.informative
    assert any("DIFFER on phi" in s for s in g.lines)


def test_s6_fails_on_iters_alone():
    g = grade_shape.g_s6(_pair(iters_c=19397))
    assert g.passed is False


def test_s6_fails_on_z_all_alone():
    g = grade_shape.g_s6(_pair(z_c=4.0000001))
    assert g.passed is False


def test_s6_is_vacuous_with_no_pair_and_vacuous_is_not_a_pass():
    g = grade_shape.g_s6({})
    assert g.passed is False and not g.informative
    assert "NOT INFORMATIVE" in g.render()[0]


def test_s6_is_vacuous_when_the_only_pairs_are_cap_cut():
    """A cap-cut cell carries no identity claim, so a re-run with no converged pair
    cannot show the protocol is the same -- VACUOUS, not PASS."""
    m = _pair(phi_d=0.6212, phi_c=0.6185, iters_d=200000, iters_c=237746, conv_d=False)
    g = grade_shape.g_s6(m)
    assert g.passed is False and not g.informative


def test_s6_cap_cut_pair_is_not_required_to_be_identical_but_must_run_past_the_cap():
    m = _pair()
    cut = _pair(phi_d=0.6212, phi_c=0.6185, iters_d=200000, iters_c=237746, conv_d=False)
    for k, v in cut.items():
        for r in v.values():
            r["cv"] = 0.5
        m[k + "_cut"] = v
    g = grade_shape.g_s6(m)
    assert g.passed is True
    assert any("1/1 ran past the old cap" in s for s in g.lines)
    # and the discriminating half: a 'raised-cap' record that stopped AT the old cap
    # means the cap was never in force for it
    bad = _pair(phi_d=0.6212, phi_c=0.6212, iters_d=200000, iters_c=200000, conv_d=False)
    for k, v in bad.items():
        m[k + "_bad"] = v
    g = grade_shape.g_s6(m)
    assert g.passed is False


def test_s6_sentence_names_the_registered_fields():
    g = grade_shape.g_s6(_pair())
    for f in ("phi", "iters", "z_all", "p_hat"):
        assert f in g.sentence
    assert grade_shape.S6_FIELDS == ("phi", "iters", "z_all", "p_hat")


def test_prereg_5_10_registers_the_cap_the_launcher_uses():
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert "### 5.10" in text
    assert "`PJ_MAX_ITER = 1000000`" in text
    assert "must be BIT-IDENTICAL to the default-cap records" in text
    assert "the raised-cap record is used and the fact printed" in text
    with open(os.path.join(HERE, "scripts", "shape_rerun.sh"), encoding="utf-8") as fh:
        sh = fh.read()
    assert "export PJ_MAX_ITER=1000000" in sh
    assert "TAG=shpcap" in sh
    assert "shpcap" not in grade_shape.GATE_TAGS, "shpcap is PRODUCTION, not a probe"


def test_prereg_5_11_records_the_seed3_abort_from_the_probe_files():
    """5.11 quotes the crash probe and its controls; every number it quotes must be the one
    on disk, and the missing cell must really be missing from the production records."""
    with open(os.path.join(HERE, "PREREG_SHAPE_2026-09-18.md"), encoding="utf-8") as fh:
        text = fh.read()
    assert text.count("### 5.11") == 1
    assert "No replacement seed is added" in text
    assert "t_crit = 12.706" in text
    crash = open(os.path.join(HERE, "out", "shape", "probe_seed3_crash.txt"), encoding="utf-8-sig").read()
    ctrl = open(os.path.join(HERE, "out", "shape", "probe_seed3_controls.txt"), encoding="utf-8").read()
    n_inside = re.search(r"N: (\d+) cross-grain pairs with a centre inside", crash).group(1)
    assert f"this one has {n_inside} cross-grain pairs" in text
    assert "a body would cross a third of the cell" in crash, "the probe must reach the abort"
    ctrl_inside = [int(x) for x in re.findall(r"N: (\d+) cross-grain pairs with a centre inside", ctrl)]
    assert len(ctrl_inside) == 7
    assert f"Centre-inside pairs {min(ctrl_inside)}-{max(ctrl_inside)}" in text
    assert min(ctrl_inside) <= int(n_inside) <= max(ctrl_inside), "5.11 says the census does not discriminate"
    assert "a body would cross" not in ctrl, "a control that aborted would refute 5.11"
    # the cell 5.11 says is missing IS missing, and only that one, in the lognormal alpha-2 CV-0.5 rung
    prod, _, _ = grade_shape.load_records(report=lambda *a, **k: None)
    seeds = sorted({r["seed"] for r in prod if r["family"] == "lognormal"
                    and abs(r["alpha"] - 2.0) < 1e-9 and abs(r["cv"] - 0.5) < 1e-9})
    assert seeds == [1, 2], seeds


def test_phi_j_table_records_the_caps_of_every_fitted_ladder():
    recs = []
    for p in grade_shape.P_HAT_LADDER_SHAPE:
        cap = 1000000 if p < 1e-3 else 200000
        r = _cap_record(cap, True, 1000, phi=0.62 - 5.0 * p, cv=0.5, p_hat=p)
        recs.append(r)
    lines = []
    tab = grade_shape.phi_j_table(recs, lines.append)
    (seeds,) = tab.values()
    assert seeds[1][5] == (200000, 1000000)
    assert any("[200000, 1000000]: 1" in s for s in lines)

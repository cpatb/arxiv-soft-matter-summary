"""test_paper_numbers.py -- the manuscript quotes nothing it did not print.

Five checks, each of which can fail:

  1. `extract_numbers.py` regenerates `numbers.tex` and `numbers.json` BYTE-IDENTICALLY
     (the committed macros are what the script prints today, not a stale copy). The
     artifacts are restored before the assertion so a failure does not destroy them.
  2. `main.tex` carries no typed result. After stripping comments and the arguments of
     \\cite, \\ref, \\label, \\includegraphics and the preamble, every numeral left must
     be one of a short allowlist of DESIGN TOKENS (the CV grid, the friction value, the
     "d_10" subscript, "95 %") and no numeral anywhere may have three or more
     significant digits. A stray "0.0621" in the prose fails here.
  3. Three load-bearing macros equal the value the ARTIFACT prints, read with a second,
     independent regex (not the one in extract_numbers.py).
  4. The figures regenerate byte-identically (`make_paper_figures.py --check`), so the
     PDFs the manuscript includes are what the scripts draw from the committed shards.
  5. An analysed "cell" is a cell and not a record, and every `*N` macro is a seed count.
     Added 2026-09-17 after a raised-cap re-run was found being counted as an extra seed.

Run from anywhere:  python -m pytest paper/test_paper_numbers.py -q
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")
MAIN_TEX = os.path.join(HERE, "main.tex")
NUMBERS_TEX = os.path.join(HERE, "numbers.tex")
NUMBERS_JSON = os.path.join(HERE, "numbers.json")
FIG = os.path.join(HERE, "figures")

# Numerals allowed to be typed into main.tex. Each is a DESIGN token, not a result:
# the CV grid end points and step, the mu_p and rattler-k values written as prose
# anchors ("k = 2"), the "10" of d_{10}, the 95 of "95 %", and small counting
# words that appear as numerals ("three", "two" are spelled out; "1" and "2" occur in
# equation numbering and "(1 +- CV)"). Adding to this list is a manual, reviewed act.
ALLOWED_NUMERALS = {"0", "0.1", "0.3", "0.5", "1", "2", "3", "10", "95"}


def _md5(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def _read(path: str) -> str:
    with open(path, encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------------------------
# 1. the macros regenerate byte-identically
# ---------------------------------------------------------------------------------
def test_numbers_tex_regenerates_byte_identically():
    assert os.path.exists(NUMBERS_TEX), "paper/numbers.tex is committed input; it is missing"
    assert os.path.exists(NUMBERS_JSON)
    before = {p: _md5(p) for p in (NUMBERS_TEX, NUMBERS_JSON)}
    keep = tempfile.mkdtemp()
    for p in before:
        shutil.copy2(p, os.path.join(keep, os.path.basename(p)))
    try:
        r = subprocess.run([sys.executable, os.path.join(HERE, "extract_numbers.py")],
                           capture_output=True, text=True, cwd=ROOT)
        after = {p: _md5(p) for p in before}
    finally:
        # restore FIRST, then judge: a failing regeneration must not eat the reference
        for p in before:
            shutil.copy2(os.path.join(keep, os.path.basename(p)), p)
        shutil.rmtree(keep, ignore_errors=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "all recomputed values matched their artifacts" in r.stdout
    assert after == before, "extract_numbers.py no longer reproduces the committed numbers.tex/json"


# ---------------------------------------------------------------------------------
# 2. no typed result digits in main.tex
# ---------------------------------------------------------------------------------
def _stripped_main_tex() -> str:
    src = _read(MAIN_TEX)
    # preamble is design (class options, package versions); results start at \begin{document}
    assert src.count(r"\begin{document}") == 1
    body = src.split(r"\begin{document}", 1)[1]
    body = re.sub(r"(?<!\\)%.*", "", body)                       # comments (not \%)
    for cmd in ("cite", "ref", "label", "includegraphics", "bibliography", "input",
                "bibliographystyle", "hypersetup"):
        body = re.sub(rf"\\{cmd}\*?(\[[^\]]*\])?\{{[^}}]*\}}", f"\\\\{cmd}{{}}", body)
    body = re.sub(r"\\date\{[^}]*\}", "", body)
    return body


def _numerals(body: str) -> list[str]:
    # `\to0.5` is a macro name followed by the numeral 0.5: separate every control word
    # from what follows it, so the numeral is read whole rather than from its second digit
    body = re.sub(r"(\\[A-Za-z]+)", r"\1 ", body)
    return re.findall(r"(?<![A-Za-z{0-9.])[0-9]+(?:\.[0-9]+)?", body)


def _offending(numerals: list[str]) -> tuple[list[str], list[str]]:
    outside = sorted({n for n in numerals if n not in ALLOWED_NUMERALS})
    sig3 = [n for n in numerals if len(n.replace(".", "").lstrip("0")) >= 3]
    return outside, sig3


def test_main_tex_has_no_typed_result_digits():
    numerals = _numerals(_stripped_main_tex())
    assert numerals, "the regex found no numerals at all; the stripping is over-eager"
    outside, sig3 = _offending(numerals)
    assert not outside, f"typed numerals in main.tex outside the design allowlist: {outside}"
    # unconditional: nothing with >= 3 significant digits, allowlist or not
    assert not sig3, f"numerals with three or more significant digits typed into main.tex: {sig3}"


def test_digit_check_catches_a_typed_result():
    # the discriminating negative: the same scan on the manuscript with one result typed
    # into the prose must fail on BOTH rules, or the positive test above proves nothing
    body = _stripped_main_tex().replace(r"\dPhiJLogn\pm", r"+0.0621\pm", 1)
    assert "+0.0621" in body
    outside, sig3 = _offending(_numerals(body))
    assert outside == ["0.0621"], outside
    assert sig3 == ["0.0621"], sig3
    # and a two-significant-digit stray that is not a design token is also caught
    outside, _ = _offending(_numerals(_stripped_main_tex() + " a ratio of 4.5 "))
    assert outside == ["4.5"], outside


def test_every_macro_used_in_main_tex_is_defined_and_no_digit_in_macro_names():
    macros = dict(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{([^}]*(?:\{[^}]*\}[^}]*)*)\}",
                             _read(NUMBERS_TEX)))
    assert len(macros) >= 200, f"only {len(macros)} macros parsed from numbers.tex"
    assert all(re.fullmatch(r"[A-Za-z]+", k) for k in macros)
    body = _read(MAIN_TEX)
    used = set(re.findall(r"\\([A-Za-z]+)", body))
    # the manuscript must actually draw on the macro file, and its symbol macros must not
    # shadow a value macro (\muP = 0.3 is a value; the symbol is \mup)
    own = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", body))
    assert not own & set(macros), f"main.tex redefines numbers.tex macros: {sorted(own & set(macros))}"
    # 127 of 234 at the 2026-09-16 draft; the unused ones are per-rung table values
    assert len(used & set(macros)) >= 100, f"main.tex uses only {len(used & set(macros))} macros"
    # a typo in a macro name is an undefined control sequence at build time; if a log is
    # present it must be clean (the build is the authoritative check, this reads its verdict)
    log = os.path.join(HERE, "main.log")
    if os.path.exists(log):
        with open(log, encoding="utf-8", errors="replace") as fh:
            txt = fh.read()
        assert "Undefined control sequence" not in txt
        assert "Citation" not in txt or "undefined" not in txt.split("Citation", 1)[1][:60]


# ---------------------------------------------------------------------------------
# 3. load-bearing macros equal the artifact, by an independent regex
# ---------------------------------------------------------------------------------
def _macro(name: str) -> str:
    m = re.findall(rf"\\newcommand\{{\\{name}\}}\{{([^}}]*)\}}", _read(NUMBERS_TEX))
    assert len(m) == 1, (name, m)
    return m[0]


def test_headline_shift_equals_the_sweep_artifact():
    txt = _read(os.path.join(OUT, "sweep_descriptive.txt"))
    m = re.findall(r"^\s+lognormal\s+([+-]\d\.\d{6})\s+(\d\.\d{6})\s+3\s+([+-]\d\.\d{4})\s+YES$", txt, re.M)
    assert len(m) == 1, m
    recomputed, se, committed = m[0]
    assert _macro("dPhiJLogn") == committed
    assert _macro("dPhiJLogn") == f"{float(recomputed):+.4f}"
    assert _macro("dPhiJLognSe") == f"{float(se):.4f}"
    assert json.load(open(NUMBERS_JSON, encoding="utf-8"))["dPhiJLogn"] == committed


def test_matched_pair_equals_the_mixture_artifact():
    txt = _read(os.path.join(OUT, "mixture.txt"))
    m = re.findall(r"^\s+MEASURED\s+dphi = ([+-]\d\.\d{4}) \+- (\d\.\d{4})\s+\(([+-]\d+\.\d) scale", txt, re.M)
    assert len(m) == 1, m
    dphi, se, scale = m[0]
    assert _macro("pairDPhi") == dphi
    assert _macro("pairDPhiSe") == se
    assert _macro("pairScale") == scale.lstrip("+")
    p = re.findall(r"filling account ([+-]\d\.\d{4})$", txt, re.M)
    assert len(p) == 1, p
    assert _macro("pairFillingPred") == p[0]


def test_fabric_slope_equals_the_primary_k2_block():
    txt = _read(os.path.join(OUT, "grade_fabric.txt"))
    # the artifact has TWO slope tables (k = 2 PRIMARY, k = 3 reported); the paper's
    # macro must be the PRIMARY one, and the two must differ so the check can fail
    blocks = re.findall(r"T_nr at k = (\d)\s+\((?:PRIMARY|reported)\)\n.*\n\s+lognormal\s+([+-]\d\.\d{4})\s+(\d\.\d{4})",
                        txt)
    assert [b[0] for b in blocks] == ["2", "3"], blocks
    k2, k3 = blocks
    assert k2[1] != k3[1], "k = 2 and k = 3 slopes coincide; this test could not discriminate"
    assert _macro("fabSlopeLogn") == f"{float(k2[1]):+.3f}"
    assert _macro("fabSlopeLognSe") == f"{float(k2[2]):.3f}"


def test_pair_skewness_macros_match_a_direct_moment_calculation():
    # anchor OUTSIDE the closed form: build the two-point number distribution and take
    # its third standardised moment directly (numpy), then compare with the macros
    sys.path.insert(0, ROOT)
    sys.path.insert(0, HERE)
    import numpy as np
    from psd import BinaryPSD
    from make_paper_figures import binary_number_skewness
    for xs, macro in ((1 / 28, "pairSkewA"), (0.25, "pairSkewB")):
        p = BinaryPSD(d10=1.0, ratio=3.0, x_small=xs)
        w = np.array([p.n_small, 1.0 - p.n_small])
        d = np.array([p.d_small, p.d_large])
        mu = (w * d).sum()
        m2 = (w * (d - mu) ** 2).sum()
        m3 = (w * (d - mu) ** 3).sum()
        direct = m3 / m2 ** 1.5
        assert abs(binary_number_skewness(3.0, xs) - direct) < 1e-12
        assert float(_macro(macro)) == round(direct, 2), (macro, direct)
    # the discriminating half: B is skewed and positive (its fines are the majority),
    # A is symmetric; a routine returning 0 for both would pass neither line below
    assert float(_macro("pairSkewA")) == 0.0
    assert float(_macro("pairSkewB")) > 2.0
    # hand check from the summary: n_s = 0.9 at d = 1 and 0.1 at d = 3 gives mean 1.2,
    # variance 0.36 and third moment 0.576, i.e. skewness 0.576 / 0.216
    assert abs(binary_number_skewness(3.0, 0.25) - 0.576 / 0.216) < 1e-9


def test_dropped_mixture_cell_macros_match_the_gate():
    sys.path.insert(0, ROOT)
    from phase12_mixture import MIN_PHASE_COUNT, N_PROD
    from psd import BinaryPSD
    assert _macro("mixMinPhaseCount") == str(MIN_PHASE_COUNT)
    assert _macro("mixNProd") == str(N_PROD)
    r, xs = float(_macro("mixDroppedRatio")), float(_macro("mixDroppedXs"))
    p = BinaryPSD(d10=1.0, ratio=r, x_small=xs)
    n_minor = min(p.n_small, 1 - p.n_small) * N_PROD
    assert n_minor < MIN_PHASE_COUNT, n_minor
    assert _macro("mixDroppedNminor") == f"{n_minor:.1f}"
    # the manuscript's sentence must agree with the figure's own sidecar
    md = _read(os.path.join(FIG, "fig3_binaries.md"))
    assert f"`X_s = {xs:.1f}` is omitted" in md
    assert f"({n_minor:.1f} expected large grains against a floor of {MIN_PHASE_COUNT})" in md


# ---------------------------------------------------------------------------------
# 4. figures regenerate byte-identically
# ---------------------------------------------------------------------------------

def test_figures_regenerate_byte_identically():
    """⚠ 2026-09-19: this used to save/restore and compare only the PDFs, so a stale
    committed `.md` caption sidecar or `.csv` passed and was then silently overwritten by
    this very test -- leaving a dirty working tree as the only signal. `save()` writes
    four files per figure and all four are now saved, compared and restored."""
    exts = (".pdf", ".png", ".csv", ".md")
    arts = sorted(p for p in os.listdir(FIG) if p.endswith(exts))
    expect = sorted(f"{s}{e}" for s in ("fig1_phiJ", "fig2_backbone", "fig3_binaries")
                    for e in exts)
    assert arts == expect, arts
    before = {p: _md5(os.path.join(FIG, p)) for p in arts}
    keep = tempfile.mkdtemp()
    for p in arts:
        shutil.copy2(os.path.join(FIG, p), os.path.join(keep, p))
    try:
        r = subprocess.run([sys.executable, os.path.join(HERE, "make_paper_figures.py"), "--check"],
                           capture_output=True, text=True, cwd=ROOT)
        after = {p: _md5(os.path.join(FIG, p)) for p in arts}
    finally:
        for p in arts:
            shutil.copy2(os.path.join(keep, p), os.path.join(FIG, p))
        shutil.rmtree(keep, ignore_errors=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert "byte-identity: PASS" in r.stdout
    changed = sorted(p for p in arts if after[p] != before[p])
    assert not changed, (
        f"the committed figure artifacts are not what make_paper_figures.py draws today: "
        f"{changed}")


# ---------------------------------------------------------------------------------
# 5. REGRESSION, 2026-09-17: an analysed "cell" is a cell, not a record
# ---------------------------------------------------------------------------------
# The shard record key carries the iteration cap (`phase2_phiJ.py` appends `_it{cap}` at
# any non-default MAX_ITER, deliberately, so a raised-cap run can never POOL with a
# default-cap one). `sweep_descriptive.load_shards` deduplicated on that key, so a cell
# run under both caps survived TWICE and every standard error over it was sqrt(n_dup)
# too small -- 35 of 252 records were one cell counted twice, and six published standard
# errors were 1.581x too small while six `*N` macros read 6 for 3 seeds.
#
# These two tests fail if that regresses. Note WHY the old suite could not catch it: the
# duplicated records are bit-identical, so a mean-based assertion stayed green. ⚠ But
# "no MEAN moved", which stood here until 2026-09-19, is wrong as stated: the duplication
# is balanced (2:2:2) in 11 of the 12 affected groups and UNBALANCED in
# `(uniform, CV 0.4, p_hat = 3e-05)`, 5 records over 3 seeds, where four means did move
# (`out/sweep_descriptive.txt` prints them). No published macro reads that group, because
# every published mean is at the p_hat = 1e-3 rung. So the general lesson stands -- a
# defect in the second moment needs a test on the second moment -- but a mean-based test
# on THAT group could have caught this one.

def _selected_cells():
    sys.path.insert(0, ROOT)
    import sweep_descriptive as sd
    sd.L.clear()
    return sd.load_shards()


def test_load_shards_returns_one_record_per_physical_cell():
    import collections
    sel = _selected_cells()
    assert sel, "load_shards returned nothing; the test below would pass vacuously"
    counts = collections.Counter(
        ((r.get("family") or "lognormal"), r["cv"], r["seed"], r["p_hat_target"],
         r.get("n_requested"), r.get("mu_p"))
        for r in sel)
    dup = {k: n for k, n in counts.items() if n > 1}
    assert not dup, (
        f"{len(dup)} physical cells appear more than once in load_shards' output; a mean "
        f"or standard error over these is not a mean over seeds. Examples: "
        f"{sorted(dup.items())[:3]}")
    # and the collapse must have had something to do -- if this ever becomes an equality
    # the raised-cap records have left the tree and the guard above is no longer testing
    # the thing it was written for.
    # ⚠ 2026-09-19: this used to parse n_raw out of the COMMITTED out/sweep_descriptive.txt
    # and compare it with a LIVE load. Add or remove a shard without regenerating the
    # artifact and the two sides describe different populations, so the guard could fire
    # (or pass) for a reason unrelated to the collapse. Both sides now come from the same
    # in-process load; the artifact is checked separately, as an artifact.
    import sweep_descriptive as sd
    n_raw = sd.LAST_LOAD["n_selected_records"]
    assert sd.LAST_LOAD["n_cells"] == len(sel), "LAST_LOAD disagrees with the returned list"
    assert n_raw > len(sel), (
        f"no record was collapsed ({n_raw} == {len(sel)}); either the raised-cap records "
        f"are gone, in which case delete this guard, or the collapse stopped running")
    m = re.search(r"N after N=1500, mu_p=0\.3, converged: (\d+)",
                  open(os.path.join(OUT, "sweep_descriptive.txt"), encoding="utf-8").read())
    assert m is not None, "the record-count line is no longer in out/sweep_descriptive.txt"
    assert int(m.group(1)) == n_raw, (
        f"the committed artifact says {m.group(1)} records and this load sees {n_raw}; "
        "out/sweep_descriptive.txt is stale -- regenerate it")


def test_every_seed_count_macro_is_a_seed_count():
    """Every `*N` macro equals the number of distinct seeds behind it, AND the group it
    reads holds one record per seed.

    ⚠ 2026-09-19. The first half alone could not detect the regression it names: `seeds`
    is a SET taken from the already-collapsed `_selected_cells()`, so it is invariant
    under duplication and re-introducing the `load_shards` defect would leave this test
    comparing `3 == 3`. The discriminating quantity is the RECORD count in the group
    against the seed-set size, which is checked below on the same in-process load.
    """
    sel = _selected_cells()
    numbers = json.load(open(NUMBERS_JSON, encoding="utf-8"))
    fams = {"Logn": "lognormal", "Unif": "uniform", "Bidi": "bidisperse"}
    cvs = {"Zero": 0.0, "Three": 0.3, "Half": 0.5}
    checked = 0
    for key in ("rat", "znr", "zall"):
        for tag, fam in fams.items():
            for ctag, cv in cvs.items():
                name = f"{key}{tag}{ctag}N"
                if name not in numbers:
                    continue
                group = [r for r in sel
                         if (r.get("family") or "lognormal") == fam and r["cv"] == cv
                         and r["p_hat_target"] == 1e-3]
                seeds = {r["seed"] for r in group}
                assert int(numbers[name]) == len(seeds), (
                    f"{name} = {numbers[name]} but only {len(seeds)} distinct seeds "
                    f"({sorted(seeds)}) stand behind it")
                # the half that can actually fail: one record per seed in the group
                assert len(group) == len(seeds), (
                    f"{name}'s group holds {len(group)} records over {len(seeds)} seeds; "
                    f"a mean over it weights seed {sorted(s for s in seeds if sum(1 for r in group if r['seed'] == s) > 1)} "
                    "more than the others")
                checked += 1
    assert checked == 27, f"expected 27 seed-count macros, checked {checked}"

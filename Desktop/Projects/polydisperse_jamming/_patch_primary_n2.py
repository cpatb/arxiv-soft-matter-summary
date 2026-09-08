"""One-shot: restrict the remaining three N-blind consumers to PRIMARY_N."""
from pathlib import Path

HERE = Path(__file__).resolve().parent

FAMILIES_OLD = '''    with open(fp, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:'''
FAMILIES_NEW = '''    with open(fp, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    # C3 compares FAMILIES at fixed CV. phase2_summary.csv gained an n_requested column
    # on 2026-09-08, so without this the same (family, cv, mu_p) appears twice at two
    # system sizes and the family comparison silently spans both. The N = 3000 coverage
    # does not span the CV range this comparison needs, so the size is fixed explicitly.
    n_before = len(rows)
    rows = [r for r in rows if int(float(r.get("n_requested", PRIMARY_N))) == PRIMARY_N]
    print(f"N RESTRICTION (C3 family comparison at fixed CV): keeping N = {PRIMARY_N}; "
          f"{len(rows)} of {n_before} summary rows, {n_before - len(rows)} dropped.")
    for r in rows:'''

DIAG_OLD = '''                if r.get("key") and not r.get("failed"):'''
DIAG_NEW = '''                if (r.get("key") and not r.get("failed")
                        and r.get("n_requested") == PRIMARY_N):'''

DIAG_NOTE = '''
# System size is fixed here on purpose. out/ holds N = 1500 and N = 3000 cells; this
# diagnostic asks about the CONVERGENCE / LADDER structure at one size, and pooling two
# would make a system-size difference read as ladder scatter. See N_GROUPING_2026-09-08.md.
from analyse_phase2 import PRIMARY_N  # noqa: E402
'''


def patch(name: str, subs, add_import: bool) -> None:
    p = HERE / name
    s = p.read_text(encoding="utf-8")
    for old, new in subs:
        n = s.count(old)
        assert n == 1, f"{name}: expected 1 match, found {n} for {old[:60]!r}"
        s = s.replace(old, new, 1)
    if add_import and "from analyse_phase2 import PRIMARY_N" not in s:
        marker = "def load("
        assert s.count(marker) >= 1, name
        i = s.index(marker)
        s = s[:i] + DIAG_NOTE.strip() + "\n\n\n" + s[i:]
    p.write_text(s, encoding="utf-8")
    print(f"patched {name}")


patch("analyse_families.py", [(FAMILIES_OLD, FAMILIES_NEW)], add_import=False)
# analyse_families imports nothing from analyse_phase2 yet
fp = HERE / "analyse_families.py"
s = fp.read_text(encoding="utf-8")
if "from analyse_phase2 import PRIMARY_N" not in s:
    anchor = "from psd import lognormal_skewness, make_psd"
    assert s.count(anchor) == 1
    s = s.replace(anchor, anchor + "\nfrom analyse_phase2 import PRIMARY_N", 1)
    fp.write_text(s, encoding="utf-8")
    print("added PRIMARY_N import to analyse_families.py")

patch("diag_convergence.py", [(DIAG_OLD, DIAG_NEW)], add_import=True)
patch("diag_ladder.py", [(DIAG_OLD, DIAG_NEW)], add_import=True)

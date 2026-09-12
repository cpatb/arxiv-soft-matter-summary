"""check_prepare_refactor.py — prove the 2026-09-12 `_prepare` extraction moved code.

Stages A-C (insert, compress, relax) were lifted verbatim out of
`shear_at_constant_volume` into `_prepare` so that the new constant-PRESSURE driver runs
the identical preparation. "Verbatim" is a claim, and the thing that settles it is
running the OLD code and the NEW code on the same cell and requiring bit-identical
output.

The old code comes out of git (`--ref` defaults to the commit before the refactor), is
written beside this file under a different module name, and is deleted afterwards. That
is a stronger check than comparing against a committed RECORD, because a record was
produced on a different day by a different process and only its own campaign's
parameters are reproducible; here both arms run in one interpreter, back to back.

⚠ Read on a CONTINUOUS channel. `phi` is driven onto a target to 1e-5 and `z_ss` is
quantized by the contact count, so both can agree while the trajectory differs
(2026-08-31 lesson). `p_ss_hat`, `tau_ss_hat` and `mu_macro` cannot hide a difference.

    python check_prepare_refactor.py          # ~2 min, needs Yade
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_PATH = "Desktop/Projects/polydisperse_jamming/shear.py"
TMP_MODULE = "shear_pre_refactor"

# One cheap cell. Small and short because the extraction can only affect stages A-C and
# those run once per cell regardless of how long the shear is.
CELL = {"cv": 0.0, "mu_p": 0.3, "phi_target": 0.605, "n": 600, "seed": 1,
        "gamma_max": 1.5, "shear_steps": 40_000}
CONTINUOUS = ("p_ss_hat", "tau_ss_hat", "mu_macro", "steady_split", "i_ss")
QUANTIZED = ("phi", "z_ss", "z_nr_ss", "rattler_frac_ss", "gdot", "dt", "phi_drift")


def fetch_old(ref: str) -> str:
    """Write the pre-refactor shear.py beside this file, under another module name.

    Prefers the SIDECAR `shear_pre_refactor_src.txt`, because the copy this runs from in
    WSL (~/pj) is not a git checkout -- it is a plain file copy of the Windows tree, so
    `git show` is unavailable exactly where Yade is. Generate the sidecar on the Windows
    side with `make_prerefactor_sidecar.sh` and sync it.
    """
    sidecar = os.path.join(HERE, "shear_pre_refactor_src.txt")
    if os.path.exists(sidecar):
        src = open(sidecar, encoding="utf-8").read()
        print(f"using the sidecar {sidecar} ({len(src)} bytes)")
    else:
        top = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=HERE,
                             capture_output=True, text=True, check=True).stdout.strip()
        src = subprocess.run(["git", "show", f"{ref}:{REPO_PATH}"], cwd=top,
                             capture_output=True, text=True, check=True).stdout
    if "def _prepare(" in src:
        raise SystemExit(f"{ref} already contains _prepare -- pass an earlier --ref")
    if "def shear_at_constant_volume(" not in src:
        raise SystemExit(f"{ref}:{REPO_PATH} does not look like shear.py")
    path = os.path.join(HERE, f"{TMP_MODULE}.py")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(src)
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD",
                    help="git ref holding the PRE-refactor shear.py")
    args = ap.parse_args()

    path = fetch_old(args.ref)
    print(f"pre-refactor shear.py from {args.ref} -> {path}")
    try:
        old = importlib.import_module(TMP_MODULE)
        import shear as new
        assert not hasattr(old, "_prepare"), "the old module must not have _prepare"
        assert hasattr(new, "_prepare"), "the new module must have _prepare"
        print(f"cell: {CELL}")
        a = json.loads(old.shear_at_constant_volume(**CELL, verbose=True).to_json())
        b = json.loads(new.shear_at_constant_volume(**CELL, verbose=True).to_json())
    finally:
        if os.path.exists(path):
            os.remove(path)
        pyc = os.path.join(HERE, "__pycache__", f"{TMP_MODULE}.cpython-312.pyc")
        if os.path.exists(pyc):
            os.remove(pyc)

    bad = []
    print(f"\n{'field':18s}{'pre-refactor':>24s}{'post':>24s}  identical")
    for group, fields in (("CONTINUOUS (the gate)", CONTINUOUS),
                          ("quantized (reported)", QUANTIZED)):
        print(f"-- {group}")
        for f in fields:
            x, y = a.get(f), b.get(f)
            same = (x == y)
            print(f"{f:18s}{x!r:>24}{y!r:>24}  {'yes' if same else 'NO'}")
            if not same and group.startswith("CONTINUOUS"):
                bad.append((f, x, y))

    if bad:
        print("\n🔴 THE MOVE WAS NOT VERBATIM:")
        for f, x, y in bad:
            rel = abs(y - x) / abs(x) if x else float("inf")
            print(f"   {f}: {x!r} -> {y!r}  (rel {rel:.3e})")
        print("   Nothing measured through either driver is comparable with the "
              "records on disk until this is explained.")
        return 1
    print("\n✅ BIT-IDENTICAL on every continuous channel: the extraction moved code "
          "and changed no trajectory.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

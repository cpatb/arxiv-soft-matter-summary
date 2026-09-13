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


def fetch_old(ref: str, sidecar_name: str, old_has_prepare: bool) -> str:
    """Write the pre-change shear.py beside this file, under another module name.

    Prefers the SIDECAR (default `shear_pre_refactor_src.txt`), because the copy this
    runs from in WSL (~/pj) is not a git checkout -- it is a plain file copy of the
    Windows tree, so `git show` is unavailable exactly where Yade is. Generate the
    sidecar on the Windows side (`git show <ref>:<path> > sidecar`) and sync it.

    `old_has_prepare` selects which change is being checked: False for the 2026-09-12
    `_prepare` extraction (the old source must NOT define `_prepare`), True for any
    later change to the default path (the old source is HEAD and does).
    """
    sidecar = os.path.join(HERE, sidecar_name)
    if os.path.exists(sidecar):
        src = open(sidecar, encoding="utf-8").read()
        print(f"using the sidecar {sidecar} ({len(src)} bytes)")
    else:
        top = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=HERE,
                             capture_output=True, text=True, check=True).stdout.strip()
        src = subprocess.run(["git", "show", f"{ref}:{REPO_PATH}"], cwd=top,
                             capture_output=True, text=True, check=True).stdout
    if ("def _prepare(" in src) != old_has_prepare:
        raise SystemExit(f"{ref}: `_prepare` present={'def _prepare(' in src}, "
                         f"expected {old_has_prepare} -- wrong --ref/--sidecar")
    if "def shear_at_constant_volume(" not in src:
        raise SystemExit(f"{ref}:{REPO_PATH} does not look like shear.py")
    path = os.path.join(HERE, f"{TMP_MODULE}.py")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(src)
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="HEAD",
                    help="git ref holding the PRE-change shear.py")
    ap.add_argument("--sidecar", default="shear_pre_refactor_src.txt",
                    help="file beside this script holding the PRE-change shear.py")
    ap.add_argument("--old-has-prepare", action="store_true",
                    help="the old source already defines _prepare (checking a change "
                         "LATER than the extraction, e.g. the 2026-09-12 mu_prep knob)")
    ap.add_argument("--fprep-smoke", action="store_true",
                    help="also run the NEW code with mu_prep=0 on the same cell and "
                         "report the friction-switch read-back and p_hat_prep")
    args = ap.parse_args()

    path = fetch_old(args.ref, args.sidecar, args.old_has_prepare)
    print(f"pre-change shear.py from {args.ref} -> {path}")
    c = None
    try:
        old = importlib.import_module(TMP_MODULE)
        import shear as new
        assert hasattr(old, "_prepare") == args.old_has_prepare, "old-module _prepare"
        assert hasattr(new, "_prepare"), "the new module must have _prepare"
        print(f"cell: {CELL}")
        a = json.loads(old.shear_at_constant_volume(**CELL, verbose=True).to_json())
        b = json.loads(new.shear_at_constant_volume(**CELL, verbose=True).to_json())
        if args.fprep_smoke:
            # The knob's NON-default setting, on the channel the campaign reads
            # (2026-08-30 lesson: a byte-identity gate on the default alone proves
            # nothing about the knob).
            c = json.loads(new.shear_at_constant_volume(
                **CELL, mu_prep=0.0, n_samples=600, verbose=True).to_json())
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
    print("\n✅ BIT-IDENTICAL on every continuous channel: the change touched no "
          "trajectory on the default path.")
    if c is not None:
        print("\n-- mu_prep = 0 smoke arm (the knob's non-default setting) --")
        for f in ("protocol", "mu_prep", "friction_switch_ok", "p_hat_prep",
                  "z_nr_prep", "rattler_frac_prep", "n_samples", "p_ss_hat",
                  "tau_ss_hat", "mu_macro", "z_nr_ss", "steady", "i_ss"):
            print(f"{f:20s}{c.get(f)!r}")
        print(f"samples recorded: {len(c['gamma'])} (z_nr_series {len(c['z_nr_series'])})")
        for note in c["notes"]:
            print("   note:", note)
        if not c["friction_switch_ok"]:
            print("🔴 friction switch read-back FAILED")
            return 1
        if c["protocol"] != "cv-shear-v3-fprep" or b["protocol"] != "cv-shear-v2":
            print("🔴 protocol strings wrong:", b["protocol"], c["protocol"])
            return 1
        # The arm must be LIVE: same cell, different preparation, different p_ss --
        # a knob that produces the default's number is inert, not verified.
        if c["p_ss_hat"] == b["p_ss_hat"]:
            print("🔴 the mu_prep arm returned the default's p_ss_hat exactly: inert")
            return 1
        print(f"p_ss_hat default {b['p_ss_hat']:.4e} vs frictionless-prep "
              f"{c['p_ss_hat']:.4e} (ratio {c['p_ss_hat'] / b['p_ss_hat']:.3f}); "
              f"p_hat_prep {c['p_hat_prep']:.3e} vs default's stage-C p (not recorded "
              f"on v2). Whether the critical state forgets the preparation is the "
              f"campaign's G3, not this smoke test's.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

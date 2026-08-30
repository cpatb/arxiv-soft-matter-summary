"""Clone a solved case's SETUP at a different reference wind speed.

WHY THIS EXISTS
---------------
A theory of throwability as a function of wind speed needs the CFD field at
every wind speed. Re-running OpenFOAM per speed is unaffordable. The usual
escape is Reynolds-number independence: for a bluff-body flow at Re >> 10^4 the
steady RANS solution is claimed to be self-similar, so

    U(x; Uref) = (Uref/Uref0) * U(x; Uref0)
    k(x; Uref) = (Uref/Uref0)^2 * k(x; Uref0)
    omega(x; Uref) = (Uref/Uref0) * omega(x; Uref0)

That is an ASSUMPTION, not an identity, for this case: molecular nu = 1.5e-5 is
finite and the `buildings` patch uses the SMOOTH-wall nutkWallFunction, whose
y+ dependence carries an explicit Reynolds number. The ground's
atmNutkWallFunction is roughness-based and is scale-free; the buildings' is not.

So this script builds the cases that TEST the assumption instead of asserting
it. Everything except the three inlet `Uref` entries and the initial-field
magnitudes is byte-identical to the source case.

Usage:  python3 makeScaleCase.py cases/dir_090_gridL1 10.0
"""
import json
import re
import shutil
import sys
from pathlib import Path

SRC_UREF = 5.0          # every template in cases/ is built at Uref = 5 m/s


def main():
    src = Path(sys.argv[1])
    uref = float(sys.argv[2])
    s = uref / SRC_UREF
    tag = f"u{uref:g}".replace(".", "p")
    dst = src.parent / f"{src.name}_{tag}"
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)

    # Copy the SETUP only. Time directories, processor*/ and logs are results of
    # the source run and must not be inherited -- a stale 0/ would silently
    # start the new case from the old speed's field.
    for sub in ("0.orig", "constant", "system"):
        shutil.copytree(src / sub, dst / sub)
    for extra in ("EMPTY_FETCH",):
        if (src / extra).exists():
            shutil.copy2(src / extra, dst / extra)

    # case_info.json is REQUIRED, not optional: sample_disc_altitude.py reads
    # theta_wind_deg and geometry_rotation_deg from it to rotate the sampled
    # velocities back into site coordinates. The first run of this script
    # copied a file called `case_meta.json`, which does not exist, so the
    # solve and the 800-iteration average both completed and the sampling
    # stage died -- 15 minutes of compute with nothing readable at the end.
    info = json.loads((src / "case_info.json").read_text())
    info["u_ref"] = uref
    for key, power in (("u_star", 1), ("k_inlet_ref", 2), ("omega_inlet_ref", 1)):
        info[key] = info[key] * s ** power
    info["scale_test"] = dict(source_case=src.name, source_u_ref=SRC_UREF,
                              factor=s)
    (dst / "case_info.json").write_text(json.dumps(info, indent=2))
    shutil.rmtree(dst / "constant" / "polyMesh", ignore_errors=True)

    # --- the three inlet Uref entries ---------------------------------------
    n_uref = 0
    for f in ("U", "k", "omega"):
        p = dst / "0.orig" / f
        t = p.read_text()
        t, n = re.subn(r"(\n\s*Uref\s+)5\s*;", rf"\g<1>{uref:g};", t)
        n_uref += n
        p.write_text(t)
    if n_uref != 3:
        raise SystemExit(f"expected 3 Uref entries, patched {n_uref}")

    # --- initial fields, scaled so the solve starts in the right ballpark ----
    pU = dst / "0.orig" / "U"
    t = pU.read_text()
    t, n = re.subn(r"internalField\s+uniform\s+\(5 0 0\)\s*;",
                   f"internalField   uniform ({uref:g} 0 0);", t)
    if n != 1:
        raise SystemExit("U internalField not patched")
    pU.write_text(t)

    for f, power in (("k", 2), ("omega", 1)):
        p = dst / "0.orig" / f
        t = p.read_text()
        m = re.search(r"internalField\s+uniform\s+([0-9.eE+-]+)\s*;", t)
        if not m:
            raise SystemExit(f"{f} internalField not found")
        val = float(m.group(1)) * s ** power
        t = t[:m.start()] + f"internalField   uniform {val:.6g};" + t[m.end():]
        p.write_text(t)

    (dst / "SCALE_TEST").write_text(
        f"cloned from {src.name} at Uref = {uref:g} m/s (s = {s:g} x {SRC_UREF:g})\n"
        "Only the inlet Uref entries and the initial-field magnitudes differ.\n")
    print(f"wrote {dst}  (Uref = {uref:g}, s = {s:g})")


if __name__ == "__main__":
    main()

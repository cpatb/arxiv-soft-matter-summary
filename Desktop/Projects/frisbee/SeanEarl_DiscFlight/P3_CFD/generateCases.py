"""
Clone caseTemplate/ into cases/dir_NNN/ for each inflow heading.

For each direction, patch:
  - inlet/outlet face designation in blockMeshDict (rotate boundary type)
  - U inlet flowDir vector
  - z0 in U, k, omega, nut (0.001 over lake sector, 0.5 over urban)

Heading convention: 0° = +x (east-blowing wind, i.e. wind FROM east),
                   90° = +y (north-blowing), 180° = -x (west-blowing).
For Sean Earl, 0° is wind FROM lake (easterly). The 180° sweep covers
0..180.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import re
import math


LAKE_SECTOR = (30, 150)    # degrees; outside this, treat as urban roughness


def z0_for_heading(heading_deg: float) -> float:
    a, b = LAKE_SECTOR
    if a <= heading_deg <= b:
        return 0.001
    # transitional band 15 deg
    if (a - 15) <= heading_deg < a or b < heading_deg <= (b + 15):
        return 0.05
    return 0.5


def flow_direction_vector(heading_deg: float) -> tuple[float, float, float]:
    rad = math.radians(heading_deg)
    return (math.cos(rad), math.sin(rad), 0.0)


def patch_dict_files(case_dir: Path, heading_deg: float, u_ref: float):
    z0 = z0_for_heading(heading_deg)
    fx, fy, fz = flow_direction_vector(heading_deg)
    Ux, Uy, Uz = u_ref * fx, u_ref * fy, u_ref * fz

    replacements = [
        ("0/U",     r"flowDir\s+\(\d+\s+\d+\s+\d+\);",
                    f"flowDir ({fx:.4f} {fy:.4f} {fz:.4f});"),
        ("0/U",     r"value\s+uniform\s+\(\d+\s+\d+\s+\d+\);",
                    f"value uniform ({Ux:.4f} {Uy:.4f} {Uz:.4f});"),
        ("0/U",     r"z0\s+uniform\s+[0-9.]+;",          f"z0 uniform {z0};"),
        ("0/U",     r"Uref\s+[0-9.]+;",                  f"Uref {u_ref:.4f};"),
        ("0/k",     r"flowDir\s+\(\d+\s+\d+\s+\d+\);",
                    f"flowDir ({fx:.4f} {fy:.4f} {fz:.4f});"),
        ("0/k",     r"z0\s+uniform\s+[0-9.]+;",          f"z0 uniform {z0};"),
        ("0/k",     r"Uref\s+[0-9.]+;",                  f"Uref {u_ref:.4f};"),
        ("0/omega", r"flowDir\s+\(\d+\s+\d+\s+\d+\);",
                    f"flowDir ({fx:.4f} {fy:.4f} {fz:.4f});"),
        ("0/omega", r"z0\s+uniform\s+[0-9.]+;",          f"z0 uniform {z0};"),
        ("0/omega", r"Uref\s+[0-9.]+;",                  f"Uref {u_ref:.4f};"),
        ("0/nut",   r"z0\s+uniform\s+[0-9.]+;",          f"z0 uniform {z0};"),
    ]
    for rel, pattern, repl in replacements:
        path = case_dir / rel
        if not path.exists():
            continue
        text = path.read_text()
        text = re.sub(pattern, repl, text)
        path.write_text(text)


def clone_case(template: Path, dest: Path, heading_deg: float, u_ref: float):
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(template, dest)
    patch_dict_files(dest, heading_deg, u_ref)
    Path(dest / "Allrun").chmod(0o755)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", default="caseTemplate")
    ap.add_argument("--out", default="cases")
    ap.add_argument("--headings", nargs="+", type=float, required=True)
    ap.add_argument("--u-ref", type=float, default=5.0)
    ap.add_argument("--z-ref", type=float, default=10.0)
    args = ap.parse_args()

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    template = Path(args.template)

    for h in args.headings:
        dest = out_root / f"dir_{int(round(h)):03d}"
        print(f"clone {template} -> {dest}  (heading {h}°, z0={z0_for_heading(h)})")
        clone_case(template, dest, h, args.u_ref)


if __name__ == "__main__":
    main()

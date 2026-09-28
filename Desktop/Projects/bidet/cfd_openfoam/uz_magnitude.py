"""uz_magnitude.py <case_dir> [time] -- max |Ux|, |Uy|, |Uz| in a written U field.

The wedge is one cell thick, so Uz (out of plane) should be ~0. A stalled Uz residual is
either a real out-of-plane mode (Uz comparable to the in-plane velocity somewhere) or a
normalisation artefact of a near-zero field. This reads the field to tell them apart.
"""
import os
import re
import sys

import numpy as np

case = sys.argv[1]
times = sorted((d for d in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", d) and d != "0"),
               key=float)
t = sys.argv[2] if len(sys.argv) > 2 else times[-1]
txt = open(os.path.join(case, t, "U")).read()
m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)\s*\n?\(", txt)
n = int(m.group(1))
body = txt[m.end():]
vals = np.array(re.findall(r"\(([^()]+)\)", body[: body.find("\n)\n") + 3])[:n],
                dtype=object)
U = np.array([[float(x) for x in v.split()] for v in vals])
assert U.shape == (n, 3), U.shape
mag = np.linalg.norm(U[:, :2], axis=1)
print(f"{os.path.basename(case)} t={t} N={n}  max|Ux| {abs(U[:,0]).max():.4g}  "
      f"max|Uy| {abs(U[:,1]).max():.4g}  max|Uz| {abs(U[:,2]).max():.4g}  "
      f"max|Uz|/max|U_inplane| {abs(U[:,2]).max()/mag.max():.3e}  "
      f"rms Uz {np.sqrt((U[:,2]**2).mean()):.3e}")
# Uz on a wedge at angle a is expected ~ Uy*tan(a) in the cell-centre frame? report ratio
i = np.argmax(abs(U[:, 2]))
print(f"   worst cell {i}: U = {U[i]}")

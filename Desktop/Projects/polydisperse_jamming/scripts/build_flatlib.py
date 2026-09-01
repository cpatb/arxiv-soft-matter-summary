"""Build a flat directory of the shared libraries the SYSTEM CANNOT ALREADY RESOLVE.

Run inside WSL.

Why not just flatten everything: shadowing a system library with an extracted copy of
the nominally same version broke the Qt5 link (`libQt5OpenGL.so.5: undefined symbol
QOpenGLFunctions_1_1::versionProfile`), because the extracted Qt5Gui came from a
differently-configured build than the one the system's Qt5Core/Qt5Gui pair expects.
The system's own Qt5 5.15.18 set is self-consistent; ours is not guaranteed to be.

So: a library goes into flatlib ONLY if `ldconfig -p` cannot find it. That keeps the
system's coherent sets intact and fills exactly the gaps (openblas, gl2ps, glut,
QGLViewer, libQt5OpenGL, the boost 1.90 set, ...).

Symlinks point at RESOLVED real files, so a dangling or self-referential link in the
extracted tree cannot propagate.
"""
import os
import subprocess

P = os.path.expanduser("~/opt/yade")
FLAT = os.path.join(P, "flatlib")
os.makedirs(FLAT, exist_ok=True)


def system_sonames() -> set[str]:
    """Every soname the system loader can already resolve, per `ldconfig -p`."""
    out = subprocess.run(["ldconfig", "-p"], capture_output=True, text=True).stdout
    names = set()
    for line in out.splitlines()[1:]:
        line = line.strip()
        if not line or "=>" not in line:
            continue
        names.add(line.split()[0])
    return names


SYS = system_sonames()
print(f"system can resolve {len(SYS)} sonames")

seen: dict[str, str] = {}
skipped_system = 0
dangling = 0

for root, _dirs, files in os.walk(os.path.join(P, "usr", "lib")):
    if "flatlib" in root:
        continue
    try:
        entries = set(files) | set(os.listdir(root))
    except OSError:
        continue
    for fn in sorted(entries):
        if ".so" not in fn:
            continue
        if fn in SYS:
            skipped_system += 1
            continue
        src = os.path.join(root, fn)
        real = os.path.realpath(src)
        if not os.path.isfile(real):
            dangling += 1
            continue
        if fn in seen:
            continue
        seen[fn] = real
        dst = os.path.join(FLAT, fn)
        if not os.path.lexists(dst):
            os.symlink(real, dst)

print(f"flatlib entries: {len(seen)}  "
      f"(skipped {skipped_system} already-resolvable, {dangling} dangling)")
for k in ("libQt5Core.so.5", "libQt5Gui.so.5", "libQt5OpenGL.so.5",
          "libopenblas.so.0", "libgl2ps.so.1.4", "libglut.so.3",
          "libQGLViewer-qt5.so.2", "libQt5Xml.so.5", "libOpenGL.so.0",
          "libGLU.so.1", "libboost_python314.so.1.90.0",
          "libboost_graph.so.1.90.0"):
    where = seen.get(k)
    tag = where if where else ("(system)" if k in SYS else "MISSING")
    print(f"  {k:35s} {tag}")

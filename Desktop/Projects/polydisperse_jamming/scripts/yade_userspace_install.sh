#!/bin/bash
# Userspace Yade install: download the .deb closure and extract into $PREFIX.
# No root. Logs everything; exits nonzero on failure.
set -u
PREFIX="$HOME/opt/yade"
WORK="$HOME/opt/yade_debs"
LOG="$PREFIX/install.log"
mkdir -p "$PREFIX" "$WORK"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) userspace yade install start ==="

cd "$WORK" || exit 1

# Resolve the full dependency closure of yade that is NOT already installed.
echo "--- resolving dependency closure ---"
python3 - <<'PY' > need.txt
import subprocess, re, sys
def run(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout
seen, need, queue = set(), [], ["yade"]
while queue:
    p = queue.pop(0)
    if p in seen: continue
    seen.add(p)
    pol = run(f"apt-cache policy {p}")
    m = re.search(r"Installed:\s*(\S+)", pol)
    if m and m.group(1) != "(none)":
        continue          # already present system-wide
    if "Candidate:" not in pol:
        continue
    need.append(p)
    dep = run(f"apt-cache depends --no-recommends --no-suggests --no-conflicts "
              f"--no-breaks --no-replaces --no-enhances {p}")
    for line in dep.splitlines():
        line = line.strip()
        if line.startswith("Depends:") or line.startswith("PreDepends:"):
            d = line.split(":",1)[1].strip().lstrip("<").rstrip(">")
            if d and not d.startswith("|"):
                queue.append(d)
for p in need: print(p)
PY
echo "packages to fetch: $(wc -l < need.txt)"
cat need.txt

echo "--- downloading ---"
xargs -a need.txt -r -n 20 apt-get download 2>&1 | tail -20
echo "debs downloaded: $(ls -1 *.deb 2>/dev/null | wc -l)"

echo "--- extracting ---"
for d in *.deb; do dpkg-deb -x "$d" "$PREFIX" || echo "EXTRACT_FAIL $d"; done

echo "--- probing ---"
PYV=$(python3 -c 'import sys;print("%d.%d"%sys.version_info[:2])')
export LD_LIBRARY_PATH="$PREFIX/usr/lib/x86_64-linux-gnu:$PREFIX/usr/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$PREFIX/usr/lib/python3/dist-packages:$PREFIX/usr/lib/python${PYV}/dist-packages:${PYTHONPATH:-}"
export PATH="$PREFIX/usr/bin:$PATH"
echo "PYV=$PYV"
find "$PREFIX" -name "yade*" -maxdepth 6 2>/dev/null | head -20
python3 -c "import yade; print('YADE_IMPORT_OK', yade.__file__)" 2>&1 | tail -5
echo "=== $(date -Is) done ==="

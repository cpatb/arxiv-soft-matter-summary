#!/bin/bash
# Userspace Yade install, v2: per-package download (batch failures no longer poison
# the whole set), virtual packages resolved to a concrete provider.
set -u
PREFIX="$HOME/opt/yade"; WORK="$HOME/opt/yade_debs"; LOG="$PREFIX/install2.log"
mkdir -p "$PREFIX" "$WORK"; exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) v2 start ==="
cd "$WORK" || exit 1

python3 - <<'PY' > need2.txt
import subprocess, re
def run(c): return subprocess.run(c, shell=True, capture_output=True, text=True).stdout
def concrete(p):
    """Map a virtual package to a real provider, or return p if already real."""
    pol = run(f"apt-cache policy {p}")
    if "Candidate:" in pol and "Candidate: (none)" not in pol:
        return p
    sp = run(f"apt-cache showpkg {p}")
    if "Reverse Provides:" in sp:
        blk = sp.split("Reverse Provides:")[1].strip().splitlines()
        for line in blk:
            cand = line.split()[0]
            if "Candidate: (none)" not in run(f"apt-cache policy {cand}"):
                return cand
    return None
seen, need, queue = set(), [], ["yade"]
while queue:
    p = queue.pop(0)
    if p in seen: continue
    seen.add(p)
    pol = run(f"apt-cache policy {p}")
    m = re.search(r"Installed:\s*(\S+)", pol)
    if m and m.group(1) != "(none)":
        continue
    c = concrete(p)
    if c is None: 
        print(f"#UNRESOLVED {p}"); continue
    if c != p and c in seen: continue
    need.append(c); seen.add(c)
    dep = run(f"apt-cache depends --no-recommends --no-suggests --no-conflicts "
              f"--no-breaks --no-replaces --no-enhances {c}")
    for line in dep.splitlines():
        line = line.strip()
        if line.startswith(("Depends:", "PreDepends:")):
            d = line.split(":",1)[1].strip().lstrip("<").rstrip(">")
            if d and not d.startswith("|"): queue.append(d)
for p in need: print(p)
PY
grep -v '^#' need2.txt > need2.clean; echo "to fetch: $(wc -l < need2.clean)"
grep '^#' need2.txt || true

echo "--- downloading one at a time ---"
ok=0; fail=0
while read -r p; do
  [ -z "$p" ] && continue
  if apt-get download "$p" >/dev/null 2>&1; then ok=$((ok+1)); else echo "DL_FAIL $p"; fail=$((fail+1)); fi
done < need2.clean
echo "downloaded ok=$ok fail=$fail ; debs on disk: $(ls -1 *.deb 2>/dev/null | wc -l)"
ls -1 *.deb | grep -i '^yade' || echo "!! NO YADE DEB !!"

echo "--- extracting ---"
for d in *.deb; do dpkg-deb -x "$d" "$PREFIX" 2>/dev/null || echo "EXTRACT_FAIL $d"; done

echo "--- probing ---"
PYV=$(python3 -c 'import sys;print("%d.%d"%sys.version_info[:2])')
export LD_LIBRARY_PATH="$PREFIX/usr/lib/x86_64-linux-gnu:$PREFIX/usr/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$PREFIX/usr/lib/python3/dist-packages:$PREFIX/usr/lib/python${PYV}/dist-packages:${PYTHONPATH:-}"
export PATH="$PREFIX/usr/bin:$PATH"
echo "PYV=$PYV"
echo "yade python dirs:"; find "$PREFIX" -type d -name "yade*" 2>/dev/null | head
echo "yade binaries:";    find "$PREFIX" -name "yade*" -type f 2>/dev/null | head
python3 -c "import yade; print('YADE_IMPORT_OK', yade.__version__ if hasattr(yade,'__version__') else '', yade.__file__)" 2>&1 | tail -8
echo "=== $(date -Is) v2 done ==="

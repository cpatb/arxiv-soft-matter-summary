#!/bin/bash
# Iteratively resolve missing shared libraries for the userspace Yade extraction.
# Each round: try `import yade`, read the missing .so from the ImportError, find the
# providing package with apt-file/dpkg -S fallback to apt-cache search, download, extract.
set -u
P="$HOME/opt/yade"; W="$HOME/opt/yade_debs"; LOG="$P/fixlibs.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== $(date -Is) fix-libs start ==="
cd "$W" || exit 1
export LD_LIBRARY_PATH="$P/usr/lib/x86_64-linux-gnu/yade-double:$P/usr/lib/x86_64-linux-gnu:$P/usr/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$P/usr/lib/x86_64-linux-gnu/yade-double/py:${PYTHONPATH:-}"

# Candidate providers for the libs Yade is known to need that are virtual in apt.
declare -A HINT=(
  [libopenblas.so.0]=libopenblas0-pthread
  [liblapack.so.3]=liblapack3
  [libblas.so.3]=libblas3
)

for round in $(seq 1 25); do
  err=$(python3 -c "import yade" 2>&1)
  if [ -z "$err" ] || echo "$err" | grep -q "^Using python version"; then
    if python3 -c "import yade" 2>/dev/null; then echo "ROUND $round: IMPORT OK"; break; fi
  fi
  miss=$(echo "$err" | grep -oP '(?<=ImportError: )[^:]+\.so[^:]*(?=: cannot open)' | head -1)
  if [ -z "$miss" ]; then
    echo "ROUND $round: no missing-.so pattern; remaining error:"; echo "$err" | tail -15; break
  fi
  echo "ROUND $round: missing $miss"
  pkg="${HINT[$miss]:-}"
  if [ -z "$pkg" ]; then
    base=$(echo "$miss" | sed 's/\.so.*//')
    soname=$(echo "$miss" | sed 's/^lib//; s/\.so\./)/' )
    # guess: libfoo.so.5  ->  libfoo5 ; libfoo.so.5.1 -> libfoo5
    guess=$(echo "$miss" | sed -E 's/^(lib.*)\.so\.([0-9]+).*$/\1\2/')
    guess2=$(echo "$miss" | sed -E 's/^(lib.*)\.so\.([0-9]+).*$/\1\2t64/')
    for g in "$guess" "$guess2" "${base}"; do
      if apt-cache policy "$g" 2>/dev/null | grep -q "Candidate: [^(]"; then pkg="$g"; break; fi
    done
  fi
  if [ -z "$pkg" ]; then echo "  !! no provider guess for $miss -- STOP"; break; fi
  echo "  provider: $pkg"
  if apt-get download "$pkg" >/dev/null 2>&1; then
    newdeb=$(ls -t ${pkg}_*.deb 2>/dev/null | head -1)
    echo "  extracting $newdeb"
    dpkg-deb -x "$newdeb" "$P" || echo "  EXTRACT_FAIL"
    # also pull that package's own deps
    for d in $(apt-cache depends --no-recommends --no-suggests --no-conflicts --no-breaks \
               --no-replaces --no-enhances "$pkg" | grep -oP '(?<=Depends: ).*'); do
      if ! apt-cache policy "$d" 2>/dev/null | grep -q "Installed: [^(]"; then
        apt-get download "$d" >/dev/null 2>&1 && dpkg-deb -x $(ls -t ${d}_*.deb 2>/dev/null|head -1) "$P" 2>/dev/null
      fi
    done
  else
    echo "  DL_FAIL $pkg -- STOP"; break
  fi
done

echo "--- final probe ---"
python3 -c "
import yade
from yade import pack, utils
print('YADE_IMPORT_OK')
print('version:', yade.libVersions.getVersion('yade') if hasattr(yade,'libVersions') else 'n/a')
" 2>&1 | tail -12
echo "=== $(date -Is) fix-libs done ==="

#!/bin/bash
# Environment for the userspace Yade install. Source this before any Yade run.
#
#   source ~/pj/scripts/env.sh
#
# Why this exists: `sudo` on this WSL needs a password, so Yade 2026.1.0-3build1 was
# installed by extracting its .deb dependency closure (101 packages) into ~/opt/yade
# with `dpkg-deb -x`. Three things had to be got right and each one silently broke the
# import when it was wrong:
#
#   1. PYTHONPATH points at   .../yade-double/py   -- NOT dist-packages.
#   2. flatlib holds ONLY the sonames `ldconfig -p` cannot already resolve. Shadowing
#      the system's self-consistent Qt5 5.15.18 set with our extracted copies produced
#      `libQt5OpenGL.so.5: undefined symbol QOpenGLFunctions_1_1::versionProfile`.
#      Rebuild it with scripts/build_flatlib.py, never by flattening everything.
#   3. IPython is the only Python dep the system lacks; numpy/scipy/matplotlib/pandas
#      are already system-wide. Installing the whole stack into ~/opt/pylibs pulled a
#      numpy 2.5.2 that conflicted with the system numba. Install IPython ALONE.

export YADE_PREFIX="$HOME/opt/yade"
export LD_LIBRARY_PATH="$YADE_PREFIX/flatlib:$YADE_PREFIX/usr/lib/x86_64-linux-gnu/yade-double:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$YADE_PREFIX/usr/lib/x86_64-linux-gnu/yade-double/py:$HOME/opt/pylibs:${PYTHONPATH:-}"
export PJ_ROOT="/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming"

# Deterministic, reproducible, and quiet.
#
# Yade DOES honour OMP_NUM_THREADS (verified: 1/2/3 -> O.numThreads 1/2/3; it is
# read-only after import, so it must be set in the environment).
#
# OMP_WAIT_POLICY=PASSIVE matters more than the thread count here. With the default
# active policy, idle OpenMP threads SPIN, so 8 workers x 2 threads showed a load
# average of 42 on 16 physical cores and the sweep crawled at 0.43 cells/min. For a
# many-small-cells sweep the right configuration is many SINGLE-threaded workers:
# no OpenMP overhead and perfect parallel efficiency across cells.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_WAIT_POLICY=PASSIVE
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PYTHONHASHSEED=0
export MPLBACKEND=Agg

pj_python() { python3 "$@"; }

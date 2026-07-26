#!/bin/bash
# Run all cases under the supplied parent directory sequentially.
# Use `xargs -P` if you have headroom for parallel runs.
set -euo pipefail

ROOT=${1:?usage: runSweep.sh cases/}

for d in "$ROOT"/dir_*/; do
    echo "==== $d ===="
    pushd "$d" >/dev/null
    ./Allrun 2>&1 | tee log.allrun
    popd >/dev/null
done

echo "All cases complete."

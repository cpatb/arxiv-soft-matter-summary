#!/bin/bash
cd "$(dirname "$0")"
bash averagePass.sh 11 800 cases/dir_090 cases/dir_090_coarse cases/dir_210 cases/dir_270
echo "### verifying averaged fields exist ###"
for c in cases/dir_090 cases/dir_090_coarse cases/dir_210 cases/dir_270; do
    latest=$(ls -d $c/[0-9]* 2>/dev/null | sort -t/ -k2 -g | tail -1)
    if [ -f "$latest/UMean" ]; then echo "OK   $c -> $latest/UMean"; else echo "MISSING UMean in $c (latest=$latest)"; fi
done
echo "### starting downstream ###"
bash ../runDownstream.sh 20 150

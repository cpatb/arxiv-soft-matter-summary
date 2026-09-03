#!/bin/bash
set -u
cd ~/pj || exit 1
source scripts/env.sh
rm -f out/diag_confine.jsonl
stdbuf -oL python3 -u diag_confine.py > out/dc2.log 2>&1
echo "exit=$?" >> out/dc2.log

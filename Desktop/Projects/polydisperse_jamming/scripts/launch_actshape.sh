#!/bin/bash
# actshape -- PREREG_ACTSHAPE_2026-09-22.md. 54 cells, per-contact dumps.
#
# Step counts are per (CV, p) so gdot is matched across CV (the rate-matched ladder) and
# halved at p = 1e-4 so I stays under I_MAX. Workers per group are sized so every group
# finishes in ~80 min (costs from the rate-matched ladder: ~410 s per 102k steps).
#
# ⚠ Foreground: the caller keeps this wsl.exe alive in a background tool call. Nothing
# spawned by a wsl.exe call survives it exiting. Restartable: skips keys already written.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
DUMP=~/pj/out/actshape_dump
mkdir -p "$DUMP"

pids=()
group () {  # cv  p-stars  steps  nworkers  nick
  local cv="$1" ps="$2" steps="$3" nw="$4" nick="$5"
  for i in $(seq 0 $((nw - 1))); do
    LOG=~/pj/out/actshape_${nick}_w${i}.log
    stdbuf -oL python3 -u phase9_constp.py \
        --cvs "$cv" --mus 0.0,0.3 --p-stars "$ps" --phi-starts 0.55 \
        --gains 0.15 --seeds 1,2,3 --n 1000 --tag "actshape_${nick}" \
        --shear-steps "$steps" --shard "$i" --nshards "$nw" \
        --dump-dir "$DUMP" > "$LOG" 2>&1 &
    pids+=($!)
  done
  echo "group $nick: cv=$cv p=$ps steps=$steps workers=$nw"
}

group 0.0 3e-4,1e-3 102000 1 c00a
group 0.0 1e-4      204000 1 c00b
group 0.2 3e-4,1e-3 200000 2 c20a
group 0.2 1e-4      400000 2 c20b
group 0.5 3e-4,1e-3 480000 4 c50a
group 0.5 1e-4      960000 6 c50b

echo "actshape: ${#pids[@]} workers"
date -u +"start %Y-%m-%dT%H:%M:%SZ"
fail=0
for pid in "${pids[@]}"; do wait "$pid" || fail=$((fail+1)); done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
echo "ALL DONE, $fail nonzero exits"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp ~/pj/out/phase9_actshape_*.jsonl "$SRC/out/"
mkdir -p "$SRC/out/actshape_dump"
cp "$DUMP"/*.npz "$SRC/out/actshape_dump/"
cp ~/pj/out/actshape_*_w*.log "$SRC/out/"
echo "synced: $(cat "$SRC"/out/phase9_actshape_*.jsonl | wc -l) records, $(ls "$SRC"/out/actshape_dump | wc -l) dumps"

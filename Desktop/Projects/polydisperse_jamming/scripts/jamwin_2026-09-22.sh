#!/bin/bash
# jamwin_2026-09-22.sh -- PREREG_JAMWIN_2026-09-22.md (polydisperse_jamming).
# Does the activation shape delta decide whether an impacting drop jams, with the sign
# sigma_jam = sigma* [ln(1/f_c)]^(-1/delta) predicts? 27 cells + 2 identity cells, L9.
#
# Run from PowerShell by file path, attached, in a background tool call:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/jamwin_2026-09-22.sh
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/jamwin_2026-09-22
REPO="/mnt/c/Users/compl/Github/Drop Impact/basilisk-drop-impact"
OLD=$W/campaigns/deltastiff_2026-09-03/DS_d1.3169_L9
NPAR=${NPAR:-6}
CELL_TIMEOUT=${CELL_TIMEOUT:-14400}
mkdir -p "$ROOT/build"

# --- build from the repo source (with JAM_FHAT) --------------------------------
tr -d '\r' < "$REPO/wc_drop_impact.c" > "$ROOT/build/wc_drop_impact.c"
( cd "$ROOT/build" && qcc -O2 -Wall -disable-dimensions wc_drop_impact.c -o wcdrop -lm \
    > compile.log 2>&1 ) || { echo "FATAL: build failed"; tail -20 "$ROOT/build/compile.log"; exit 1; }
BIN=$ROOT/build/wcdrop
MD5=$(md5sum "$BIN" | cut -c1-12)
grep -q "JAM_FHAT" "$ROOT/build/wc_drop_impact.c" || { echo "FATAL: source lacks JAM_FHAT"; exit 1; }
echo "=== jamwin build ok bin=$MD5 src=$(md5sum "$ROOT/build/wc_drop_impact.c" | cut -c1-12) ==="

U0=2.0; TH=5; TEND=6; CFL=0.10; L=9
PHI0=0.62; PHIM=0.496

cell () {  # label  phi  delta  sigma*  jam_fhat(or "")  [preset]
  local lab=$1 phi=$2 d=$3 ss=$4 jf=$5 o="$ROOT/$1"
  [ -f "$o/RUN_STATUS.json" ] && { echo "  skip $lab"; return 0; }
  mkdir -p "$o"; cd "$o" || return 1
  local s; s=$(date +%s)
  local extra=()
  [ -n "$jf" ] && extra+=("JAM_FHAT=$jf")
  if [ "$ss" = "DS" ]; then       # identity cells: deltastiff's exact configuration
    timeout "$CELL_TIMEOUT" env LEVEL=$L TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 PHI=0.44 \
        WC_CLOCK=1 THETA0=$TH DELTA=1.3169 "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$U0" cs_fit "$L" > out.log 2> err.log
  else
    timeout "$CELL_TIMEOUT" env LEVEL=$L TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 \
        PHI="$phi" PHI0=$PHI0 PHIM=$PHIM SIGMA_STAR="$ss" DELTA="$d" \
        WC_CLOCK=1 THETA0=$TH "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$U0" cs_fit "$L" > out.log 2> err.log
  fi
  local rc=$? e=$(( $(date +%s) - s ))
  rm -f "$o"/dump-* "$o"/*.mp4
  local dat n=0 lastt=0
  dat=$(ls "$o"/spread_*.dat 2>/dev/null | head -1)
  if [ -n "$dat" ]; then
    n=$(grep -vc '^#' "$dat"); lastt=$(grep -v '^#' "$dat" | tail -1 | awk '{print $2}')
  fi
  local ban; ban=$(grep -a -m1 '^# U0=' err.log | tr -d '"')
  local jb; jb=$(grep -a -o 'JAM_FHAT=[0-9.e+-]*' err.log | head -1 | cut -d= -f2)
  cat > "$o/RUN_STATUS.json" <<EOJ
{"label":"$lab","campaign":"jamwin_2026-09-22","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt jam_fhat=${jb:-unset}"
}

# --- G0 identity first, alone -----------------------------------------------------
echo "=== G0 identity $(date -Is) ==="
( cell G0a_unset x x DS "" ) & ( cell G0b_half x x DS 0.5 ) & wait
for g in G0a_unset G0b_half; do
  a=$(md5sum "$ROOT/$g"/spread_*.dat | cut -c1-32); b=$(md5sum "$OLD"/spread_*.dat | cut -c1-32)
  if [ "$a" = "$b" ]; then echo "G0 $g PASS md5 $a"; else echo "G0 $g FAIL $a != $b"; fi
done

# --- the 27-cell grid -------------------------------------------------------------
LIST=()
for arm in A:0.589:0.25 B:0.527:0.75 N:0.465:0.5; do
  IFS=: read -r nm phi jf <<< "$arm"
  for ss in 30 300 3000; do for d in 0.7 1.3169 2.6; do
    LIST+=("${nm}_s${ss}_d${d} $phi $d $ss $jf")
  done; done
done
echo "=== grid ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== jamwin complete $(date -Is) ==="

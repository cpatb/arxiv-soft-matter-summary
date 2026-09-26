#!/bin/bash
# dsttype -- PREREG_DSTTYPE_2026-09-26.md: 20 cells, EPS_MIN 0.03, L9, the cvmap binary. Cells from out/dsttype_cells.txt.
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/dsttype_2026-09-26
BUILD=$W/campaigns/jamwin_2026-09-22/build
CELLS=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/dsttype_cells.txt
NPAR=${NPAR:-10}
CELL_TIMEOUT=${CELL_TIMEOUT:-14400}
mkdir -p "$ROOT"
BIN=$BUILD/wcdrop
[ -x "$BIN" ] || { echo "FATAL: $BIN missing"; exit 1; }
MD5=$(md5sum "$BIN" | cut -c1-12)
[ "$MD5" = "db413e87dd21" ] || { echo "FATAL: binary md5 $MD5 != the G0-gated db413e87dd21"; exit 1; }
echo "=== dsttype reusing G0-gated binary $MD5 ==="
TH=5; TEND=6; CFL=0.10; L=9; PHI0=0.62; EPS=0.03

cell () {  # label phi delta sigma* jam_fhat phim u0
  local lab=$1 phi=$2 d=$3 ss=$4 jf=$5 pm=$6 u0=$7 o="$ROOT/$1"
  [ -f "$o/RUN_STATUS.json" ] && { echo "  skip $lab"; return 0; }
  mkdir -p "$o"; cd "$o" || return 1
  local s; s=$(date +%s)
  timeout "$CELL_TIMEOUT" env LEVEL=$L TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 \
      PHI="$phi" PHI0=$PHI0 PHIM="$pm" SIGMA_STAR="$ss" DELTA="$d" \
      WC_CLOCK=1 THETA0=$TH JAM_FHAT="$jf" EPS_MIN=$EPS OMP_NUM_THREADS=1 \
      "$BIN" "$u0" cs_fit "$L" > out.log 2> err.log
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
{"label":"$lab","campaign":"dsttype_2026-09-26","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt"
}

mapfile -t LIST < <(tr -d '\r' < "$CELLS")
echo "=== dsttype ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  [ -z "$c" ] && continue
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== dsttype complete $(date -Is) ==="

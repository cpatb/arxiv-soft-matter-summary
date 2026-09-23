#!/bin/bash
# kmermapRG -- PREREG_KMERMAPE03 Amendment 2 (chi from the p_hat 1e-3 rung).
# Run via scripts/run_kmermapRG.sh from PowerShell.
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/kmermapRG_2026-09-23
BUILD=$W/campaigns/jamwin_2026-09-22/build
REPO="/mnt/c/Users/compl/Github/Drop Impact/basilisk-drop-impact"
OLD=$W/campaigns/deltastiff_2026-09-03/DS_d1.3169_L9
NPAR=${NPAR:-8}
CELL_TIMEOUT=${CELL_TIMEOUT:-14400}
mkdir -p "$ROOT"

BIN=$BUILD/wcdrop
[ -x "$BIN" ] || { echo "FATAL: $BIN missing"; exit 1; }
MD5=$(md5sum "$BIN" | cut -c1-12)
[ "$MD5" = "db413e87dd21" ] || { echo "FATAL: binary md5 $MD5 != the G0-gated db413e87dd21"; exit 1; }
echo "=== kmermapRG reusing G0-gated binary $MD5 ==="

U0=2.0; TH=5; TEND=6; CFL=0.10; L=9
PHI0=0.62; PHIM=0.496

cell () {  # label  phi  delta  sigma*  jam_fhat  phim  u0  [level]
  local lab=$1 phi=$2 d=$3 ss=$4 jf=$5 pm=${6:-$PHIM} u0=${7:-$U0} lvl=${8:-$L} eps=${9:-} o="$ROOT/$1"
  [ -f "$o/RUN_STATUS.json" ] && { echo "  skip $lab"; return 0; }
  mkdir -p "$o"; cd "$o" || return 1
  local s; s=$(date +%s)
  local extra=()
  [ -n "$jf" ] && extra+=("JAM_FHAT=$jf")
  [ -n "$eps" ] && extra+=("EPS_MIN=$eps")
  if [ "$ss" = "DS" ]; then       # identity cells: deltastiff's exact configuration
    timeout "$CELL_TIMEOUT" env LEVEL=$lvl TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 PHI=0.44 \
        WC_CLOCK=1 THETA0=$TH DELTA=1.3169 "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$u0" cs_fit "$lvl" > out.log 2> err.log
  else
    timeout "$CELL_TIMEOUT" env LEVEL=$lvl TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 \
        PHI="$phi" PHI0=$PHI0 PHIM="$pm" SIGMA_STAR="$ss" DELTA="$d" \
        WC_CLOCK=1 THETA0=$TH "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$u0" cs_fit "$lvl" > out.log 2> err.log
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
{"label":"$lab","campaign":"kmermapRG_2026-09-23","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt jam_fhat=${jb:-unset}"
}

# --- kmermapRG (PREREG_KMERMAPE03 Amendment 2): the p_hat 1e-3 rung chi, doublets at x 0.91 ---
LIST=()
LIST+=("rg_d_x0.91_U2.0_e03 0.5642 1.1163 235.8 0.8806 0.55664 2.0 9 0.03")
LIST+=("rg_d_x0.91_U4.0_e03 0.5642 1.1163 235.8 0.8806 0.55664 4.0 9 0.03")
LIST+=("rg_d_x0.91_U2.0_e01 0.5642 1.1163 235.8 0.8806 0.55664 2.0 9")
LIST+=("rg_d_x0.91_U4.0_e01 0.5642 1.1163 235.8 0.8806 0.55664 4.0 9")
echo "=== kmermapRG ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== kmermapRG complete $(date -Is) ==="

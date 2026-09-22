#!/bin/bash
# kmerdrop_2026-09-22.sh -- PREREG_KMERDROP_2026-09-22.md. The doublet fraction w enters the
# DNS only through chi(w) (DEM-measured); x = 0.85, sigma* = 300 Pa, delta = 1.3169 fixed.
# 5 cells, L9. Reuses the G0-gated jamwin binary.
# Does the activation shape delta decide whether an impacting drop jams, with the sign
# sigma_jam = sigma* [ln(1/f_c)]^(-1/delta) predicts? 27 cells + 2 identity cells, L9.
#
# Run from PowerShell by file path, attached, in a background tool call:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/jamwin_2026-09-22.sh
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/kmerdrop_2026-09-22
BUILD=$W/campaigns/jamwin_2026-09-22/build
REPO="/mnt/c/Users/compl/Github/Drop Impact/basilisk-drop-impact"
OLD=$W/campaigns/deltastiff_2026-09-03/DS_d1.3169_L9
NPAR=${NPAR:-6}
CELL_TIMEOUT=${CELL_TIMEOUT:-14400}
mkdir -p "$ROOT"

BIN=$BUILD/wcdrop
[ -x "$BIN" ] || { echo "FATAL: $BIN missing"; exit 1; }
MD5=$(md5sum "$BIN" | cut -c1-12)
[ "$MD5" = "db413e87dd21" ] || { echo "FATAL: binary md5 $MD5 != the G0-gated db413e87dd21"; exit 1; }
echo "=== kmerdrop reusing G0-gated binary $MD5 ==="

U0=2.0; TH=5; TEND=6; CFL=0.10; L=9
PHI0=0.62; PHIM=0.496

cell () {  # label  phi  delta  sigma*  jam_fhat  phim
  local lab=$1 phi=$2 d=$3 ss=$4 jf=$5 pm=${6:-$PHIM} o="$ROOT/$1"
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
        PHI="$phi" PHI0=$PHI0 PHIM="$pm" SIGMA_STAR="$ss" DELTA="$d" \
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
{"label":"$lab","campaign":"kmerdrop_2026-09-22","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt jam_fhat=${jb:-unset}"
}

# --- the 21-cell grid -------------------------------------------------------------
LIST=()
# label phi delta sigma* jam_fhat PHIM ; PHI0 = 0.62, x = 0.85 -> PHI = 0.527
for spec in "w0.00:0.54932:1.0" "w0.25:0.53258:1.0" "w0.50:0.51708:0.9036" "w0.75:0.50282:0.7937" "w1.00:0.48856:0.7075"; do
  IFS=: read -r nm pm jf <<< "$spec"
  LIST+=("${nm} 0.527 1.3169 300 $jf $pm")
done
echo "=== kmerdrop ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== kmerdrop complete $(date -Is) ==="

#!/bin/bash
# cvmapE03 -- PREREG_CVMAP_2026-09-24.md Amendment 2: 18 grid cells + W/A at every x0, EPS_MIN 0.03, L9.
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/cvmapE03_2026-09-24
BUILD=$W/campaigns/jamwin_2026-09-22/build
REPO="/mnt/c/Users/compl/Github/Drop Impact/basilisk-drop-impact"
OLD=$W/campaigns/deltastiff_2026-09-03/DS_d1.3169_L9
NPAR=${NPAR:-20}
CELL_TIMEOUT=${CELL_TIMEOUT:-14400}
mkdir -p "$ROOT"

BIN=$BUILD/wcdrop
[ -x "$BIN" ] || { echo "FATAL: $BIN missing"; exit 1; }
MD5=$(md5sum "$BIN" | cut -c1-12)
[ "$MD5" = "db413e87dd21" ] || { echo "FATAL: binary md5 $MD5 != the G0-gated db413e87dd21"; exit 1; }
echo "=== cvmapE03 reusing G0-gated binary $MD5 ==="

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
{"label":"$lab","campaign":"cvmapE03_2026-09-24","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt jam_fhat=${jb:-unset}"
}

# --- cells ---
LIST=()
LIST+=("c0.0_x0.920_U2.0 0.57040 1.3169 300.0 0.8952 0.56459 2.0 9 0.03")
LIST+=("c0.0_x0.920_U4.0 0.57040 1.3169 300.0 0.8952 0.56459 4.0 9 0.03")
LIST+=("c0.2_x0.920_U2.0 0.55950 1.2353 287.1 1.0000 0.56553 2.0 9 0.03")
LIST+=("c0.2_x0.920_U4.0 0.55950 1.2353 287.1 1.0000 0.56553 4.0 9 0.03")
LIST+=("c0.5_x0.920_U2.0 0.51908 1.0733 220.5 1.0000 0.56855 2.0 9 0.03")
LIST+=("c0.5_x0.920_U4.0 0.51908 1.0733 220.5 1.0000 0.56855 4.0 9 0.03")
LIST+=("c0.0_x0.932_U2.0 0.57784 1.3169 300.0 0.7609 0.56459 2.0 9 0.03")
LIST+=("c0.0_x0.932_U4.0 0.57784 1.3169 300.0 0.7609 0.56459 4.0 9 0.03")
LIST+=("c0.2_x0.932_U2.0 0.56679 1.2353 287.1 0.9768 0.56553 2.0 9 0.03")
LIST+=("c0.2_x0.932_U4.0 0.56679 1.2353 287.1 0.9768 0.56553 4.0 9 0.03")
LIST+=("c0.5_x0.932_U2.0 0.52586 1.0733 220.5 1.0000 0.56855 2.0 9 0.03")
LIST+=("c0.5_x0.932_U4.0 0.52586 1.0733 220.5 1.0000 0.56855 4.0 9 0.03")
LIST+=("c0.0_x0.950_U2.0 0.58900 1.3169 300.0 0.5595 0.56459 2.0 9 0.03")
LIST+=("c0.0_x0.950_U4.0 0.58900 1.3169 300.0 0.5595 0.56459 4.0 9 0.03")
LIST+=("c0.2_x0.950_U2.0 0.57774 1.2353 287.1 0.7759 0.56553 2.0 9 0.03")
LIST+=("c0.2_x0.950_U4.0 0.57774 1.2353 287.1 0.7759 0.56553 4.0 9 0.03")
LIST+=("c0.5_x0.950_U2.0 0.53601 1.0733 220.5 1.0000 0.56855 2.0 9 0.03")
LIST+=("c0.5_x0.950_U4.0 0.53601 1.0733 220.5 1.0000 0.56855 4.0 9 0.03")
LIST+=("W_x0.920_U4.0 0.51908 1.3169 300.0 1.0000 0.56855 4.0 9 0.03")
LIST+=("A_x0.920_U4.0 0.57040 1.0733 220.5 0.8952 0.56459 4.0 9 0.03")
LIST+=("W_x0.932_U4.0 0.52586 1.3169 300.0 1.0000 0.56855 4.0 9 0.03")
LIST+=("A_x0.932_U4.0 0.57784 1.0733 220.5 0.7609 0.56459 4.0 9 0.03")
LIST+=("W_x0.950_U4.0 0.53601 1.3169 300.0 1.0000 0.56855 4.0 9 0.03")
LIST+=("A_x0.950_U4.0 0.58900 1.0733 220.5 0.5595 0.56459 4.0 9 0.03")
echo "=== cvmapE03 ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== cvmapE03 complete $(date -Is) ==="

#!/bin/bash
# kmerdecomp -- PREREG_KMERDECOMP_2026-09-23.md. Clone of kmermap_2026-09-23.sh (same binary and settings);
# only ROOT, the campaign name and LIST differ. Run via scripts/run_kmerdecomp.sh from PowerShell.
set -u
source /home/cpb/work/.project_config
W=/home/cpb/work
ROOT=$W/campaigns/kmerdecomp_2026-09-23
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
echo "=== kmerdecomp reusing G0-gated binary $MD5 ==="

U0=2.0; TH=5; TEND=6; CFL=0.10; L=9
PHI0=0.62; PHIM=0.496

cell () {  # label  phi  delta  sigma*  jam_fhat  phim  u0
  local lab=$1 phi=$2 d=$3 ss=$4 jf=$5 pm=${6:-$PHIM} u0=${7:-$U0} o="$ROOT/$1"
  [ -f "$o/RUN_STATUS.json" ] && { echo "  skip $lab"; return 0; }
  mkdir -p "$o"; cd "$o" || return 1
  local s; s=$(date +%s)
  local extra=()
  [ -n "$jf" ] && extra+=("JAM_FHAT=$jf")
  if [ "$ss" = "DS" ]; then       # identity cells: deltastiff's exact configuration
    timeout "$CELL_TIMEOUT" env LEVEL=$L TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 PHI=0.44 \
        WC_CLOCK=1 THETA0=$TH DELTA=1.3169 "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$u0" cs_fit "$L" > out.log 2> err.log
  else
    timeout "$CELL_TIMEOUT" env LEVEL=$L TEND=$TEND CFLNUM=$CFL SHAPE_MODE=1 \
        PHI="$phi" PHI0=$PHI0 PHIM="$pm" SIGMA_STAR="$ss" DELTA="$d" \
        WC_CLOCK=1 THETA0=$TH "${extra[@]}" OMP_NUM_THREADS=1 \
        "$BIN" "$u0" cs_fit "$L" > out.log 2> err.log
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
{"label":"$lab","campaign":"kmerdecomp_2026-09-23","binary_md5":"$MD5","phi":"$phi",
 "delta":"$d","sigma_star":"$ss","jam_fhat_requested":"$jf","jam_fhat_banner":"${jb:-unset}",
 "banner":"$ban","rc":$rc,"wall_s":$e,"n_rows":$n,"last_t_over_tau":${lastt:-0},
 "finished":"$(date -Is)"}
EOJ
  echo "  done $lab rc=$rc ${e}s rows=$n t_end=$lastt jam_fhat=${jb:-unset}"
}

# --- kmerdecomp: 4 channel-isolation cells (PREREG_KMERDECOMP_2026-09-23.md) ---------
LIST=()
LIST+=("win_x0.88_U2.0 0.5456 1.3169 300 0.8766 0.53512 2.0")
LIST+=("win_x0.88_U4.0 0.5456 1.3169 300 0.8766 0.53512 4.0")
LIST+=("win_x0.85_U4.0 0.527 1.3169 300 1.0 0.53512 4.0")
LIST+=("act_x0.85_U4.0 0.527 1.1163 235.8 1.0 0.56315 4.0")
# Amendment 1 (added after the four cells above finished; they are skipped via RUN_STATUS)
LIST+=("act_x0.85_U1.0 0.527 1.1163 235.8 1.0 0.56315 1.0")
LIST+=("act_x0.85_U2.0 0.527 1.1163 235.8 1.0 0.56315 2.0")
echo "=== kmerdecomp ${#LIST[@]} cells, NPAR=$NPAR, start $(date -Is) ==="
running=0
for c in "${LIST[@]}"; do
  ( cell $c ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait -n; running=$((running - 1)); fi
done
wait
echo "=== kmerdecomp complete $(date -Is) ==="

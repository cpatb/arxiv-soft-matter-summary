#!/bin/bash
# shape_rerun_status.sh -- what the raised-cap re-run is doing. Invoke by file path.
cd ~/pj/out || exit 1
echo "=== workers by invocation (family alphas phats nshards) ==="
ps -eo args | grep '[p]hase2_shape.py' | sed 's/.*--family //' | cut -d' ' -f1,3,5,17 | sort | uniq -c
echo "=== rerun log ==="
cat shape_rerun.log 2>/dev/null
echo "=== shpcap records on disk ==="
wc -l phase2_shpcap_*.jsonl 2>/dev/null
echo "=== per-log: cells total / failures / last progress line ==="
for f in shpcap_*_w*.log; do
  [ -e "$f" ] || continue
  printf '%-32s lines=%-5s ' "$f" "$(grep -c . "$f")"
  printf 'N_cells=%s ' "$(grep -o 'N_cells_total=[0-9]*' "$f" | tr '\n' ',')"
  printf 'failed=%s ' "$(grep -c 'CELL FAILED\|Traceback' "$f")"
  grep -h 'iter\|elapsed\|conv=' "$f" | tail -1 | cut -c1-110
  echo
done
echo "=== identity pair: default-cap vs raised-cap phi (once both exist) ==="
python3 - <<'PYEOF'
import glob, json
recs = {}
for f in glob.glob("phase2_shp_lognormal_shard*.jsonl") + glob.glob("phase2_shpcap_lognormal_shard*of02.jsonl"):
    for line in open(f):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("failed") or r["cv"] != 0.0 or r["seed"] != 1 or r["p_hat_target"] != 1e-3:
            continue
        recs.setdefault((r["alpha"], r["max_iter"]), r)
for a in (1.0, 2.0):
    d = recs.get((a, 200000)); c = recs.get((a, 1000000))
    if d and c:
        same = all(d[k] == c[k] for k in ("phi", "iters", "z_all", "p_hat"))
        print(f"  alpha={a}: phi {d['phi']!r} vs {c['phi']!r} iters {d['iters']} vs {c['iters']} -> {'BIT-IDENTICAL' if same else 'DIFFER'}")
    else:
        print(f"  alpha={a}: default={'yes' if d else 'no'} raised={'yes' if c else 'no'} (not both yet)")
PYEOF

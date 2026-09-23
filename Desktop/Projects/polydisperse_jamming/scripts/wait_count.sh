#!/bin/bash
# Block until the jsonl files matching $1 hold at least $2 lines. Echoes what it read.
PAT=$1; N=$2
while true; do
  c=$(cat $PAT 2>/dev/null | wc -l)
  if [ "$c" -ge "$N" ]; then echo "read $c >= $N for $PAT"; break; fi
  sleep 30
done

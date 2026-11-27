#!/bin/bash
# Core budget for the long queue = 11 - A1 - C2 - D lane (4) - railgun lane (3)
cd "$(dirname "$0")"
while true; do
  a1=$(( $(grep -c "\] start" runq_A1.log) - $(grep -c "\] done" runq_A1.log) ))
  c2=$(( $(grep -c "\] start" runq_C2.log) - $(grep -c "\] done" runq_C2.log) ))
  d=0
  r=0; pgrep -f "rail2.list|railgun/coax_.*/run.sh" >/dev/null && r=3
  echo $((11 - 2*a1 - 2*c2 - d - r)) > cores_long.txt
  sleep 60
done

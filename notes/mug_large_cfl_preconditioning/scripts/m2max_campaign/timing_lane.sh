#!/bin/bash
# Exclusive timing lane: wait until all other queues and runs are finished, then run cases one at a time
cd "$(dirname "$0")"
while pgrep -f "runq.py" >/dev/null || pgrep -f "/run.sh" >/dev/null || pgrep -f "MUG_sph_tilt|MUG-railgun|test_alfven" >/dev/null; do sleep 60; done
pkill -f governor.sh
echo "[timing] start $(date)" >> runq_T.log
../../../venv/bin/python runq.py --cores 8 --mem 54 --timeout 14400 --list T.list >> runq_T.log 2>&1
[ -f Trail.list ] && ../../../venv/bin/python runq.py --cores 8 --mem 54 --timeout 14400 --list Trail.list >> runq_T.log 2>&1
echo "[timing] end $(date)" >> runq_T.log

#!/bin/bash
# SPDX-License-Identifier: LGPL-3.0-only
# Usage: ./watchdog.sh <rss_limit_kB>   (run alongside run_sph.py)
# kill MUG_sph_tilt if RSS exceeds limit (kB); record peak
lim=$1; peak=0
while pgrep -x MUG_sph_tilt >/dev/null || sleep 5; do
  r=$(ps -C MUG_sph_tilt -o rss= | sort -n | tail -1); [ -z "$r" ] && { pgrep -f run_sph.py >/dev/null || break; sleep 2; continue; }
  [ $r -gt $peak ] && peak=$r
  if [ $r -gt $lim ]; then pkill -x MUG_sph_tilt; echo "KILLED rss=$r"; fi
  echo $peak > watchdog_peak.txt; sleep 2
done

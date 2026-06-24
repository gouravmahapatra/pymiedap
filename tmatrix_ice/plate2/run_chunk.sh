#!/bin/bash
# run_chunk.sh PANEL OUTFILE TODOFILE  -- processes sizes until ~35s, then pauses
cd "$HOME/build"
P=$1; OUT=$2; TODO=$3
START=$SECONDS
while :; do
  x=$(head -1 "$TODO" 2>/dev/null)
  [ -z "$x" ] && { echo "ALL DONE for $OUT"; break; }
  echo "$x $P" | ./amplq_p2 >/dev/null 2>&1
  awk -v x=$x '{printf "%s %s %s\n", x,$1,$4}' ampl.out >> "$OUT"
  tail -n +2 "$TODO" > "$TODO.tmp" && mv "$TODO.tmp" "$TODO"
  [ $((SECONDS-START)) -ge 35 ] && { echo "paused: $(wc -l < $TODO) sizes left in $TODO"; break; }
done

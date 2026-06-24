#!/bin/bash
cd "$HOME/build"; JOBS=jobs.txt; START=$SECONDS
while :; do
  line=$(head -1 "$JOBS" 2>/dev/null); [ -z "$line" ] && { echo "ALL DONE"; break; }
  set -- $line; x=$1; eps=$2; np=$3; tag=$4
  echo "$x $eps $np" | ./tmq_p3 >/dev/null 2>&1
  awk -v x=$x 'f&&NF==7&&$1~/^[0-9]/{printf "%s %s %s\n",x,$1,(-$6/$2)} /F11.*F22.*F44/{f=1}' test >> plate3_$tag.dat
  tail -n +2 "$JOBS" > "$JOBS.tmp" && mv "$JOBS.tmp" "$JOBS"
  [ $((SECONDS-START)) -ge 3 ] && { echo "paused: $(wc -l < $JOBS) jobs left"; break; }
done

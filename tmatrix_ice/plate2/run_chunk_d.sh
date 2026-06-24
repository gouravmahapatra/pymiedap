#!/bin/bash
cd "$HOME/build"
OUT=panelD.dat; TODO=todoD.txt
START=$SECONDS
while :; do
  x=$(head -1 "$TODO" 2>/dev/null)
  [ -z "$x" ] && { echo "ALL DONE"; break; }
  echo "$x" | ./tmq_p2d >/dev/null 2>&1
  awk -v x=$x 'f&&NF==7&&$1~/^[0-9]/{printf "%s %s %s\n",x,$1,(-$6/$2)} /F11.*F22.*F44/{f=1}' test >> "$OUT"
  tail -n +2 "$TODO" > "$TODO.tmp" && mv "$TODO.tmp" "$TODO"
  [ $((SECONDS-START)) -ge 35 ] && { echo "paused: $(wc -l < $TODO) left"; break; }
done

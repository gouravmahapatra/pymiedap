#!/bin/bash
# run_amplq.sh PANEL OUTFILE XMIN XMAX XSTEP
cd "$HOME/build"
P=$1; OUT=$2; XMIN=$3; XMAX=$4; XSTEP=$5
for x in $(seq $XMIN $XSTEP $XMAX); do
  echo "$x $P" | ./amplq_p2 >/dev/null 2>&1
  awk -v x=$x '{printf "%s %s %s\n", x, $1, $4}' ampl.out >> "$OUT"
done

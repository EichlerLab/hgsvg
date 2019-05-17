#!/usr/bin/env bash

set -euo pipefail

h1=$1
h2=$2
out=$3
additional=$4
SRC_DIR=`dirname "$(readlink -f "$0")"`

set -v 
if [ $# -gt 4 ]; then
		mergeFraction=$5
else
		mergeFraction=0.5
fi


for hapfile in $h1 $h2; do 
		for op in insertion deletion; do
				cat $hapfile | bioawk -c hdr -vop=$op 'BEGIN{OFS="\t";} { if (NR==1 || $svType == op) { print; } }' | tr " " "\t"> $hapfile.$op
		done
done
for op in insertion deletion; do 
		$SRC_DIR/MergeHaplotypes.sh $h1.$op $h2.$op $out.$op "$additional" $mergeFraction
done

head -1 $out.deletion > $out
cat $out.insertion $out.deletion | \
		grep -v "^#" | bedtools sort >> $out

				



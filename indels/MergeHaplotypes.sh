#!/usr/bin/env bash

set -euo pipefail

set -v
a=$1
b=$2
r=$3
head -1 $1 | awk '{ print $0"\thap";}' > $r
bedtools intersect -header -a $a -b $b -f 0.8 -u -r -wa | awk '{ if (NR > 1) print $0"\tHOM";}' > $r.a
bedtools intersect -header -a $b -b $a -f 0.8 -u -r | awk '{ if (NR > 1) print $0"\tHOM";}' > $r.b
bedtools intersect -header -v -a $b -b $a -f 0.8  -r | awk '{ if (NR > 1) print $0"\tHAP1";}' > $r.hb
bedtools intersect -header -v -a $a -b $b -f 0.8 -r | awk '{ if (NR > 1) print $0"\tHAP0";}'> $r.ha

cat $r.a $r.hb $r.ha | bedtools sort >> $r





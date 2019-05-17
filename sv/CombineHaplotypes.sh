#!/usr/bin/env bash

set -euo pipefail

a=$1
b=$2
base=$3
bedtools intersect -a $a -b $b  -wao -f 0.5 -r | awk '{ if ($17 != ".") print;}' | bedtools groupby -c 4 -o first -full | cut -f 1-16 > $base.diploid.a.bed
bedtools intersect -a $b -b $a  -wao -f 0.5 -r | awk '{ if ($17 != ".") print;}' | bedtools groupby -c 4 -o first -full |  cut -f 1-16 > $base.diploid.b.bed
bedtools intersect -v -a $a -b $base.diploid.a.bed -wa -f 0.9 -r > $base.hap1.bed
bedtools intersect -v -a $b -b $base.diploid.b.bed -wa -f 0.5 -r > $base.hap2.bed

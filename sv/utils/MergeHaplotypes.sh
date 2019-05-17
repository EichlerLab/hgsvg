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


#
# Mark the values in haplotype 0 that are overlapping r% with hap 1
#w
nh1=`wc -l $h1 | awk '{ print $1; }'`;
nh2=`wc -l $h2 | awk '{ print $1; }'`;
echo $nh1
echo $nh2
echo $h1
echo $h2
if [ $nh1 -gt 1 ]  && [ $nh2 -gt 1 ]; then 
{ $SRC_DIR/HeaderMod.py --source $h1 $h2 --append "svOp" --index ; \
  bedtools intersect -a $h1 -b $h2  -r -f $mergeFraction -wao; } | \
  bioawk -c hdr '{ if ($chrom_2 != ".") print; }' | \
	tr -d "#" | \
  $SRC_DIR/Select.py --cols chrom tStart tEnd   $additional | \
	$SRC_DIR/SafeGroupby.sh "bedtools groupby -header -g 1-5 -c 3 -o first -full"  | \
	awk 'BEGIN{OFS="\t";} { if (NR> 1) {NF-=1}; print; }' >  $out.0.r.bed
else
		cat $h1 | tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd  $additional 		 > $out.0.r.bed
fi		

#
# Mark the values in haplotype 1 that are overlapping r% with hap 0.
# These are counted as the homozygous sites.
#
nf=`head -1 $h1 | awk '{ print NF;}'`
if [ $nh1 -gt 1 ] && [ $nh2 -gt 1 ]; then 

{ $SRC_DIR/HeaderMod.py --source $h2 $h1 --append "svOp" --index ; \
  bedtools intersect -b $h1 -a $h2 -r -f $mergeFraction -wao ;} | \
  bioawk -c hdr '{ if ($chrom_2 != ".") print; }' | \
  tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd $additional | \
  	$SRC_DIR/SafeGroupby.sh " bedtools groupby -header -c 3 -o first -full " | \
	awk 'BEGIN{OFS="\t";} {  if (NR > 1) {NF-=1}; print; }' > $out.1.r.bed
else
		echo "cat $h1 | $SRC_DIR/Select.py --cols chrom tStart tEnd  $additional 		 > $out.1.r.bed"
		cat $h1 | tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd  $additional 		 > $out.1.r.bed
fi		

#
# Remove calls marked as homozygous from hap0
#
nHom=`wc -l $out.0.r.bed | awk '{ print $1;}'`
nHap=`wc -l $h1 | awk '{ print $1;}'`

if [ $nHom -gt 1 ] && [ $nHom -lt $nHap ]; then
nf=`head -1 $h1 | awk '{ print NF;}'`
  bedtools intersect -v -a $h1 -b $out.0.r.bed  -header -r -f 1.0 | \
		 $SRC_DIR/SafeGroupby.sh "bedtools groupby -header -g 1-3 -c 1 -o first -full " | cut -f 1-$nf > $out.h0.bed
else
		head -1 $h1 > $out.h0.bed
fi		

#
# Remove calls marked as homozygous from hap1
#

nHom=`wc -l $out.1.r.bed | awk '{ print $1;}'`
nHap=`wc -l $h2 | awk '{ print $1;}'`
if [ $nHom -gt 1 ] && [ $nHom -lt $nHap ] ; then
nf=`head -1 $h2 | awk '{ print NF;}'`
	bedtools intersect  -v -a $h2 -b $out.1.r.bed  -header -r -f 1.0 | \
    $SRC_DIR/SafeGroupby.sh "	bedtools groupby -header -g 1-3 -c 1 -o first -full " | cut -f 1-$nf > $out.h1.bed
else
		head -1 $h2 > $out.h1.bed
fi



cp $out.1.r.bed $out.hom.bed
head -1 $out.h0.bed | awk '{ print $0"\thap"; }' > $out.h0.bed.hap
grep -v "^#" $out.h0.bed | awk '{ print $0"\tHAP0";}' >> $out.h0.bed.hap
head -1 $out.h1.bed | awk '{ print $0"\thap"; }' > $out.h1.bed.hap
grep -v "^#" $out.h1.bed | awk '{ print $0"\tHAP1";}' >> $out.h1.bed.hap
head -1 $out.hom.bed | awk '{ print $0"\thap"; }' > $out.hom.bed.hap
grep -v "^#" $out.hom.bed | awk '{ print $0"\tHOM";}' >> $out.hom.bed.hap


{ cat $out.h0.bed.hap | tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd hap $additional ; \
  cat $out.h1.bed.hap | tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd hap $additional | grep -v "^#" ; \
  cat $out.hom.bed.hap | tr -d "#" | $SRC_DIR/Select.py --cols chrom tStart tEnd hap $additional | grep -v "^#" ; } | \
		bedtools sort -header > $out

#	rm -f $out.0.r.bed
#	rm -f $out.h0.bed
#	rm -f $out.h1.bed
#	rm -f $out.1.r.bed
#	rm -f $out.hom.bed
#	rm -f $out.h0.bed.hap $out.h1.bed.hap $out.hom.bed.hap

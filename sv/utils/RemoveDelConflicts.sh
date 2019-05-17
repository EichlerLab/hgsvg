#!/usr/bin/env bash

set -euo pipefail

cat $1 | egrep "^#|deletion"  > $1.all_del
cat $1 | egrep "^#|deletion" | egrep "^#|HAP0|HOM" > $1.hap0_del
cat $1 | egrep "^#|deletion" | egrep "^#|HAP1|HOM" > $1.hap1_del
scripts=$(dirname $0)


cat $1 | egrep "^#|insertion" | $scripts/ToPoint.sh | egrep "^#|HOM" | bedtools intersect -header -v -a stdin -b $1.all_del | $scripts/FromPoint.sh > $1.no_hom_conflict

cat $1 | egrep "^#|insertion" | $scripts/ToPoint.sh | egrep "^#|HAP0" | bedtools intersect -header -v -a stdin -b $1.hap0_del | $scripts/FromPoint.sh > $1.no_hap0_conflict

cat $1 | egrep "^#|insertion" | $scripts/ToPoint.sh | egrep "^#|HAP1" | bedtools intersect -header -v -a stdin -b $1.hap1_del | $scripts/FromPoint.sh > $1.no_hap1_conflict

cat $1.all_del | bioawk -c hdr '{ if (NR==1 || $hap == "HOM") print;}' | $scripts/SimpleRMDup.py > $1.hom_del

cat $1.all_del | bioawk -c hdr '{ if ($hap != "HOM") print;}' >  $1.hap_del


cat <( head -1 $1) <( grep -v "^#" $1.hom_del ) <( grep -v "^#" $1.hap_del ) <( grep -v "^#" $1.no_hom_conflict ) <( grep -v "^#" $1.no_hap0_conflict ) <( grep -v "^#" $1.no_hap1_conflict ) | tr " " "\t"  | bedtools sort -header > $2

#rm -f $1.hom_del $1.hap0_del $1.hap1_del $1.no_hap0_conflict $1.no_hap1_conflict $1.no_hom_conflict



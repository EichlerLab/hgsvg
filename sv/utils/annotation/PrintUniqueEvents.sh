#!/usr/bin/env bash
BASE=$(dirname "$0")
calls=$1
outDir=$2
mkdir -p $outDir
set -v
#cutoffs for repeat content determined from manual inspection of results
cat $calls | bioawk -c hdr '{ if (NR==1 || ($fracTR >= 0.8)) {if (NR > 1) { $svAnn = "TandemRepeat";}} print;}' | tr " " "\t"> $calls.TandemAnnotated
cat $calls.TandemAnnotated | bioawk -c hdr '{ if (NR==1 || $svAnn != "TandemRepeat") print;}' > $calls.not_trf
cat $calls.TandemAnnotated | bioawk -c hdr '{ if (NR==1 || $svAnn == "TandemRepeat") print;}' > $outDir/TandemRepeat.bed
cat $calls.not_trf | $BASE/FixMasked.py | bioawk -c hdr '{ if (NR==1||$svRep >= 0.7) print;}' > $calls.repeat
cat $calls.not_trf | $BASE/FixMasked.py | bioawk -c hdr '{ if (NR==1||$svRep < 0.7) print;}' > $calls.not_repeat

$BASE/PrintUniqueEvents.py $1.repeat  --prefix AluY --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $outDir/1.bed > $outDir/AluY.simple.bed
$BASE/PrintUniqueEvents.py $outDir/1.bed --prefix AluS --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $outDir/2.bed > $outDir/AluS.simple.bed
$BASE/PrintUniqueEvents.py $outDir/2.bed --minSTR 1  --maxNotPrefix 0 --remainder $outDir/4.bed > $outDir/STR.bed
$BASE/PrintUniqueEvents.py $outDir/4.bed --prefix L1HS  --maxNotPrefix 0 --remainder $outDir/5.bed > $outDir/L1HS.simple.bed
$BASE/PrintUniqueEvents.py $outDir/5.bed --prefix Alu  --minPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $outDir/6.bed > $outDir/Alu.Mosaic.bed
$BASE/PrintUniqueEvents.py $outDir/6.bed --prefix Alu  --minSTR 1 --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/7.bed > $outDir/Alu.STR.bed
$BASE/PrintUniqueEvents.py $outDir/7.bed --prefix ALR   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/8.bed > $outDir/ALR.bed
$BASE/PrintUniqueEvents.py $outDir/8.bed --prefix SVA   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/9.bed > $outDir/SVA.simple.bed
$BASE/PrintUniqueEvents.py $outDir/9.bed --prefix HERV   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/10.bed > $outDir/HERV.simple.bed
$BASE/PrintUniqueEvents.py $outDir/10.bed --prefix L1P   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/11.bed > $outDir/L1P.bed
$BASE/PrintUniqueEvents.py $outDir/11.bed --prefix BSR/Beta   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/12.bed > $outDir/Beta.bed
$BASE/PrintUniqueEvents.py $outDir/12.bed --prefix HSAT   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/13.bed > $outDir/HSAT.bed
$BASE/PrintUniqueEvents.py $outDir/13.bed --prefix MER   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/14.bed > $outDir/MER.bed
$BASE/PrintUniqueEvents.py $outDir/14.bed --prefix L1   --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/15.bed > $outDir/L1.bed
$BASE/PrintUniqueEvents.py $outDir/15.bed --prefix LTR  --minPrefix 1 --maxNotPrefix 0 --remainder $outDir/16.bed > $outDir/LTR.bed
$BASE/PrintUniqueEvents.py $outDir/16.bed --max 1 --remainder $outDir/17.bed > $outDir/Singletons.bed

mv $outDir/17.bed $outDir/Complex.bed
mv $calls.not_repeat $outDir/Complex.NotRepeat.bed
rm -f $outDir/[0-9]*.bed

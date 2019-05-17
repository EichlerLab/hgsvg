#!/usr/bin/env bash

set -euo pipefail

samFile=$1
dotsFile=$2
gapsFile=$3
pdfFile=$4
refRegion=$5
queryRegion=$6
/net/eichler/vol5/home/mchaisso/projects/mcutils/bin/samToDotPlot $samFile > $dotsFile
/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/utils/GapBedToBed6.py $gapsFile $gapsFile.bed6

Rscript /net/eichler/vol5/home/mchaisso/projects/mcutils/src/RenderPlot.R --dots $dotsFile --output $pdfFile  --targettrack $gapsFile.bed6 --xlabel $queryRegion --ylabel  $refRegion --width 2

rm -f $gapsFile.ins
rm -f $gapsFile.del






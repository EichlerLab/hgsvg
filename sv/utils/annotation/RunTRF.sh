#!/usr/bin/env sh
BASE=$(dirname "$0")

$BASE/GapBedToFasta.py $1 $1.fasta
~/software/bin/trf $1.fasta  2 7 7 80 10 20 500 -m -ngs -h  > $1.trf_annot

$BASE/AnnotateWithTRF.py $1 $1.trf_annot $2

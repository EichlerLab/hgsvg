#!/usr/bin/env bsah

REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py $REF $1 --flank 1000 --flankIdentity 0.70 --minFraction 0.6 --minLength 50 --outFile gaps.bed


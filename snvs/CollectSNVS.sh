#!/usr/bin/env bash

set -euo pipefail

for bam  in `cat $1 `; do
		base=$(basename $bam)
		/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/Phasing/BamToNucFreq.sh $bam $2 $base.vcf /net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
done

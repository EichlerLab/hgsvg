#!/usr/bin/env bash

set -euo pipefail

mkdir -p hap1
mkdir -p hap2
/net/eichler/vol5/home/mchaisso/projects/HGSVG/scripts/FilterSamByHaplotype.py hap_samfiles.fofn 1 | samtools view -bS - | samtools sort -T $TMPDIR -o hap1/alignments.bam
samtools index hap1/alignments.bam

/net/eichler/vol5/home/mchaisso/projects/HGSVG/scripts/FilterSamByHaplotype.py hap_samfiles.fofn 2 | samtools view -bS - | samtools sort -T $TMPDIR -o hap2/alignments.bam
samtools index hap2/alignments.bam

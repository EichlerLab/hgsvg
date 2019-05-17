#!/usr/bin/env bash

set -euo pipefail

bedtools sort -i inversions.h0.bed  | bedtools merge > inversions.h0.merged.bed
bedtools sort -i inversions.h1.bed  | bedtools merge > inversions.h1.merged.bed

bedtools intersect -a inversions.h0.merged.bed -b inversions.h1.merged.bed -wa > inversions.overlap.h0.bed
bedtools intersect -b inversions.h0.merged.bed -a inversions.h1.merged.bed -wa > inversions.overlap.h1.bed

bedtools intersect -a inversions.h0.merged.bed -b inversions.overlap.h0.bed -v > inversions.h0-only.bed
bedtools intersect -a inversions.h1.merged.bed -b inversions.overlap.h1.bed -v > inversions.h1-only.bed

cat inversions.h0-only.bed | awk '{ print $1"\t"$2"\t"$3"\tHET-H0";}' > $1.inversions.haplotype-annotated.bed
cat inversions.h1-only.bed | awk '{ print $1"\t"$2"\t"$3"\tHET-H1";}' >> $1.inversions.haplotype-annotated.bed
cat inversions.overlap.h0.bed | awk '{ print $1"\t"$2"\t"$3"\tHOM";}' >> $1.inversions.haplotype-annotated.bed


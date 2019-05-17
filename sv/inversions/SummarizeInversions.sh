#!/usr/bin/env bash

set -euo pipefail

bedtools sort -i inversions.h0.bed | bedtools merge > inversions.h0.merged.bed
bedtools sort -i inversions.h1.bed | bedtools merge > inversions.h1.merged.bed


bedtools intersect -a inversions.h0.merged.bed -b inversions.h1.merged.bed -v -wa  | awk '{ print $3-$2;}' | stats.py > inv.h0.stats
bedtools intersect -b inversions.h0.merged.bed -a inversions.h1.merged.bed -v -wa | awk '{ print $3-$2;}' | stats.py > inv.h1.stats
bedtools intersect -a inversions.h0.merged.bed -b inversions.h1.merged.bed -wa |  awk '{ print $3-$2;}' | stats.py > inv.hom.stats


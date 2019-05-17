#!/usr/bin/env bash

set -euo pipefail

/bin/rm -f unaligned.fofn
mkdir -p unaligned
for f in `cat $1`; do
    b=$(basename $f)
    n=`echo $b | tr "\." "\t" | cut -f 2`
    samtools view -f 4 $f -o unaligned/unaligned.$n.bam
    echo "unaligned/unaligned.$n.bam" >> unaligned.fofn
    samtools index unaligned/unaligned.$n.bam
done

samtools merge -1 -@2 -b unaligned.fofn unaligned.bam
samtools index unaligned.bam



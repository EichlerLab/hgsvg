#!/usr/bin/env python

import sys
if len(sys.argv) != 3:
    print("Usage: SplitGapsByHaplotype.py gaps.bed out_base")
    sys.exit(0)
gapsFile = open(sys.argv[1])
h1Name = sys.argv[2] + ".1.bed"
h2Name = sys.argv[2] + ".2.bed"
h1Out = open(h1Name, 'w')
h2Out = open(h2Name, 'w')

for line in gapsFile:
    v = line.split()
    contig = v[7]
    hap = contig.split("/")[2]
    if hap == "1":
        h1Out.write(line)
    elif hap == "2":
        h2Out.write(line)



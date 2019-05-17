#!/usr/bin/env python

import sys
import pysam
import Overlap
import re
if len(sys.argv) != 4:
	print("usage: FixOverlaps.py overlaps.txt assemblies.fasta overlaps.out.txt")
	sys.exit(0)

overlaps = Overlap.ReadOverlapFile(sys.argv[1])
assemblies = pysam.FastaFile(sys.argv[2])
outputOverlaps = open(sys.argv[3], 'w')
nre=re.compile("^(N*)[^N]*(N*)$")

def GetN(seq):
    m= nre.match(seq)
    if m is None:
        return [0,len(seq)]
    else:
        g =m.groups()
        return [len(g[0]), len(seq) - len(g[1])]

for ovp in overlaps:
    try:
        aSeq = assemblies.fetch(ovp.a)
        bSeq = assemblies.fetch(ovp.b)
    except IndexError:
        sys.stderr.write("WARNING one of " + ovp.a + " and " + ovp.b + " was missing from the contig list\n")
        continue

    ovp.aRead = GetN(aSeq)
    ovp.bRead = GetN(bSeq)

    ovp.Print(outputOverlaps)
outputOverlaps.close()

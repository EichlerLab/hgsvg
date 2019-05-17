#!/usr/bin/env python
import argparse
import pysam
import os
import sys

ap = argparse.ArgumentParser(description="Given a list of assembly sam files, print records into distinct sam files")
ap.add_argument("--alignments", help="FOFN of alignments.")
ap.add_argument("--header", help="Header to use.")
ap.add_argument("--base", help="Base name for haplotypes", default="alignments")

args = ap.parse_args()

headerFile = open(args.header)
header=headerFile.readlines()
h0=open(args.base + ".h0.sam", 'w')
h1=open(args.base + ".h1.sam", 'w')
h0.write(''.join(header))
h1.write(''.join(header))
alignments = open(args.alignments)
for alnFileName in alignments:
    if os.stat(alnFileName.rstrip()).st_size == 0:
        continue
    alnFile = open(alnFileName.rstrip())
    
    for aln in alnFile:
        vals=aln.split()
        hap = None
        if vals[0].count("/") == 2:
            titleVals = vals[0].split("/")
            if titleVals[2] == "1":
                h0.write(aln)
                hap=0
            elif  titleVals[2] == "2":
                h1.write(aln)
                hap=1
            sys.stderr.write(vals[0]  + "\t" + str(hap) + "\n")
        else:
            for v in vals:
                if v == "HA:i:1":
                    hap=0
                    h0.write(aln)
                    break                
                elif v == "HA:i:2":
                    hap=1
                    h1.write(aln)
                    continue
        if hap is None:
            sys.stderr.write("ERROR, did not find haplotype for " + vals[0] + "\n")
            



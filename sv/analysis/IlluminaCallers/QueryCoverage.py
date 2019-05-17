#!/usr/bin/env python

import argparse
import tempfile
import os
import sys
import pysam
import subprocess
from multiprocessing import Process, Lock, Semaphore, Pool

ap = argparse.ArgumentParser(description="Determine how many pb reads overlap an Illumina call")
ap.add_argument("--fofn", help="FOFN of pb-bams", required=True)
ap.add_argument("--region", help="Region for coverage")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of threads.", default=1,type=int)
ap.add_argument("--op", help="Guide coverage by operation", default="DEL")
ap.add_argument("--bionano", help="Compute bionano estimage of coverage -- min coverage in window specified by svLen for deletions, and max for insertions", action="store_true", default="False")
ap.add_argument("--header", help="Write header", default=False,action="store_true")
ap.add_argument("--svlen", help="Write average coverage in this length.", default=None,type=int)
args = ap.parse_args()

tmpdir="."
if "TMPDIR" in os.environ:
	 tmpdir=os.environ["TMPDIR"]
readsFofn = open(args.fofn)
blines = readsFofn.readlines()
bamFiles = [pysam.AlignmentFile(bline.rstrip(), 'rb') for bline in blines]
bamFileIndex = {}
for b in range(0,len(blines)):
    vals = blines[b].split(".")
    for v in vals:
        if "chr" in v:
            bamFileIndex[v]= b
            break

if args.header:
    outFile.write("coverage\n")

def ProcessRegion(rgn):
    (regionChrom,coords) = rgn.split(":")
    cstr = coords.split("-")
    
    regionStart = int(cstr[0])
    regionEnd   = int(cstr[1])

    
    tempFileNames = []
    regionLength = regionEnd-regionStart
    cov = [0]*(regionLength)

    if regionChrom not in bamFileIndex:
        sys.stderr.write("Missing " + regionChrom + "\n")
        return 0
    b = bamFileIndex[regionChrom]


    covall = bamFiles[b].count_coverage(regionChrom, regionStart, regionEnd, quality_threshold=0)

    for i in range(0,len(covall[0])):
        cov[i] = covall[0][i] + covall[1][i] + covall[2][i]  + covall[3][i]

    if args.svlen is None:
        print(cov)
    else:
        if len(cov) < args.svlen:
            print("{}\t{}\t{}\t{:2.2f}".format(regionChrom, regionStart, regionStart+1, sum(cov)/args.svlen))
        else:
            wCov = sum(cov[0:args.svlen])
            print("{}\t{}\t{}\t{:2.2f}".format(regionChrom, regionStart, regionStart+1, wCov/float(args.svlen)))
            for i in range(1,len(cov)-args.svlen-1):
                wCov-=cov[i-1]
                wCov+=cov[i+args.svlen-1]
                print("{}\t{}\t{}\t{:2.2f}".format(regionChrom, regionStart+i, regionStart+i+1, wCov/float(args.svlen)))
                
#    print "{:2.2f}".format(float(total)/len(cov))

ProcessRegion(args.region)

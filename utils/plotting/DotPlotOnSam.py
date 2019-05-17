#!/usr/bin/env python

import sys
import subprocess
import bisect
import re
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
if (len(sys.argv) < 6):
    print("usage DotPlotOnSam.py samFile chrom start end ref ")
    sys.exit(0)
    
samFile = open(sys.argv[1])
chrom = sys.argv[2]
start = int(sys.argv[3])
end   = int(sys.argv[4])
ref   = sys.argv[5]
l = end-start
rgnStart = start-l
rgnEnd   = end+l
center = (end+start)/2


index = 0
for line in samFile:
    vals = line.split()
    ctgChrom = vals[2]
    ctgStart = int(vals[3])

    cigar = vals[5]
    cigarOps = re.findall(r'(\d+)(\w)', cigar)
    query = vals[9]
    qOffset = start - ctgStart 
    rPos = [0]*len(query)
    r = ctgStart
    qi = 0
    refLen = 0
    soft = 0
    for i in range(0,len(cigarOps)):
        cigarOp = (int(cigarOps[i][0]), cigarOps[i][1])

        if (cigarOp[1] == 'S'):
            qi += cigarOp[0]
            if (i == 0):
                soft = cigarOp[0]
        elif (cigarOp[1] == 'M' or cigarOp[1] == '=' or cigarOp[1] == 'X'):
            for j in range(0,cigarOp[0]):
                rPos[qi] = r
                r+=1
                qi+=1
            refLen += cigarOp[0]
        elif (cigarOp[1] == 'D'):
            r += cigarOp[0]
        elif (cigarOp[1] == 'I'):
            for j in range(0,cigarOp[0]):
                rPos[qi] = r
                qi+=1
            refLen+= cigarOp[0]
    
    # Now find the region to dotplot
    refStart = bisect.bisect_left(rPos, rgnStart)
    refEnd   = bisect.bisect_left(rPos, rgnEnd)
    ctgEnd = ctgStart + refLen


    qStart=max(0,soft+qOffset-l)
    qEnd  = min(len(query), soft+qOffset+3*l)
    qSeq = query[qStart:qEnd]
    print(qOffset)
    print(l)
    print(str(qStart) + "\t" + str(qEnd))
    rec = SeqRecord.SeqRecord(Seq.Seq(qSeq), id="Inversion", name="", description="")
    tempInv= open("tmp.inv.fasta", 'w')
    SeqIO.write(rec, tempInv, "fasta")
    refFile = open("ref.fasta", 'w')

    getRef = "samtools faidx {} {}:{}-{} ".format(ref, chrom, center-3*l, center+3*l)
    p = subprocess.Popen(getRef.split(), stdout=refFile)
    p.wait()
    dots = open("ref.dots", 'w')
    dotPlot = "/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/dotPlot tmp.inv.fasta ref.fasta 15"
    p = subprocess.Popen(dotPlot.split(), stdout=dots)
    p.wait()
    renderPlot = "Rscript /net/eichler/vol5/home/mchaisso/projects/mcutils/src/RenderPlot.R -d ref.dots --xlabel {}:{}-{}.{} --ylabel {}:{}-{} --output {}.{}-{}.{}.pdf".format(chrom,start,end,index,chrom,start,end, chrom,start,end,index)
    p = subprocess.Popen(renderPlot.split())
    p.wait()
    index+=1

    

#!/usr/bin/env python

import pysam
import tempfile
import subprocess
import argparse
import re
import os
import sys
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO

ap = argparse.ArgumentParser(description="Realign regions of SV")
ap.add_argument("--asm", help="Assembled genome", required=True)
ap.add_argument("--ref", help="Target genome", required=True)
ap.add_argument("--asmRegion", help="Assembled region", required=True)
ap.add_argument("--refRegion", help="Target region", required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--window", help="Region to align around window", type=int, default=10000)
ap.add_argument("--tmpdir", help="Where to place temporary files.", default=None)
ap.add_argument("--index", help="Index of SV.",default=0)
args = ap.parse_args()

if args.tmpdir is None:
    if "TMPDIR" not in os.environ or os.environ["TMPDIR"] == "":
        print("ERROR. The TEMPDIR variable must be set or --tmpdir specified on as a command  argument")
        sys.exit(1)
    else:
        args.tmpdir = os.environ["TMPDIR"]

asmFile = pysam.FastaFile(args.asm)
refFile = pysam.FastaFile(args.ref)


asmFaiFile = open(args.asm + ".fai")
asmFai = {l.split()[0]: int(l.split()[1]) for l in asmFaiFile}

refFaiFile = open(args.ref + ".fai")
refFai = {l.split()[0]: int(l.split()[1]) for l in refFaiFile}

regionRe=re.compile("([^:]*):(\d+)-(\d+)")
def GetRegion(regionStr):
    null=(None, None, None)

    m = regionRe.match(regionStr)

    if m is None:
        return null
    else:
        g=m.groups()

        if len(g) != 3:
            return null
        else:
            return (g[0], int(g[1]), int(g[2]))

def ExpandRegion(contig, start, end, fai, window):
    start = max(0, start - window)
    end   = min(fai[contig], end + window)
    return (start, end)

(aContig, aStart, aEnd) = GetRegion(args.asmRegion)
(rContig, rStart, rEnd) = GetRegion(args.refRegion)

args.window = max(args.window, 2*(max(aEnd-aStart, rEnd - rStart)))

(aStartExp, aEndExp) = ExpandRegion(aContig, aStart, aEnd, asmFai, args.window)
(rStartExp, rEndExp) = ExpandRegion(rContig, rStart, rEnd, refFai, args.window)
aSeq = asmFile.fetch(aContig, aStartExp, aEndExp)
rSeq = refFile.fetch(rContig, rStartExp, rEndExp)

aFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".fasta", delete=False, mode='w')
rFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".fasta", delete=False, mode='w')
sFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".sam", delete=False, mode='w')

def WriteSeq(fh,seq,seqName):
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq),id=seqName,name="",description=""),fh,"fasta")

WriteSeq(aFile, aSeq, "asm")
WriteSeq(rFile, rSeq, "ref")
aFile.close()
rFile.close()
blasrCmd="/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/blasr"
alnCommand="{} {} {} -sam -bestn 1 -sdpTupleSize 13 -out {}  ".format(blasrCmd, aFile.name, rFile.name, sFile.name)
devnull=open(os.devnull)
subprocess.call(alnCommand.split(), stderr=devnull)
#sys.stderr.write( alnCommand + "\n")
samFile = pysam.AlignmentFile(sFile.name, "r")
# There should only be one alignment
minRatio = 1.0
delta = (aEnd - aStart ) - (rEnd - rStart)

for aln in samFile.fetch():
    # negative delta is deletion, positive is insertion


#    print((aEnd, aStart, rEnd , rStart))
    verified=False
    nMapped=0
    if aln.cigartuples is None:
        break
#    sys.stderr.write(str((aln.query_alignment_start, aln.query_alignment_end, aln.reference_end, aln.reference_start))+"\n")
    nMapped = aln.reference_end-aln.reference_start
    for t in aln.cigartuples:
        
        if t[0] == 1 or t[0] == 2:
            oplen = t[1]
            if t[0] == 2:
                oplen = -oplen

            ratio = abs(1-delta / float(oplen))
            minRatio = min(minRatio, ratio)
            if  ratio < 0.8:
                verified=True
                break
    sys.stderr.write("Index {} mapped {}\n".format(args.index, nMapped))
command = "rm -f {} {} {}".format(aFile.name, rFile.name, sFile.name)
subprocess.call(command.split())
print(args.asmRegion + "\t" + str(verified) + "\t{}\t{:2.2f}".format(abs(delta), minRatio) + "\t" + str(nMapped) + "\t" + str(len(aSeq)) + "\t" + str(len(rSeq)) + "\t" + str(args.index))

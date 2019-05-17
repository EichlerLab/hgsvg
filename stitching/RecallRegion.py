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
ap.add_argument("--asm", help="Assembled genome", required=False)
ap.add_argument("--ref", help="Target genome", required=True)
ap.add_argument("--asmRegion", help="Assembled region", required=True)
ap.add_argument("--refRegion", help="Target region", required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--window", help="Region to align around window", type=int, default=10000)
ap.add_argument("--tmpdir", help="Where to place temporary files.", default=None)
ap.add_argument("--ngaps", help="Number of gaps originally in alignment", default=None, type=int)
ap.add_argument("--gapLines", help="Original gap lines", default=None)
ap.add_argument("--index", help="Index of SV.",default=0)
ap.add_argument("--keep", help="Keep temporary files", default=False, action='store_true')
ap.add_argument("--dotplot", help="Produce a dotplot of this region.", default=None)
ap.add_argument("--indels", help="Store indels from this region.", default=None)
ap.add_argument("--ngmlr", help="Use the ngmlr method to reacall", default=0,type=int)
ap.add_argument("--header", help="Output values from this header.", default=None, required=True)
ap.add_argument("--minAlignLength", help="Call gaps if the alignment of the regions is at least this length.", type=int, default=0)
ap.add_argument("--blasr", help="default blasr to run, defaults to system", default="blasr")
args = ap.parse_args()

if args.tmpdir is None:
    if "TMPDIR" not in os.environ or os.environ["TMPDIR"] == "":
        print "ERROR. The TEMPDIR variable must be set or --tmpdir specified on as a command  argument"
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

useHeader = args.header.split("__HEADER__")


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

#args.window = max(args.window, 2*(max(aEnd-aStart, rEnd - rStart)))
rRegion="{}:{}-{}".format(rContig,rStart,rEnd)
aRegion="{}:{}-{}".format(aContig,aStart,aEnd)

(aStartExp, aEndExp) = (aStart, aEnd) 
(rStartExp, rEndExp) = (rStart, rEnd) 

if aStartExp < 0 or rStartExp < 0:
    sys.stderr.write("ERROR with regions" + str((aStartExp,aEndExp, rStartExp, rEndExp)))

aSeq = asmFile.fetch(aContig, aStartExp, aEndExp)
rSeq = refFile.fetch(rContig, rStartExp, rEndExp)
aFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".asm.fasta", delete=False, mode='w')
rFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".ref.fasta", delete=False, mode='w')
sFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".sam", delete=False, mode='w')

def WriteSeq(fh,seq,seqName):
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq),id=seqName,name="",description=""),fh,"fasta")

WriteSeq(aFile, aSeq, "asm")
WriteSeq(rFile, rSeq, "ref")
aFile.close()
rFile.close()
sFile.close()
method = "blasr"
scriptDir=os.path.dirname(os.path.realpath(sys.argv[0]))
if args.ngaps > args.ngmlr and len(aSeq) < 100000 and len(rSeq) < 100000:
    alnCommand=scriptDir+"/RunNGMLR.sh {} {} {}  ".format(aFile.name, rFile.name, sFile.name)
    method = "ngmlr"
else:
    alnCommand="{} {} {} -indelRate 3 -sdpTupleSize 7 -sdpIns 5 -sdpDel 5  -maxAnchorsPerPosition 10 -maxMatch 25 -bestn 1 -sam -out {}  ".format(args.blasr, aFile.name, rFile.name, sFile.name)

 
devnull=open(os.devnull)
subprocess.call(alnCommand.split(), stderr=devnull)
#sys.stderr.write( alnCommand + "\n" + "\n")
#samFile = pysam.AlignmentFile(sFile.name, "r")
# There should only be one alignment
minRatio = 1.0
gFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".bed", delete=False, mode='w')

faiCommand="samtools faidx {}".format(rFile.name)

subprocess.call(faiCommand.split())
#gFile.close()
pgCommand=scriptDir+"/../sv/utils/PrintGaps.py {} {} --outFile {} --condense 20 --maxMasked 10 ".format(rFile.name, sFile.name, gFile.name)
subprocess.call(pgCommand.split())


if args.dotplot is not None:
    dotsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".dots", delete=False, mode='w')
    dotsFile.close()    
    dotPlotCommand = scriptDir+"/stitching/DotPlotRealignedRegion.sh {} {} {} {} {} {}".format(sFile.name, dotsFile.name, gFile.name, args.dotplot, rRegion, aRegion )
    sys.stderr.write(dotPlotCommand + "\n")
    subprocess.call(dotPlotCommand.split())



if args.indels is not None:
    tempIndelsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".bed", delete=False, mode='w')
    pgIndelsCommand= scriptDir+"/../sv/utils/PrintGaps.py {} {} --outFile {} --condense 20 --maxMasked 10 --maxLength 50 --minLength 2 --minAlignmentLength {} ".format(rFile.name, sFile.name, tempIndelsFile.name, args.minAlignLength)
    sys.stderr.write(pgIndelsCommand+"\n\n")
    subprocess.call(pgIndelsCommand.split())
    tempIndelsFile.close()
    sys.stderr.write( "opening indels file " + args.indels + "\n")
    indelFile = open(args.indels,'w')
    tmpIndelFile = open(tempIndelsFile.name)
    for line in tmpIndelFile:
        if line[0] == '#':
            indelFile.write(line)
        else:
            gap = line.split()
            gap[0] = rContig
            gap[1] = str(int(gap[1]) + rStart)
            gap[2] = str(int(gap[2]) + rStart)
            gap[7] = aContig
            gap[8] = str(int(gap[8]) + aStart)
            gap[9] = str(int(gap[9]) + aStart)
            indelFile.write("\t".join(gap) + "\n")
    
    command = "rm -f {} ".format(tempIndelsFile.name )
    subprocess.call(command.split())

gaps = ""
pgFile = open(gFile.name)
nGaps =0
header={}
cols={}
for gapLine in pgFile:
    gap = gapLine.rstrip().split()

    #
    # Compute the coordinate of the SV in the context of the genome and not just the local alignment.  You would want to fix strand before this.
    #
    if gap[0][0] == "#":
        headerVals = gapLine[1:].split()
        header = {headerVals[i]: i for i in range(0,len(headerVals))}
        continue        
    if gap[0][0] != "#":
        #
        #  Check for reverse strand, and flip coordinates if so
        #
        if "strand" in header and gap[header["strand"]] == "1":
            qStart = int(gap[8])
            qEnd   = int(gap[9])
            qLen = aEndExp - aStartExp
            gap[8] = str(qLen - qEnd)
            gap[9] = str(qLen - qStart)
            
        gap[0] = rContig
        gap[1] = str(int(gap[1]) + rStartExp)
        gap[2] = str(int(gap[2]) + rStartExp)
        gap[7] = aContig
        gap[8] = str(int(gap[8]) + aStartExp)
        gap[9] = str(int(gap[9]) + aStartExp)
    gaps +="\t".join(gap) +"\n"
    nGaps+=1
#sys.stderr.write("recalling got " + str(npg) + "\n")
if args.keep is False:
    command = "rm -f {} {} {}.fai {} {}".format(aFile.name, rFile.name, rFile.name, sFile.name, gFile.name)
    subprocess.call(command.split())
else:
    sys.stderr.write("input: " + str(args.ngaps) + " called " + str(nGaps) + "\n")
    sys.stderr.write(alnCommand)


if nGaps <= args.ngaps or args.gapLines is None:
    for gap in gaps.split("\n"):
        printGapLineVals = []
        vals = gap.split("\t")
        if len(gap) > 0:
            for i in range(0,len(useHeader)):
                if useHeader[i] in header:
                    printGapLineVals.append(vals[header[useHeader[i]]])
                else:
                    printGapLineVals.append("NA")

            gap="\t".join(printGapLineVals)
            if args.ngaps >= args.ngmlr:
                print gap + "\t" + method
            else:
                print gap + "\t" + method
                    

else:
    # The number of recalled gaps is /greater/ than the original calls
    sys.stderr.write("using prev calls\n")
    if args.gapLines is not None:
        sys.stderr.write(args.gapLines+ "\n")
        gapLinesFile = open(args.gapLines)
        gapLines = ";".join(gapLinesFile.readlines())
        gapLines = gapLines.replace("%", "\t")
        gapLines = gapLines.replace(";", "\n")
        allGapLines = gapLines.split("\n")

        gapLineHeader = {}
        if len(allGapLines) > 0 and allGapLines[0][0] == "#":
            gapLineHeaderVals = allGapLines[0][1:].split("\t")
            gapLineHeader = { gapLineHeaderVals[i] : i for i in range(0,len(gapLineHeaderVals)) }
        for line in allGapLines:
            if len(line) > 0 and line[0] == "#":
                continue
            vals = line.split()
            outputVals = []
            for i in range(0,len(useHeader)):
                if useHeader[i] in gapLineHeader:
                    outputVals.append(vals[gapLineHeader[useHeader[i]]])
                else:
                    outputVals.append("NA")
            
            print "\t".join(outputVals) + "\tgenome-align"

sys.stderr.write("Done with " + args.refRegion + "\n")

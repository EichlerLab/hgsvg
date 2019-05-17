#!/usr/bin/env python

import argparse
import subprocess
import sys
import os

from multiprocessing import Process, Lock, Semaphore, Pool
import tempfile


ap = argparse.ArgumentParser(description="Realign regions in a gap bed file")
ap.add_argument("--asm", help="Assembled genome", required=False, default=None)
ap.add_argument("--ref", help="Target genome", required=False, default=None)
ap.add_argument("--gaps", help="Gaps bed file.", required=True)
ap.add_argument("--split", help="Just split input into N files, do not align.", type=int, default=None)
ap.add_argument("--splitDir", help="Write split input files here.", default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--window", help="Region to align around window", type=int, default=20000)
ap.add_argument("--tmpdir", help="Where to place temporary files.", default=None)
ap.add_argument("--nproc", help="Number of threads.", default=8,type=int)
ap.add_argument("--keep", help="Keep temporary files", default=False,action='store_true')
ap.add_argument("--refRegions", help="Write reference regions here", default=None)
ap.add_argument("--ngmlr", help="Use ngmlr to realign region.",default=None, type=int)
ap.add_argument("--indels", help="Store indels here.", default=None)
ap.add_argument("--indelDir", help="Store indels here.", default="indels")
ap.add_argument("--commands", help="Write commands and exit", default=None)
ap.add_argument("--blasr", help="default blasr to run, defaults to system", default="blasr")
args = ap.parse_args()

gapFile= open(args.gaps)
outFile = open(args.out,'w')
i=0
tmpdir=""
if args.tmpdir is not None:
    tmpdir= "--tmpdir " + args.tmpdir
refRegions = None
if args.refRegions is not None:
    refRegions = open(args.refRegions,'w')
    
recallCommands = []
indelFiles=[]

i=0
gaps = [l.split() for l in gapFile]
nGaps = 0
for g in gaps:
    if len(g)> 0 and g[0][0] != "#":
        nGaps+=1
if nGaps == 0:        
    sys.exit(0)
prevEnd = 0
asmFai  = None

if args.asm is not None:
    asmFaiFile = open(args.asm + ".fai")
    asmFai = {l.split()[0]: int(l.split()[1]) for l in asmFaiFile}
    
if args.ref is not None:
    refFaiFile = open(args.ref + ".fai")
    refFai = {l.split()[0]: int(l.split()[1]) for l in refFaiFile}

def ExpandRegion(contig, start, end, fai, window):
    start = max(0, start - window)
    end   = min(fai[contig], end + window)
    return (start, end)

i=0
headerLine = None
header = None
if gaps[0][0][0] == "#":
    headerLine = "\t".join(gaps[0])
    headerPre="%".join(gaps[0])
    i+=1
    headerVals = headerLine[1:].split()
    header = {headerVals[i]: i for i in range(0,len(headerVals)) }

prevREnd = 0
prevAEnd = 0
gapGroups = []
if header is not None:
    chromI = header["chrom"]
    qNameI = header["qName"]
    qStartI = header["qStart"]
    qEndI = header["qEnd"]
    svLenI = header["svLen"]
    tEndI = header["tEnd"]
    tStartI = header["tStart"]
else:
    chromI = 0
    tStartI = 1
    tEndI = 2
    qNameI = 7
    qStartI = 8
    qEndI = 9
    svLenI = 4

headerStr = "__HEADER__".join(headerVals)
#import pdb
#pdb.set_trace()
while i <len(gaps):
    j=i
    maxGap = int(gaps[i][svLenI])

    if i == 0 or (i > 0 and (gaps[i][chromI] != gaps[i-1][chromI] or gaps[i][qNameI] != gaps[i-1][qNameI])):
        prevREnd = 0
        prevAEnd = 0
    gapGroups.append([gaps[i]])
    clusterSize = 1
    gapLines = ""
    if headerVals is not None:
        gapLines="#"+"%".join(headerVals) + ";"
    
    gapLines += "%".join(gaps[i])
    #
    # as long as the next gap ends within the window of this gap, advance to include mor gaps to recall
    #
    maxGap=max(2*int(gaps[j][4]),args.window)
    maxWindow=max(2*int(gaps[j][svLenI]),args.window)
    while j + 1 < len(gaps) and \
         int(gaps[j][tEndI]) + int(gaps[j][svLenI]) + maxWindow > int(gaps[j+1][tStartI]) and \
         gaps[j][chromI] == gaps[j+1][chromI] and \
         gaps[j][qNameI] == gaps[j+1][qNameI]:
         maxGap=max(maxGap,int(gaps[j][svLenI]))

         j+=1
         clusterSize+=1
         gapLines+=";"+"%".join(gaps[j])
         rEndWindow=max(2*int(gaps[j][svLenI]),maxWindow)         
         gapGroups[-1].append(gaps[j])

#    sys.stderr.write( str(int(gaps[j][2]) + args.window - int(gaps[i][1]) + args.window ) + "\t" + str(maxGap) + "\n")
#    sys.stderr.write("Validating " + str(j-i+1) + " sites.\n")
#    sys.stderr.write("\n\nstarting at " + str(i) + "\n")
#    sys.stderr.write("\n".join(gapLines.split(";")))
    
    (rStartExp, rEndExp) = ExpandRegion(gaps[i][chromI], int(gaps[i][tStartI]), int(gaps[j][tEndI]), refFai, maxWindow)
 
    if asmFai is not None:
        (aStartExp, aEndExp) = ExpandRegion(gaps[i][qNameI], int(gaps[i][qStartI]), int(gaps[j][qEndI]), asmFai, maxWindow)
    else:
        aStartExp = int(gaps[i][qStartI]) - window
        aEndExp   = int(gaps[i][qEndI]) + window

    
    rStartExp = max(rStartExp, prevREnd)
    aStartExp = max(aStartExp, prevAEnd)

    #
    # Check to see if this region has already been recalled.
    #
    if rStartExp > rEndExp or aStartExp > aEndExp:
        i=j+1
        continue
    rRegion="{}:{}-{}".format(gaps[i][chromI], rStartExp, rEndExp)
    aRegion="{}:{}-{}".format(gaps[i][qNameI], aStartExp, aEndExp)
#    print rRegion + "\t"+ aRegion
    if refRegions is not None:
        refRegions.write("{}\t{}\t{}\t{}\t{}\n".format(gaps[i][chromI], rStartExp, rEndExp, i, j))

    keep = ""
    if args.keep:
        keep = " --keep "

    ngmlr=""
    if args.ngmlr is not None:
        ngmlr = " --ngmlr {} ".format(args.ngmlr)
    else:
        ngmlr = " --ngmlr 10000"

    indelsOption = ""
    if args.indels is not None:
        indelFile = "{}/{}.{}-{}.indels.bed".format(args.indelDir, gaps[i][chromI], rStartExp, rEndExp)
        indelsOption = " --indels {} ".format(indelFile) 
        indelFiles.append(indelFile)

    asmOption = ""
    if args.asm is not None:
        asmOption = " --asm " + args.asm
    
    refLen = rEndExp - rStartExp
    
    scriptDir=os.path.dirname(os.path.realpath(sys.argv[0]))
    command = (scriptDir+'/RecallRegion.py  {} --ref {} --asmRegion \"{}\" --window 0 --refRegion {} --index {} --ngaps {} {} --minAlignLength {} {} {} --header {} --blasr {} ##SEP##{}'.format(asmOption, args.ref, aRegion, rRegion, i, clusterSize, keep, int(0.75*refLen) , ngmlr, indelsOption, headerStr, args.blasr, gapLines))

    recallCommands.append(command)

    i=j+1
    prevREnd = rEndExp
    prevAEnd = aEndExp
#    sys.stderr.write("recall_command " + command + "\n")
    
if args.commands is not None:
    cf = open(args.commands, 'w')
    cf.write("\n".join(recallCommands) + "\n")
    cf.close()
    exit(0)
if args.split is not None:
    #split gap groups into N separate files.
    if len(gapGroups) > args.split:
        gapsPerFile = len(gapGroups)/ args.split
    else:
        gapsPerFile = len(gapGroups)
        args.split=1


    start=0
    end=0
    idx=0

    if gapsPerFile > 0:
        for i in range(0,args.split-1):
            end = start + gapsPerFile
            outFile = open(args.splitDir + "/gaps.bed.{}".format(idx),'w')
            if headerLine is not None:
                outFile.write(headerLine + "\n")
                outFile.write("\n".join(["\t".join(gap) for gapGroup in gapGroups[start:end] for gap in gapGroup])+"\n")
                outFile.close()
            idx+=1
            start = end
    outFile = open(args.splitDir + "/gaps.bed.{}".format(idx),'w')
    end=len(gapGroups)
    if headerLine is not None:
        outFile.write(headerLine + "\n")
    outFile.write("\n".join(["\t".join(gap) for gapGroup in gapGroups[start:end] for gap in gapGroup])+"\n")
    outFile.close()
    sys.exit(0)
    
        
    
def Run(command):

    commandTuple = command.split("##SEP##")
    gapLinesFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".txt", delete=False, mode='w')
    gapLinesFile.write(commandTuple[1])
    gapLinesFile.close()
    commandTuple[0] += " --gapLines " + gapLinesFile.name
#    sys.stderr.write(commandTuple[0] + "\n")    
    res = subprocess.check_output(commandTuple[0], shell=True)
    if args.keep is False:
        cleanup = "rm -f {}".format(gapLinesFile.name)
        subprocess.call(cleanup.split())
    return res.rstrip()

alns=[]
import codecs

if args.nproc > 1:
    pool = Pool(args.nproc)

    aln = pool.map(Run, recallCommands)
    pool.close()

    for a in aln:
        if len(a) > 0 and a != "\n":
            alns.append(codecs.latin_1_decode(a)[0])
else:
    for c in recallCommands:
        aln = Run(c)
        if aln != "\n":
            alns += codecs.latin_1_decode(aln)[0].split("\n")
            
oneHeader = []
first=True
outFile.write(headerLine+"\trecMethod\n")
#import pdb
for a in range(0, len(alns)):
#    alnRows = alns[a].split('\n')
    aln = alns[a]
    if first and len(aln) > 0 and len(aln) > 0 and aln[0] == "#":
#        pdb.set_trace()
        oneHeader.append(aln)
        first = False
    else:
        if len(aln) > 0 and aln[0] != '#':
#            pdb.set_trace()
            oneHeader.append(aln)

#import pdb
#pdb.set_trace()
outFile.write("\n".join(oneHeader)+"\n")
    

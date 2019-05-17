#!/usr/bin/env python

import argparse
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
import pysam
import tempfile
import subprocess
import time
import math
from multiprocessing import Process, Lock, Semaphore, Pool

import re
import os
import sys
import codecs

ap = argparse.ArgumentParser(description="Given a gaps.bed file, splice variants into the reference, remap reads to this, and count gaps. ")
ap.add_argument("--gaps",   help="Gaps file (can be insertions and deletions).", required=True)
ap.add_argument("--ref",    help="Reference file.", required=True)
ap.add_argument("--paIndex",  help="When other jobs are running, which one is this", default=0,type=int),
ap.add_argument("--paNumber",  help="How many other jobs are running", default=1,type=int),
ap.add_argument("--reads",  help="FOFN of aligned reads.", required=False)
ap.add_argument("--blasr",  help="which blasr to use", required=True)
ap.add_argument("--flank",  help="Use this amount of the reference when aligning reads", default=900,type=int)
ap.add_argument("--window", help="Collect all gaps within this window.", type=int, default=1000)
ap.add_argument("--minGap", help="Minimum gap size.", default=30, type=int)
ap.add_argument("--maxSize", help="Maximum size of region.", default=None, type=int)
ap.add_argument("--maxCoverage", help="Maximum coverage to allow.", default=100, type=int)
ap.add_argument("--maxInsertion", help="Maximum length database to generate", default=5000,type=int)
ap.add_argument("--out",    help="Output file.", default="/dev/stdout")
ap.add_argument("--tmpdir", help="Temporary directory",default=None)
ap.add_argument("--nproc",  help="Number of threads.", default=4,type=int)
ap.add_argument("--keep",   help="Keep temporary files", default=False, action="store_true")
ap.add_argument("--split",  help="Split input gaps.", default=None, type=int)
ap.add_argument("--splitDir",  help="Split input gaps into this dir.", default=None)
ap.add_argument("--mssm",   help="Parse mssm lines.", default=False, action="store_true")
ap.add_argument("--falcon", help="Parse in falcon format", default=False, action="store_true")
ap.add_argument("--count",  help="Just count clusters, don't do any validation.", default=None)
ap.add_argument("--genotypeVcf", help="Use SNVs to call genotype.", default=None)
ap.add_argument("--sample", help="Sample to select from genotype", default=None)
ap.add_argument("--maxMergedSV", help="Maximum merged SV", default=10000,type=int)
ap.add_argument("--separate-haplotypes", help="Separate haplotypes", dest="separateHaplotypes", default=False, action="store_true")

args = ap.parse_args()
if args.count is None:
    outFile =open(args.out,'w')

if args.tmpdir is None:

    if "TMPDIR" not in os.environ or os.environ["TMPDIR"] == "":
        print("ERROR. The TEMPDIR variable must be set or --tmpdir specified on as a command  argument\n")
        sys.exit(1)
    else:
        args.tmpdir = os.environ["TMPDIR"]

def WaitOnFile(fn):
    cmd=["lsof",fn]
    while True:
        res=subprocess.Popen(cmd, stdout=subprocess.PIPE)
        l=codecs.latin_1_decode(res.stdout.read())[0]
        res.wait()
        if l == '':
            return
        sys.stderr.write("Waiting on temporary file.\n")
        
def SafeFetch(f, c, s, e):
    s1 = f.fetch(c,s,e)
    s2 = f.fetch(c,s,e)
    while s1 != s2:
        s1 = f.fetch(c,s,e)
        s2 = f.fetch(c,s,e)
    return s1
        
    
def CountRefCoverage(allLines, checkHaplotype=False):
#    samFile = pysam.AlignmentFile(samFileName,'r')
    refs = {}
    refs["re"] = 0
    refs["db"] = 0
    lines=allLines.split("\n")
    if checkHaplotype == True:
        haps = {}
        haps["re"]={}
        haps["db"]={}        
        haps["re"][0] = 0
        haps["re"][1] = 0
        haps["re"]['u'] = 0        
        haps["db"][0] = 0
        haps["db"][1] = 0
        haps["db"]['u'] = 0        

    for line in lines:
        if len(line) > 0 and line[0] == '@':
            continue
        vals = line.split()
        if len(vals) < 2:
            continue
        rn = vals[2]

        if rn not in refs:
            refs[rn] = 0
        refs[rn]+=1
        if checkHaplotype:
            if read.query_name[-2:] == "/0":
                hap=0
            elif read.query_name[-2:] == "/1":
                hap=1
            elif read.query_name[-2:] == "/u":
                hap='u'
            if rn in haps:
                haps[rn][hap]+=1

    if checkHaplotype == False:
        return refs
    else:
        return haps

def CountGaps(samFileName, rStart, rEnd):
    samFile = pysam.AlignmentFile(samFileName,'r')
    nGaps = 0
    nReads = 0
    nGapBases=0
    for read in samFile.fetch():
        nReads+=1
        
        if read.cigartuples is not None:
            refPos = 0
            for t in read.cigartuples:
                
                if (rStart  <= refPos and refPos < rEnd):
                    if (t[0] == 1 or t[0] == 2) and t[1] > args.minGap:
                        nGaps+=1
                        nGapBases+=t[1]
                    if t[0] == 4:
                        nGaps+=1

                if t[0] == 0 or t[0] == 2 or t[0] == 7 or t[0] == 8:
                    refPos+=t[1]
    return nGaps, nGapBases, nReads


def CountCoverage(samFileName, rStart, rEnd):
    samFile = pysam.AlignmentFile(samFileName,'r')
    nGaps = 0
    nReads = 0
    nCovered = 0
    for read in samFile.fetch():
        nReads+=1
        covered=False
        hasGap =False
        
        if read.cigartuples is not None:
            refPos = 0
            for t in read.cigartuples:
                
                if (rStart  <= refPos and refPos < rEnd):
                    if (t[0] == 1 or t[0] == 2) and t[1] > args.minGap:
                        hasGap=True
                        break
                    if t[0] == 4:
                        hasGap=True
                        break

                if t[0] == 0 or t[0] == 2 or t[0] == 7 or t[0] == 8:
                    refPos+=t[1]
        if hasGap == False:
            nCovered+=1
    return nCovered, nReads


class Gap:
    def __init__(self):
        self.chrom=""
        self.start= 0
        self.end  = 0
        self.svType=0
        self.svLen=0
        self.hap=None
        self.rem=[]
        self.line=""

def ParseFalconLine(vals):
    gap=Gap()
    gap.chrom=vals[0]
    gap.start=int(vals[1])
    gap.end=int(vals[2])
    gap.svType="deletion"
    return gap

def ParseMSSMLine(vals):

    if len(vals) > 0:
        gap=Gap()
        gap.chrom=vals[0]
        gap.start=int(vals[1])
        gap.end=int(vals[2])
        gap.svType=vals[3]
        gap.svLen=len(vals[5])
        gap.line = "\t".join(vals)
        return gap
    else:
        return None

def ParseGapsLine(vals):
    if len(header) >= header["svSeq"]:
        chromI    = header["chrom"]
        tStartI   = header["tStart"]
        tEndI     = header["tEnd"]
        svTypeI   = header["svType"]
        svLenI    = header["svLen"]
        svSeqI    = header["svSeq"]
        remainder = 6
    else:
        chromI, tStartI, tEndI, svTypeI, svLenI, svSeq = 0,1,2,3,4,5

        
    if len(vals) >= header["svSeq"] :
        gap=Gap()
        gap.chrom=vals[chromI]
        gap.start=int(vals[tStartI])
        gap.end=int(vals[tEndI])
        gap.svType=vals[svTypeI]
        gap.seq=vals[svSeqI]
        gap.line="\t".join(vals)
        gap.svLen=int(vals[svLenI])
        if "hap" in header:
            gap.hap=vals[header["hap"]]
        else:
            gap.hap=None
        return gap
    else:
        return None


def GetEnd(v):
    if args.mssm:
        return v.start + v.svLen
    elif args.falcon:
        return v.end
    else:
        if v.svType == "deletion":
            return v.end
        else:
            return v.start+1

def GetLength(v):
    return v.svLen

def GetStart(v):
    return v.start

def GetChrom(v):
    return v.chrom

gapsFile = open(args.gaps)
headerLine = gapsFile.readline()
gapLines = gapsFile.readlines()

if len(headerLine) > 0 and headerLine[0] != '#':
    gapLines = [headerLine] + gapLines
    header=None
else:
    v = headerLine[1:].split()
    header = {v[i]: i for i in range(0,len(v))}

if args.mssm:
    gaps = [ParseMSSMLine(line.split()) for line in gapLines]
elif args.falcon:
    gaps = [ParseFalconLine(line.split()) for line in gapLines]
else:
    gaps = [ParseGapsLine(line.split()) for line in gapLines]

if args.paNumber > 1:
    sys.stderr.write("ngaps: " + str(len(gaps)) + "\n")
    gapsPerJob=int(math.ceil(float(len(gaps))/args.paNumber))
    start=args.paIndex*gapsPerJob
    end=min(len(gaps), (args.paIndex+1)*gapsPerJob)
    sys.stderr.write("Processing " + str(start) + " - " + str(end) + "\n")
    gaps=gaps[start:end]
    
temp=[]
for g in gaps:
    if g is not None:
        temp.append(g)
gaps= temp

svClusters = []
i=0

def CompareHaplotypes(a, b, care=False):
    if care == False:
        return True
    else:
        return a.hap == b.hap



while i < len(gaps):
    svClusters.append([])
    svClusters[-1].append(gaps[i])

    while i < len(gaps)-1 and \
        GetChrom(gaps[i]) == GetChrom(gaps[i+1]) and \
        GetStart(gaps[i+1])-GetEnd(gaps[i]) < args.window and \
        GetLength(gaps[i+1]) < args.maxMergedSV and \
        CompareHaplotypes(gaps[i], gaps[i+1], args.separateHaplotypes):
        sys.stderr.write("{} \t {} \t{}\t{}\n".format(i, GetStart(gaps[i+1]), GetEnd(gaps[i]), GetStart(gaps[i+1])-GetEnd(gaps[i])))
        i+=1
        svClusters[-1].append(gaps[i])
    i+=1



if args.split is not None:
    sStart = 0
    sEnd = 0
    b = len(svClusters)/args.split

    for i in range(0,args.split-1):
        sEnd+=b
        clustFile = open("{}/gaps.bed.{}".format(args.splitDir, i),'w')
        clustFile.write(headerLine)
        clustFile.write("\n".join([c.line for cl in svClusters[sStart:sEnd] for  c in cl ]))        
        clustFile.close()
        sStart = sEnd
    i=args.split-1
    clustFile = open("{}/gaps.bed.{}".format(args.splitDir, i),'w')
    clustFile.write(headerLine)
    clustFile.write("\n".join([c.line for cl in svClusters[sStart:] for  c in cl ]))
    clustFile.close()
    sys.exit(0)
        
ref = pysam.FastaFile(args.ref)

refFaiFile = open(args.ref + ".fai")
refFai = {l.split()[0]: int(l.split()[1]) for l in refFaiFile}
readsFofn = open(args.reads)
bamFiles = [pysam.AlignmentFile(line.rstrip(), 'rb') for line in readsFofn]

if args.count is not None:
    countFile = open(args.count,'w')
    countFile.write("\n".join([str(len(c)) for c in svClusters]))
    countFile.close()
    sys.exit(0)
    
def WriteSeq(fh,seq,seqName):
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq),id=seqName,name="",description=""),fh,"fasta")

sem = Semaphore(1)

def SpliceTestLine(svs):
    svIndex = 0
    for sv in svs:
        if sv.chrom not in refFai:
            continue
        if svIndex == 0:
            prefixStart = max(0, sv.start-args.flank)
            prefixEnd   = sv.start
            sem.acquire()
            gapPrefix   = SafeFetch(ref, sv.chrom, prefixStart, prefixEnd)
            sem.release()
            spliceSeqs  = [gapPrefix]
            refPos      = sv.start
        
        if args.mssm:
            spliceSeqs.append(sv[5])

            refPos = sv[2]

            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1][1]:
                    continue
                sem.acquire()
                between = SafeFetch(ref, sv[0], refPos, svs[svIndex+1][1])
                sem.release()
                spliceSeqs.append(between)
        elif args.falcon:
            # only checking deletions
            refPos = sv[2]
            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1][1]:
                    continue
            
        else:
            if sv.svType == "insertion":
                svSeq = sv.seq
                if len(svSeq) > args.maxInsertion:
                    mid=args.maxInsertion/2
                    svSeq = sv.seq[0:mid] +sv.seq[-mid:]
                    sys.stderr.write("keeping center of insertion {}\n".format(len(sv.seq)))
                spliceSeqs.append(svSeq)
                refPos = sv.start
            else:
                refPos = sv.end
        

            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1].start:
                    sys.stderr.write("Ignoring an sv {} {} {}\n".format(svs[svIndex+1].chrom, svs[svIndex+1].start, svs[svIndex+1].end))
                    svIndex+=1
                    continue
                else:
                    sem.acquire()
                    between = SafeFetch(ref, sv.chrom, refPos, svs[svIndex+1].start)
                    sem.release()
                    refPos += len(between)
                    spliceSeqs.append(between)
        svIndex +=1

    svIndex = len(svs)-1
    sem.acquire()
    suffix = SafeFetch(ref, svs[svIndex].chrom, refPos, min(refFai[svs[svIndex].chrom], refPos + args.flank))
    sem.release()
    
    spliceSeqs.append(suffix)
    dbSeq= ''.join(spliceSeqs)
    refChrom = svs[0].chrom


    refStart = max(0, svs[0].start - args.flank)
    refEnd   = min(refFai[svs[svIndex].chrom], GetEnd(svs[svIndex]) + args.flank)
    nBases = 0
    if args.maxSize is not None and refEnd - refStart > args.maxSize:
        if args.genotypeVcf is None:
            results="\n".join(["{}:{}-{}\t{}\t{}".format(sv.chrom,sv.start,sv.end,0,0) for sv in svs])
            return results
        
        
    sem.acquire()
    refSeq = SafeFetch(ref,refChrom, refStart, refEnd)
    sem.release()
    tempFileNames = []
    fSuffix="."+str(refPos) + ".fasta"
    sSuffix="."+str(refPos) + ".sam"    
    rFile     = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=fSuffix, delete=False, mode='w')
    tempFileNames.append(rFile.name)
    dbFile    = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=fSuffix, delete=False, mode='w')
    tempFileNames.append(dbFile.name)
    readsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=fSuffix, delete=False, mode='w')
    tempFileNames.append(readsFile.name)

    WriteSeq(dbFile, dbSeq, "db")
    WriteSeq(dbFile, refSeq, "re")
    WriteSeq(rFile, refSeq, "re")    
    #
    # Now collect all of the sequences.
    #
    fetchStart = svs[0].start - args.flank
    fetchEnd   = GetEnd(svs[-1]) + args.flank
    # just count one breakpoint if large event.
    if fetchEnd - fetchStart > 30000:
        sys.stderr.write("******Truncating fetch region {}\n".format(fetchEnd-fetchStart))
        fetchEnd = svs[0].start  + args.flank
    sys.stdout.write("Fetching from region " + str(fetchEnd-fetchStart) + " " + str(svs[-1].svType) + "\n")
    
    
    dbFile.close()
    rFile.close()
    nBases = 0
    if args.genotypeVcf is not None:
        #
        # This uses the SNV vcf in the argument to partition reads, and genotype by phase tag.
        #
        print("about to start genotyping\n")
        dipSamFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".dip"+sSuffix, delete=False, mode='w')
        tempFileNames.append(dipSamFile.name)
        dipSamFile.close()
        dipHandle = pysam.AlignmentFile(dipSamFile.name, 'wh', header=bamFiles[0].header)
        
        sem.acquire()        
        for b in range(0,len(bamFiles)):
            for read in bamFiles[b].fetch(sv.chrom, fetchStart, fetchEnd+1):
                dipHandle.write(read)
        sem.release()
        #
        # Now partition the file by haplotype
        #
        hap0SamFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".hap0"+sSuffix, delete=False, mode='w')
        hap1SamFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".hap1"+sSuffix, delete=False, mode='w')
        unassignedSamFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".unassigned"+sSuffix, delete=False, mode='w')        
        regionVCF = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".vars.vcf", delete=False, mode='w')
        tempFileNames += [hap0SamFile.name, hap1SamFile.name, regionVCF.name, unassignedSamFile.name]
        
        vcfStart = max(0,svs[0].start - args.flank*10)
        vcfEnd   = GetEnd(svs[-1]) + args.flank*10

        tabixCommand = "tabix -h {} {}:{}-{}".format(args.genotypeVcf, sv.chrom, vcfStart, vcfEnd)
        subprocess.call(tabixCommand.split(), stdout=regionVCF)
        regionVCF.close()        
        partitionCommand = "{}/partitionByPhasedSNVs --vcf {} --sam {} --rgn {}:{}-{} --pad 10000 --h1 {} --h2 {} --ref {} --minGenotyped 1 --sample {} --unassigned {}".format("/net/eichler/vol5/home/mchaisso/projects/pbgreedyphase", regionVCF.name, dipSamFile.name, sv.chrom, fetchStart, fetchEnd, hap0SamFile.name, hap1SamFile.name, args.ref, args.sample, unassignedSamFile.name )
        subprocess.call(partitionCommand.split())

        sams = [hap0SamFile.name, hap1SamFile.name, unassignedSamFile.name]
        haps = ["0", "1", "u"]
        sem.acquire()
        for i in range(0,3):

            samHandle = pysam.AlignmentFile(sams[i], 'r')
            for read in samHandle.fetch():
                nBases+=min(read.reference_end, fetchEnd) - max(fetchStart,read.reference_start)                    
                WriteSeq(readsFile, read.seq, read.query_name + "/" + haps[i])
        sem.release()
    else:
        sem.acquire()
        for b in range(0,len(bamFiles)):
            for read in bamFiles[b].fetch(sv.chrom, fetchStart, fetchEnd+1):
                nBases+=min(read.reference_end, fetchEnd) - max(fetchStart,read.reference_start)
                WriteSeq(readsFile, read.seq, read.query_name)
        sem.release()
    readsFile.close()
    
        
#    rsFile  = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".sam", delete=False, mode='w')
    dbsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=sSuffix, delete=False, mode='w')
    
    commandOptions = " -maxMatch 25 -sdpMaxAnchorsPerPosition 5 -sdpTupleSize 10 -sam -bestn 1 -affineOpen 5 -affineExtend 5 -nproc 8 -out /dev/stdout -minAlignLength {} ".format(int(1.5*args.flank))

    dbCommand = "{} {} {} -preserveReadTitle -clipping soft ".format(args.blasr, readsFile.name, dbFile.name, dbsFile.name) + commandOptions
    tempFileNames.append(dbsFile.name)
    

    dn = open(os.devnull)
    regionLength = fetchEnd - fetchStart
    coverage = nBases/regionLength
    sys.stderr.write("coverage: " + str(coverage)+ "\n")

    if args.genotypeVcf is not None:
        genotype=True
    else:
        genotype=False
    fs=0
    rs=0
    if coverage < args.maxCoverage:
        proc=subprocess.Popen(dbCommand.split(),stderr=dn,stdout=subprocess.PIPE)
        alnLines=codecs.latin_1_decode(proc.stdout.read())[0]
        proc.wait()
        dbsFile.close()
#        WaitOnFile(dbsFile.name)

        fs=os.path.getsize(dbsFile.name)
        rs=os.path.getsize(readsFile.name)
        cov = CountRefCoverage(alnLines, genotype)
    else:
        sys.stderr.write("Skipping event from coverage " + str(coverage) + "\n")
        cov = { "db": 0, "re": 0}
    sys.stderr.write(svs[0].svType+ " " + str(cov) + "\n")
    
    if args.genotypeVcf is False and "db" in cov and cov["db"] < 500:
        print( "spliced " + str(len(svs))        )
        print( "cov:")
        print( cov)
        print( "{}:{}-{}".format(svs[0].chrom,svs[0].start,svs[0].end))


    cleanup = "/bin/rm " + " ".join(tempFileNames)

    if args.keep is False:
        subprocess.call(cleanup.split())
    else:
        print(cleanup)

    if args.genotypeVcf is None:
        dbCov = cov["db"]
        rCov  = cov["re"]
    
        results="\n".join(["{}:{}-{}\t{}\t{}".format(sv.chrom,sv.start,sv.end,dbCov,rCov ) for sv in svs])
    else:
        results="\n".join(["{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}".format(sv.chrom,sv.start,sv.end,\
                                                                                     cov["db"][0],cov["db"][1],cov["db"]['u'],\
                                                                                     cov["re"][0],cov["re"][1],cov["re"]['u']\
                                                                     ) for sv in svs])
        sys.stdout.write(results + "\n")
    return results
if header is not None:
    if args.genotypeVcf is None:
        outFile.write("#region\tnAlt\tnRef\n")
    else:
        outFile.write("#region\tnAltH0\tnAltH1\tnAltUn\tnRefH0\tnRefH1\tnRefUn\n")
if args.nproc > 1:    
    pool = Pool(args.nproc)

    res = pool.map(SpliceTestLine, svClusters)
    pool.close()
    outFile.write("\n".join(res)+"\n")
else:
    for svCluster in svClusters:
        res = SpliceTestLine(svCluster)
        outFile.write(res+ "\n")

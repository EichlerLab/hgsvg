#!/usr/bin/env python

import tempfile
import subprocess
import argparse
import pysam
import re
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import os
import sys

ap = argparse.ArgumentParser(description="")
ap.add_argument("--a", help="First file")
ap.add_argument("--b", help="Second file")
ap.add_argument("--asm", help="Assemblies")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--tmpdir", help="Use this directory for temp.", default=".")
ap.add_argument("--minOvp", help="Minimum overlap to store.", default=10000,type=int)
ap.add_argument("--keep", help="Keep output file.",default=False,action='store_true')
ap.add_argument("--center",help="Attempt to find the central region of this length.", type=int,default=20000)
ap.add_argument("--blasr",  help="Specify alternative blasr", default="/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/blasr")
args = ap.parse_args()


def GetSeq(name, ref, tempFile):
    sys.stderr.write("Grabbing " + ref + "'" + name + "'\n")

    faidx = "samtools faidx {} {}".format(ref, name)
    sys.stderr.write(faidx + "\n")
    p = subprocess.Popen(faidx.split(), stdout=tempFile)
    p.wait()
    sys.stderr.write("done\n")

aTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=False)
bTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=False)
samTempFile = tempfile.NamedTemporaryFile(suffix=".sam", dir=args.tmpdir, delete=False)
GetSeq(args.a, args.asm, aTempFile)
GetSeq(args.b, args.asm, bTempFile)

seqASet = SeqIO.to_dict(SeqIO.parse(open(aTempFile.name), "fasta"))
seqB = SeqIO.read(open(bTempFile.name), "fasta")
nre=re.compile("^(N*)[^N]*(N*)$")


def GetN(seq):
    m= nre.match(str(seq.seq))
    if m is None:
        return (0,len(seq.seq))
    else:
        g =m.groups()
        return (len(g[0]), len(seq.seq) - len(g[1]))


command = "{} {} {} -bestn 1 -sam -out {} -clipping soft -maxMatch 20 -sdpMaxAnchorsPerPosition 3 -advanceExactMatches 10 -sdpTupleSize 13 -extend -maxExtendDropoff 50 -preserveReadTitle ".format(args.blasr, aTempFile.name, bTempFile.name, samTempFile.name)
#sys.stderr.write(command + "\n")
subprocess.call(command.split())


(bNPre, bNPost) = GetN(seqB)

#
# parse the sam file to get the cooridnates of the overlap
#
fileSize = os.path.getsize(samTempFile.name)
#sys.stderr.write("Results" + aTempFile.name + " " + bTempFile.name + " " + samTempFile.name + " " + str(fileSize) + "\n")

if fileSize > 0:
    aln = pysam.AlignmentFile(samTempFile.name, "r")


    for read in aln.fetch():

        seqA = seqASet[read.query_name]
        (aNPre, aNPost) = GetN(seqA)    
        coords = [read.qstart, read.qend, read.reference_start, read.reference_end]
        qPos = 0
                            
        if coords[0] == coords[1] and coords[0] == 0:
            continue
        else:
            if read.qend - read.qstart > args.minOvp:
                nGapBases = 0
                if read.cigartuples is not None:
                    for tup in read.cigartuples:
                        if tup[0] == 4 or tup[0] == 5:
                            qPos += tup[1]
                        elif tup[0] == 0 or tup[0] == 8 or tup[0] == 7:
                            qPos += tup[1]
                        elif tup[0] == 1 or tup[0] == 2:
                            if qPos > read.qstart and qPos < read.qend and tup[1] > 50:
                                #sys.stderr.write("read " + read.query_name + " has indel relative to " + args.b + " of len " + str(tup[1]) +"\n")
                                nGapBases += tup[1]
                            if tup[0] == 1:
                                qPos += tup[1]
                p = read.get_aligned_pairs(matches_only=True)
                if len(p) == 0:
                    continue
                l = 0
                r = len(p)-1

                while (l < r and p[r][0] - p[l][0] > args.center):
                    l+=1
                    r-=1

                print("{}\t{}\t{}\t".format(read.query_name, aNPre, aNPost) + "{}\t{}\t{}\t".format(args.b, bNPre, bNPost) +  "\t".join([str(i) for i in coords]) + "\t" + str(nGapBases) + "\t{}\t{}\t{}\t{}".format(p[l][0], p[r][0], p[l][1], p[r][1]))

if args.keep == False:
    command="rm -f {} {} {}".format(aTempFile.name, bTempFile.name, samTempFile.name)
    subprocess.call(command.split())
else:
    sys.stderr.write("kept temp files\n{}\n{}\n{}\n".format(aTempFile.name, bTempFile.name, samTempFile.name))

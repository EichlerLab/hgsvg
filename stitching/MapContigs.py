#!/usr/bin/env python

import argparse
import pysam
import tempfile
import sys
import subprocess
from multiprocessing import Process, Lock, Semaphore, Pool
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="Map contigs to the region that was patched")
ap.add_argument("--contigs", help="Input contigs file.", required=True)
ap.add_argument("--ref", help="Reference",required=True)
ap.add_argument("--tmpdir", help="tempdir", required=True)
ap.add_argument("--out", help="Output file.", required=True)
ap.add_argument("--nproc", help="Number of cores.", type=int, default=6)
ap.add_argument("--blasr",  help="Specify alternative blasr", default="/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/blasr")
args = ap.parse_args()

contigFile = open(args.contigs)
ref      = pysam.FastaFile(args.ref)
outFile = open(args.out,'w')

def Run(files):
    command = "{} {} {} -minMapQV 30 -minAlignLength 10000 -sam -out {} -clipping soft  -preserveReadTitle -alignContigs ".format(args.blasr, \
                                                                                                        files[0],\
                                                                                                        files[1],\
                                                                                                        files[2])
    print(command+"\n")
    subprocess.call(command.split())


    samFile = open(files[2])
    alnLines = []
    for line in samFile:
        if line[0] == '@':
            continue
        vals = line.split()
        vals[2] = files[3]
        vals[3] = str(int(vals[3]) + files[4])
        alnLines.append("\t".join(vals))

    command = "rm -f {} {} {}".format(files[0], files[1], files[2])
    subprocess.call(command.split())

    return '\n'.join(alnLines)

fileList = []

for contig in SeqIO.parse(contigFile, "fasta"):
    queryTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=False, mode='w')
    sys.stderr.write("writing " + str(contig) + "\n")
    SeqIO.write(contig, queryTempFile, "fasta")
    queryTempFile.close()
    
    vals = contig.id.split(".")
    chrom = vals[-3]
    refStart = int(vals[-2])
    refEnd   = int(vals[-1])
    
    seq = ref.fetch(reference=chrom, start=refStart, end=refEnd)
    refTempFile = tempfile.NamedTemporaryFile(suffix=".fasta", dir=args.tmpdir, delete=False, mode='w')
    refTempFile.write(">ref\n")
    refTempFile.write(seq+"\n")
    refTempFile.close()

    
    samTempFile = tempfile.NamedTemporaryFile(suffix=".sam", dir=args.tmpdir, delete=False)
    samTempFile.close()
    
    fileList.append([queryTempFile.name, refTempFile.name, samTempFile.name, chrom, refStart,refEnd])
    

pool = Pool(args.nproc)
ovp = pool.map(Run, fileList)
pool.close()
outFile.write("\n".join(ovp)+"\n")
    


#!/usr/bin/env python

import argparse
from multiprocessing import Process, Lock, Semaphore, Pool
import sys
import subprocess
import re

from functools import cmp_to_key  # cmp_to_key Converts a python2 cmp comparator func. to a python3 key func.


ap = argparse.ArgumentParser(description="Given an alignment bed file, run overlaps between contigs that are near each other")
ap.add_argument("bed", help="Alignment bed file.")
ap.add_argument("asm", help="Assemblies file")
ap.add_argument("--chrom", help="Process this chromosome only.", default =None)
ap.add_argument("--path", help="Path to HGSVG stitching scripts.")
ap.add_argument("--ahead", help="Allow overlap between alignments ending this far ahead of current.",\
                type=int, default=20000)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of processors.", default=1,type=int)
ap.add_argument("--tmpdir", help="Pass along tempdir.", default=".")
ap.add_argument("--blasr", help="Blasr command to run.", default=None)
ap.add_argument("--minOvp", help="Minimum overlap to store", default=10000,type=int)
ap.add_argument("--start", help="First line to start overlapping", default=0,type=int)
ap.add_argument("--n", help="N entries to compute overlap for", default=-1,type=int)
args = ap.parse_args()

def ParseBedLine(line):
    v = line.split()
    return [v[0], int(v[1]), int(v[2]), v[3]]

bedFile = open(args.bed)
outFile = open(args.out,'w')
alignments= [ParseBedLine(line) for line in bedFile]


if (args.n == -1):
    last = len(alignments)
else:
    last = min(len(alignments), args.start+args.n)
start=args.start
stepBack=start

while stepBack > 0 and \
      alignments[stepBack-1][0] == alignments[start][0]  and \
      alignments[start][1] <= alignments[stepBack][2]+args.ahead:
    stepBack-=1

start=stepBack

alignments = alignments[start:last]

commands = []
def Run(command):
    res = subprocess.check_output(command, shell=True)
    
    return res.rstrip()

if args.blasr is not None:
    blasrCommand = " --blasr " + args.blasr
else:
    blasrCommand = ""


n=0
firstTargetOffset=args.start-start
for i in range(1+firstTargetOffset,len(alignments)):
    if args.chrom is not None and alignments[i][0] != args.chrom:
        continue
    j = i-1
    while j > 0 < len(alignments) and \
          alignments[j][0] == alignments[i][0] and \
          alignments[j][2] + args.ahead > alignments[i][1]:
        j-=1

    aAln = " ".join([alignments[k][3] for k in range(j,i)])
    commands.append( "{}/OverlapTwoContigs.py --a \"{}\" --b {} --asm {} --tmpdir {} {} --minOvp {}".format(args.path, aAln, alignments[i][3], args.asm, args.tmpdir, blasrCommand, args.minOvp))

    
    n+=1

if len(commands) == 0:
    outFile.close()
    sys.exit(0)

#
# Now process alignments
#

# Setup threads
pool = Pool(args.nproc)
ovps = []
filt = []

# Error handling overlap threads
async_ex = None

def thread_error(ex):
    global async_ex
    global pool

    if async_ex is None:
        async_ex = ex

    pool.terminate()


# Run
res = pool.map_async(Run, commands, callback=ovps.append, error_callback=thread_error)

res.wait((last-start)*1000)

if async_ex is not None:
    raise async_ex


# Finish overlaps
def GetStarts(l):
    """
    Get start position from record name.
    """
    vals = l.split()

    #print('GetStarts({})'.format(l))

    v1 = re.split('[:/]', vals[0])[1].split("-")[0]
    v2 = re.split('[:/]', vals[3])[1].split("-")[0]

    return int(v1), int(v2)


def Compare(a,b):
    sa=GetStarts(a)
    sb=GetStarts(b)
    if (sa[0] != sb[0]):
        return sa[0] < sb[0]
    else:
        return sa[1] < sb[1]


if len(ovps) > 0:
    for ovp in ovps[0]:
        ovp = ovp.decode().strip()

        if ovp:
            filt.append(ovp)
    
    filt.sort(key=cmp_to_key(Compare))

    outFile.write("\n".join(filt)+"\n")
else:
    sys.stderr.write("No overlaps were found\n")
    sys.stderr.write("Input bed file had {} - {}: {} alignments".format(start, last, last-start))

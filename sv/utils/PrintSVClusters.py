#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Given a gaps.bed file, splice variants into the reference, remap reads to this, and count gaps. ")
ap.add_argument("--gaps", help="Gaps file (can be insertions and deletions).", required=True)
ap.add_argument("--window", help="Collect all gaps within this window.", type=int, default=100)
ap.add_argument("--mssm", help="Parse mssm lines.", default=False, action='store_true')
ap.add_argument("--falcon", help="Parse in falcon format", default=False, action='store_true')
ap.add_argument("--clusters", help="Clusters file.", default=None)
ap.add_argument("--minCalls", help="Write out calls if they have at least this many events.", type=int, default=0)
ap.add_argument("--out", help="Table of cluster index.", default="/dev/stdout")
args = ap.parse_args()
clusterOut = None
if args.clusters is not None:
    clusterOut = open(args.out,'w')

clusterTable = open(args.out,'w')

class SVCluster:
    def __init__(self, chrom, start, end, svType, line):
        self.chrom = chrom
        self.start = start
        self.end   = end
        self.svtype = svType
        self.line = line
        self.svLen = 0
        self.seq = ""

        self.additional = []
        
    
def ParseFalconLine(vals):
    sv = SVCluster(vals[0], int(vals[1]), int(vals[2]), "deletion", "\t".join(vals))
    sv.svLen = int(vals[2]) - int(vals[1])
    return sv


def ParseMSSMLine(vals):
    if len(vals) > 0:
        sv = SVCluter(vals[0], int(vals[1]), int(vals[2]), vals[3], "\t".join(vals))
        sv.seq = vals[5]
        sv.additional = vals[7:]
        return sv
    else:
        return None

def ParseGapsLine(vals):
    if len(vals) > 0:
        sv = SVCluster(vals[header["chrom"]],int(vals[header["tStart"]]), int(vals[header["tEnd"]]), vals[header["svType"]], "\t".join(vals))
        sv.svLen = int(vals[header["svLen"]])
        sv.additional=vals[5:]
        return sv
    else:
        return None


def GetEnd(sv):
    if args.mssm:
        return sv.start + len(sv.seq)
    elif args.falcon:
        return sv.end
    else:
        if sv.svType == "deletion":
            return sv.end
        else:
            return sv.start+1

def GetLength(sv):
    return sv.svLen


gapsFile = open(args.gaps)
headerLine = gapsFile.readline()[1:].split()
header = {headerLine[i]: i for i in range(0,len(headerLine))}

gapLines = gapsFile.readlines()

if args.mssm:
    gaps = [ParseMSSMLine(line.split()) for line in gapLines]
elif args.falcon:
    gaps = [ParseFalconLine(line.split()) for line in gapLines]
else:
    gaps = [ParseGapsLine(line.split()) for line in gapLines]

svClusters = []
i=0
inCluster  = []
clusterCount = {}
clusterIndex = 0
while i < len(gaps):
    svClusters.append([])
    svClusters[-1].append(gaps[i])
    inCluster.append(clusterIndex)
    clusterCount[clusterIndex] = 1
    while i < len(gaps)-1 and \
          gaps[i].chrom == gaps[i+1].chrom and \
          gaps[i+1].start-gaps[i].end < args.window:
          inCluster.append(clusterIndex)
          clusterCount[clusterIndex]+=1
          i+=1
    svClusters[-1].append(gaps[i])
    clusterIndex+=1
    i+=1

clusterTable.write("clusterCount\tclusterIndex\n")
for i in range(0,len(inCluster)):
    clusterTable.write(str(clusterCount[inCluster[i]]) + "\t" + str(inCluster[i]) +"\n")
clusterTable.close()


def GetRegion(cluster):
    return "{}:{}-{}/{}".format(cluster.chrom, cluster.start, cluster.start + cluster.end,cluster.svtype)

if clusterOut is not None:
    for clust in svClusters:
        if len(clust) >= args.minCalls:
            start = min([c.start for c in clust])
            end = max([c.end for c in clust])
            clusterOut.write(clust[0].chrom + "\t{}\t{}\t{}\t".format(start,end,len(clust)) + ";".join([GetRegion(c) for c in clust]) + "\n")
    clusterOut.close()

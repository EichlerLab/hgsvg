#!/usr/bin/env python



import argparse

ap = argparse.ArgumentParser(description="Summarize unions, possibly just tally missing from Illumina")
ap.add_argument("--miss", help=" just print missing from IL.", action="store_true", default=False)

args = ap.parse_args()
import sys
headerLine = sys.stdin.readline()
hv = headerLine[1:].split()
h  = { hv[i]: i for i in range(0,len(hv)) }

pbv = []
pbu = []
ilv = []
ilu = []
bnv = []
nStats = 1
for line in sys.stdin:
    vals = line.split()
    u = vals[h["union"]]
    svlen = int(vals[h["svLen"]])
    statVals = [svlen]    
    annCount = None
    if "phc_count" in h:
        annCount = int(vals[h["phc_count"]])
        statVals.append(annCount)
        nStats = 2

    if "All" in u or "PacBio" in u:
        pbv.append(statVals)
    if "All" in u or "Illumina" in u:
        ilv.append(statVals)
    if "All" in u or "BioNano" in u:
        bnv.append(statVals)

    if "PacBio" in u and "Illumina" not in u:
        pbu.append(statVals)
    if "Illumina" in u and "PacBio" not in u:
        ilu.append(statVals)

def GetStats(valList):

    if len(valList) == 0:
        return "\t".join(["0"]*nStats)
    else:
        nv = len(valList)
        valStrs = []
        for i in range(0,nStats):

            total =sum(v[i] for v in  valList)
            avg = float(total)/nv
            valStrs.append(str(total))
            valStrs.append("{:2.2f}".format(avg))
            
        return "\t".join(valStrs)


print("PB\t{}\t".format(len(pbv)) + GetStats(pbv))
print("PB-unique\t{}\t".format(len(pbu)) + GetStats(pbu))
print("IL\t{}\t".format(len(ilv)) + GetStats(ilv))
print("IL-unique\t{}\t".format(len(ilu)) + GetStats(ilu))
print( "BN\t{}\t".format(len(bnv)) + GetStats(bnv))

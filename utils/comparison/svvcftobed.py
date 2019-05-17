#!/usr/bin/env python

import sys

import argparse

ap = argparse.ArgumentParser(description="Transform a VCF with variants to BED with. The coordinates of the insertion are start, start+length, and deletion, start, end")
ap.add_argument("vcf", help="Input vcf file.")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

vcf = open(args.vcf)
colIdx = {}
outFile = open(args.out, 'w')

def GetSV(info):
    svtype = None
    svlen = 0
    svop = None
    for kv in info:
        kvp = kv.split("=")
        if (kvp[0] == "SVTYPE"):
            svtype = kvp[1]
        if (kvp[0] == "SVLEN"):
            svlen = abs(int(kvp[1]))
    if (svtype == "DUP" or svtype == "INS" or svtype=="insertion"):
        svop = "insertion"
    elif (svtype == "DEL" or svtype =="deletion"):
        svop = "deletion"
    elif (svtype == "OTHER"):
        svop = "other"
    elif (svtype == "inversion"):
        svop = "inversion"
    else:
        print("ERROR. SV type must be either DEL,DUP,or INS, got " + svtype)
        sys.exit(1)
    return (svop, svlen)

for line in vcf:
    if (line[0:2] == "##"):
        continue
    if (line[0] == "#"):
        cols=line[1:].split()
        colIdx = { cols[i]: i for i in range(0,len(cols)) }
        if ("INFO" not in colIdx):
            print("Input vcf requires INFO field.")
            sys.exit(1)
        infoIdx = colIdx["INFO"]
        continue
    vals = line.split()
    info = vals[infoIdx].split(";")
    (svop, svlen) = GetSV(info)
    outFile.write( vals[0] + "\t" + vals[1] + "\t" + str(svlen + int(vals[1])) + "\t" + svop + "\n")
    
outFile.close()

#!/usr/bin/env python


import argparse

ap = argparse.ArgumentParser(description="Given an Illumina VCF that has been intersected with pacbio and loj, collapse for most similar overlap of the same type ")
ap.add_argument("--isect", help="Intersection set.", default="/dev/stdin")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--svDist", help="Consider distance to sv in this column", type=int, default=18)
ap.add_argument("--maxSVDist", help="Consider two SVs the same if they are within this distance, and appropriate size threshold", default=1000,type=int)
ap.add_argument("--maxRatioDiff", help="PacBio/Illumina must be within 1-mrd < len < 1+mrd",type=float, default=0.25)
ap.add_argument("--svSeq", help="Include SV sequence from this field.", default=None, type=int)
ap.add_argument("--svClass", help="Index of sv class", default=12,type=int)

args = ap.parse_args()

def StoreKV(kvpString):
    kvpairs=kvpString.split(";")
    allKvp={}
    for kvStr in kvpairs:
        kv = kvStr.split("=")
        if len(kv) == 2:
            allKvp[kv[0]] = kv[1]
    return allKvp

isectFile = open(args.isect)

curSet = []
prevChrom = None
prevStart = None
nFields=4

if args.svSeq is not None:
    nFields+=1
    
def GroupSet(curSet):
    kvp = StoreKV(curSet[0][7])
    svEnd = int(kvp["END"])
    #
    svOp = curSet[0][4]
    svChrom = curSet[0][0]
    svStart = int(curSet[0][1])
    
    # Start out with simple deletions or insertions
    check = "NOT_SET"
    if svOp == "<DEL>":
        check = ["deletion", "DEL"]
    elif svOp == "<INS>" or svOp == "<DUP>":
        check = ["insertion", "INS"]

    optIndex = None
    maxFrac = 0
    for i in range(0,len(curSet)):
        
        if curSet[i][args.svClass] not in check:
            continue
        svDist = int(curSet[i][args.svDist])
        print("dist: " + str(svDist) + str(curSet[i][0:3]) + " sv: " + str(curSet[i][9:12]))
        if svDist > args.maxSVDist:
            continue
        
        isectStart = int(curSet[i][9])
        isectEnd   = int(curSet[i][10])
        ovpFrac = float(isectEnd - isectStart) / int(kvp["SVLEN"])
        #            print svOp, check, curSet[i][12], ovpFrac, kvp["SVLEN"], svChrom, svStart, svEnd, curSet[i][8:10], isectEnd-isectStart
        if ovpFrac > maxFrac and \
           1 - args.maxRatioDiff < ovpFrac and \
           1 + args.maxRatioDiff > ovpFrac:
            optIndex = i
            maxFrac = ovpFrac

    if optIndex is not None:
        return (maxFrac, curSet[optIndex])
    else:
        return (None,None)
            
prevSVLen = None
prevSVEnd = None
prevSVClass=None
for line in isectFile:
    vals = line.split()
    kvp = StoreKV(vals[7])

    if prevChrom is not None and ( vals[0] != prevChrom or \
                                   vals[1] != prevStart or \
                                   vals[4] != prevSVClass or \
                                   kvp["END"] != prevSVEnd or \
                                   kvp["SVLEN"] != prevSVLen):
        # curSet is a group of SVs that all have the same coordinate and type
        (maxFrac, grouped) = GroupSet(curSet)
        if grouped is not None:
            res = grouped[8:11] + ["{:.2}".format(maxFrac)]    
            if args.svSeq:
                res.append(grouped[args.svSeq])
            print("\t".join(res))
        else:
            print("\t".join(["."]*nFields))
        curSet = [vals]
    else:
        curSet.append(vals)
    prevChrom = vals[0]
    prevStart = vals[1]
    prevSVLen = kvp["SVLEN"]
    prevSVEnd = kvp["END"]
    prevSVClass = vals[4]
    

(maxFrac, grouped) = GroupSet(curSet)

if grouped is not None:
    res = grouped[8:11] + ["{:.2}".format(maxFrac)]    
    if args.svSeq:
        res.append(grouped[args.svSeq])
    print("\t".join(res))
else:
    print("\t".join(["."]*nFields))

#!/usr/bin/env python
import sys


import argparse

ap = argparse.ArgumentParser(description="Modify headers and write out")
ap.add_argument("--source", help="Source files.", nargs="+", required=True)
ap.add_argument("--append", help="Append these arguments", nargs="+", default=[])
ap.add_argument("--index", help="Add index to other source files after the first.", action='store_true', default=False)
args = ap.parse_args()


header=[]
index = 1
for inFileName in args.source:
    inFile = open(inFileName)
    oneHeader = inFile.readline()
    oneHeaderVals = oneHeader[1:].split()
    if index > 1 and args.index:
        
        newHeader = []
        for i in range(0,len(oneHeaderVals)):
            newHeader.append(oneHeaderVals[i] + "_" + str(index))
        oneHeaderVals = newHeader
    header = header + oneHeaderVals
    
    index +=1
header += args.append
print("#"+ "\t".join(header))

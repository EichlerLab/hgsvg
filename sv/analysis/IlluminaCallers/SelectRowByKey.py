#!/usr/bin/env python

import sys
import argparse

ap = argparse.ArgumentParser(description="Select rows from one file based on coordinate key in anohter, assumes uniqueness of coordinate keys")
ap.add_argument("--keyFile", help="Input key file, default column is 'key'")
ap.add_argument("--keyColumn", help="Column for key in key file.", default="key")
ap.add_argument("--dbFile", help="Database file.", required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

dbFile = open(args.dbFile)
outFile =open(args.out,'w')
dbHeader = []
db = {}
dbHeaderLine = dbFile.readline()
dbHeaderVals = dbHeaderLine.split()
dbHeaderVals = [i+"_db" for i in dbHeaderVals]
outFile.write("\t".join(dbHeaderVals )+ "\n")
for line in dbFile:
    vals = line.split()
    dbKey = "_".join(vals[0:3])
    db[dbKey] = line

keyFile = open(args.keyFile)
keyHeaderVals = keyFile.readline().split()
keyHeader = {keyHeaderVals[i] : i for i in range(0,len(keyHeaderVals))}
if args.keyColumn not in keyHeader:
    print("ERROR, key file " + args.keyFile + " is missing key column " + args.keyColumn)
    sys.exit(1)

missingLine = "\t".join(["NA"]*len(dbHeaderVals)) + "\n"
for line in keyFile:
    vals = line.split()
    key = vals[keyHeader[args.keyColumn]]
    if key in db:
        outFile.write(db[key])
    else:
        outFile.write(missingLine)





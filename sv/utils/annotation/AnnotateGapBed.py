#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

if (len(sys.argv) < 3):
    print "usage: AnnotateGapBed.py bedIn bedOut annotation.out"
    sys.exit(0)

ap = argparse.ArgumentParser(description="Print gap sequences to fasta files.")
ap.add_argument("bedin", help="Input bed file.")
ap.add_argument("bedout", help="Output bed file.")
ap.add_argument("dotout", help="RepeatMasker file.out annotation file.")
ap.add_argument("maskedout", help="Masked output file.", default=None)
ap.add_argument("--seqidx", help="Index of gap sequence (6)", default=6, type=int)
args = ap.parse_args()
bedFileIn = open(args.bedin, 'r')
bedFileOut = open(args.bedout, 'w')

annotations = {}

dotoutFile = open(args.dotout, 'r')
maskedDict = {}

if (args.maskedout is not None):
    maskedSequences = open(args.maskedout)
    maskedDict = SeqIO.to_dict(SeqIO.parse(maskedSequences, "fasta"))


for i in range(3):
    dotoutFile.readline()

for line in dotoutFile:
    vals = line.split()
    name = vals[4]
    rep  = vals[9]
    pre  = vals[11]
    post = vals[13]
    pre = int(pre.replace("(","").replace(")",""))
    post = int(post.replace("(","").replace(")",""))
    
    if (pre+post < 30):
        rep = rep + ":FULL"
    else:
        rep = rep + ":INC"
    
    if (name not in annotations):
        annotations[name] = []
    annotations[name].append(rep)

cols = {}
for line in bedFileIn:
    if line[0] == "#":
        hv = line[1:].split()
        for i in range(0,len(hv)):
            cols[hv[i]] = i
        
        bedFileOut.write(line.rstrip() + "\tsvAnn\tsvRep\n"  )
        continue
    vals = line.split()
    name = '/'.join(vals[0:3])
    if (name in annotations):
        annotation = ','.join(annotations[name])
    else:
        annotation = "NONE"

    repeatContent = ""
    if (name in maskedDict):
        svSeqi=cols["svSeq"]
        vals[svSeqi] =  maskedDict[name].seq.tostring()
        repeatContent = "{:2.2f}".format(float(vals[svSeqi].count("a") + vals[svSeqi].count("c") + vals[svSeqi].count("g") + vals[svSeqi].count("t"))/len(vals[svSeqi]))
        
    line = '\t'.join(vals) + "\t" + annotation + "\t" + repeatContent + '\n'
    
    bedFileOut.write(line)

        
bedFileOut.close()    
        
        

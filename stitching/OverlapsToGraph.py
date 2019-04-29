#!/usr/bin/env python
import sys
import networkx as nx
import argparse

import Overlap

ap = argparse.ArgumentParser(description="Convert overlaps to a graph, and remove transitively encoded edges.")
ap.add_argument("overlaps", help="Original overlaps file.")
ap.add_argument("--maxIndel", help="Maximum indels to consider.", type=int,default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
args = ap.parse_args()

#chr13:18849420-18929420/0       0       0       chr13:18829420-18909420/1       0       0       12222   12651   12953   13385
overlapFile = open(args.overlaps)
overlaps = []
i=1
for line in overlapFile:
    # hack to get around blank lines
    if len(line) > 1:
        overlaps.append(Overlap.Overlap(line))
    i+=1

g = nx.DiGraph()
ovpIndex = 0
t='chr5:149787757-149867757/0'
for ovp in overlaps:
    if ovp.a == t:
        print("src: " + t)
    if ovp.OverlapLength() > 0 and \
       ovp.Contained() == False and \
       ovp.Extends(wiggle=1000) == True and \
       (args.maxIndel == None or ovp.indel <= args.maxIndel):
        g.add_edge(ovp.a, ovp.b, index=ovpIndex)
        if ovp.a == t:
            print(t + " overlaps " + ovp.b)
    else:
        if ovp.a == t:
            print(t + " not overlaps " + ovp.b)
                
    ovpIndex+=1

nx.write_gml(g, args.out)

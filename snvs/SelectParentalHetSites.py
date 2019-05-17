#!/usr/bin/env python

import sys


import argparse

ap = argparse.ArgumentParser(description="Select sites in parents that are het when the child is phased and het")
ap.add_argument("--vcf", help="VCF of trio.",required=True)
ap.add_argument("--fa", help="Name of father.",required=True)
ap.add_argument("--mo", help="Name of mother.",required=True)
ap.add_argument("--ch", help="Name of child", required=True)

args = ap.parse_args()


faOut = open(args.fa + ".het.vcf",'w')
moOut = open(args.mo+ ".het.vcf",'w')

vcf = open(args.vcf)
ifa = 0
imo = 0
ich = 0
gts = ["0|1", "1|0", "0/1", "1/0"]
for line in vcf:
    if line[0:2] == "##":
        faOut.write(line)
        moOut.write(line)
        continue
    else:
        vals = line.split()
        if vals[0] == "#CHROM":
            faOut.write("\t".join(vals[0:8]) + "\t" + args.fa + "\n")
            moOut.write("\t".join(vals[0:8]) + "\t" + args.mo + "\n")            
            for i in range(1,len(vals)):
                if vals[i] == args.fa:
                    ifa = i
                elif vals[i] == args.mo:
                    imo = i
                elif vals[i] == args.ch:
                    ich = i
            if ifa == 0 or imo == 0 or ich == 0:
                print("Did not find a sample ")
                sys.exit(0)
            continue
        else:
            fa = vals[ifa].split(":")[0]
            mo = vals[imo].split(":")[0]
            ch = vals[ich].split(":")[0]
            if "|" in ch:
                if fa in gts:
                    faOut.write("\t".join(vals[0:8]) + "\t" + vals[ifa] + "\n")
                if mo in gts:
                    moOut.write("\t".join(vals[0:8]) + "\t" + vals[ifa] + "\n")
            
            
                    

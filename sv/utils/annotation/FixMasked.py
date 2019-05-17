#!/usr/bin/env python
import sys
idx=11
if (len(sys.argv) == 2):
	idx = int(sys.argv[1])

while sys.stdin:
	line = sys.stdin.readline()
        if line[0] == "#":
            vals = line[1:].split()
            h = {vals[i]:i for i in range(0,len(vals))}
            sys.stdout.write(line)
            continue
	if (line == ""):
		break
	v = line.split()
	nl = sum([v[h["svSeq"]].count(i) for i in ['a','g','c','t']])
	t = len(v[h["svSeq"]])
	v[h["svRep"]] = "{:2.2f}".format(float(nl)/t)
	sys.stdout.write('\t'.join(v) + '\n')

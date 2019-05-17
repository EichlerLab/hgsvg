#!/usr/bin/env python

import sys
for line in sys.stdin:
    if line[0] == '#':
        continue
    v = line.split()
    ps = None
    if len(v) < 9:
        continue
    if "PS" in v[8]:
        v2 = v[8].split(":")
        psi=None
        
        for j in range(0,len(v2)):
            if v2[j] == "PS":
                psi = j
                break
        if psi is not None:
            vals = v[9].split(":")
            if (vals[0] == "0|1" or vals[0] == "1|0") and psi < len(vals):
                ps=vals[psi]

        if ps is not None:
            print(v[0] + "\t" + v[1] + "\t" + ps)

if ps is not None:
    print(v[0] + "\t" + v[1] + "\t" + ps)

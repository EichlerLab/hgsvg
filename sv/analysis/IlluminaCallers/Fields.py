#!/usr/bin/env python
import sys
for line in sys.stdin:
    vals = line.split()
    nMethods = len(vals[0].split(","))
    print(str(nMethods) + "\t" + vals[1])

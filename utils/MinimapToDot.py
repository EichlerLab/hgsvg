#!/usr/bin/env python
import sys
#stitch.006.1912.45383312        45383312        795604  796317  +       chr11   135086622       4449496 4450209 452     713     255     cm:i:62
#1	3662697	51	0	0
for line in sys.stdin:
    v = line.split()
    if v[4] == "+":
	strand = 0
    else:
	strand  = 1
    print("{}\t{}\t{}\t{}\t{}".format(v[2],v[7], int(v[3]) - int(v[2]), strand, 0))
    

#! /usr/bin/env python
# compute stat on depth cov file from samtools depth
import sys
import numpy as np

f = open(sys.argv[1])
cmin = int(sys.argv[2])
cmax = int(sys.argv[3])


for line in f:
    line = line.rstrip()
    arline = line.split("\t")
    cov = int(arline[2])

    if cov > cmax or cov < cmin:
        print(line)


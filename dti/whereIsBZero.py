#!/usr/bin/env python
import sys
with open(str(sys.argv[1]), 'r') as ff:
    bval_line = ff.readlines()[0]

bval_list = bval_line.split()
print bval_list.index('0')

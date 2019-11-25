#!/usr/bin/env python

import os, sys

"""pass in a list of abs_mean.rms file path, print the path with the smallest abs_mean.rms value"""

# If mcflirt has been performed
path_list = list(sys.argv[1:])
smallest_val = None
smallest_path = None
for path in path_list:
    with open(path, 'rb') as handle:
        abs_mean_rms = float(handle.read())
        if (smallest_val != None):
            if (abs_mean_rms < smallest_val):
                smallest_val = abs_mean_rms
                smallest_path = path
        else:
            smallest_val = abs_mean_rms
            smallest_path = path
        
print smallest_path, smallest_val


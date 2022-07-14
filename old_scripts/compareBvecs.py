#!/usr/bin/env python

import os, sys
import math

path_1 = str(sys.argv[1])
path_2 = str(sys.argv[2])

with open(path_1, 'r') as fh:
    bvecs_1 = fh.readlines()
with open(path_2, 'r') as fh:
    bvecs_2 = fh.readlines()

for row in range(3):
    bvecs_1[row] = bvecs_1[row].split() # convert from space-separated string to list.
    bvecs_2[row] = bvecs_2[row].split()
    if (len(bvecs_1[row]) != len(bvecs_2[row])):
        print "Error: Different number of vectors in input files."
        sys.exit()

max_abs_diff = 0.0 # largest absolute difference between corresponding components of any of the vectors.
num_vecs = len(bvecs_1[0])

for row in range(3):
    for col in range(num_vecs):
        # Convert list of strings to list of floats
        bvecs_1[row][col] = float(bvecs_1[row][col])
        bvecs_2[row][col] = float(bvecs_2[row][col])

        abs_diff = abs( bvecs_1[row][col] - bvecs_2[row][col])
        if (abs_diff > max_abs_diff):
            max_abs_diff = abs_diff

max_len_diff_vec = 0.0
max_len_diff_flip_vec = 0.0
for col in range(num_vecs):
    # Compute the difference between the two bvecs.
    diff = [0.0, 0.0, 0.0]
    diff[0] = bvecs_1[0][col] - bvecs_2[0][col]
    diff[1] = bvecs_1[1][col] - bvecs_2[1][col]
    diff[2] = bvecs_1[2][col] - bvecs_2[2][col]

    # Compute the difference between the two bvecs after flipping one of them.
    diff_neg = [0.0, 0.0, 0.0]
    diff_neg[0] = bvecs_1[0][col] + bvecs_2[0][col]
    diff_neg[1] = bvecs_1[1][col] + bvecs_2[1][col]
    diff_neg[2] = bvecs_1[2][col] + bvecs_2[2][col]

#    print diff
    
    # Get length of the difference vector:
    len_diff_vec = math.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
    len_diff_neg_vec = math.sqrt(diff_neg[0]**2 + diff_neg[1]**2 + diff_neg[2]**2)
    len_diff_flip_vec = min(len_diff_vec, len_diff_neg_vec) # 

    if (len_diff_vec > max_len_diff_vec):
        max_len_diff_vec = len_diff_vec

    if (len_diff_flip_vec > max_len_diff_flip_vec):
        max_len_diff_flip_vec = len_diff_flip_vec
        
#        print '\033[93m' + "("+str(bvecs_1[0][col])+", "+str(bvecs_1[1][col])+", "+str(bvecs_1[2][col])+ ") - (" +str(bvecs_2[0][col])+", "+str(bvecs_2[1][col])+", "+str(bvecs_2[2][col])+") = (" + str(diff[0])+", "+str(diff[1])+", "+str(diff[2])+")" + '\033[0m'

#print "Largest absolute difference between vector components:", max_abs_diff
print "Length of largest difference vector (bvec1 - bvec2)  :", max_len_diff_vec
print "Length of largest difference vector allowing flips   :", max_len_diff_flip_vec

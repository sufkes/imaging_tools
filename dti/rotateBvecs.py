#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np

def main(affine_path, bvec_in_path, bvec_out_path):
    M = np.loadtxt(affine_path) 
    A = M[:3,:3] # rotate, stretch, shear portion of matrix?
    X = np.loadtxt(bvec_in_path)

    Y = A.dot(X)

    for col in range(X.shape[1]):
        print(np.linalg.norm(X[:, col]))
        print(np.linalg.norm(Y[:, col]))
    
    np.savetxt(bvec_out_path, Y, delimiter='  ', fmt='%.10g')
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("affine_path", help="path to affine matrix file (a 4x4 matrix with two-space-separated values, as output by FSL FLIRT).")
    parser.add_argument("bvec_in_path", help="path to input bvec file")
    parser.add_argument("bvec_out_path", help="path to output bvec file")
    
    # Define optional arguments.
#    parser.add_argument("-", "--", help="")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

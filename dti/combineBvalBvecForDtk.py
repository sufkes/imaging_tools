#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np

def main(bval_path, bvec_path, out_path):
    bval = np.loadtxt(bval_path, dtype=str)
    bvec = np.loadtxt(bvec_path, dtype=str)
    
    
    print(bvec.T.shape)
    print(bval.reshape(bval.shape[0],1).shape)
    combined = np.append(bvec.T, bval.reshape(bval.shape[0],1), axis=1)
    print(combined)
    print(combined.shape)

    np.savetxt(out_path, combined, delimiter=',', fmt='%s')
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Combine BVAL and BVEC files output by dcm2niix into a single CSV file with 4 columns and one row per DTI volume, as required by dti_tracker from the Diffusion Toolkit."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("bval_path", help="input BVAL path")
    parser.add_argument("bvec_path", help="input BVEC path")
    parser.add_argument("out_path", help="output combined BVEC & BVAL path")
    
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

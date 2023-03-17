#!/usr/bin/env python

import os
import sys
import argparse

def main():
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="path to input NIFTI file with integer labels")
    parser.add_argument("out_path", help="path to output CSV file to which coordinates of each label's centre-of-gravity')
    
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

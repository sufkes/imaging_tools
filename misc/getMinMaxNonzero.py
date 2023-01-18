#!/usr/bin/env python

import os
import sys
import argparse

import nibabel as nib
import numpy as np


def main(in_path):
    nii = nib.load(in_path)
    nii_array = nii.get_fdata()
    try:
        min_nonzero = nii_array[nii_array != 0].min()
        max_nonzero = nii_array[nii_array != 0].max()
    except ValueError: # nii_array might be all zero, in which case, return Nones.
        min_nonzero = None
        max_nonzero = None
    print(f'{min_nonzero} {max_nonzero}')
    return (min_nonzero, max_nonzero)

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Output the minimum value across nonzero voxels."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="path to NIFTI file")
    
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

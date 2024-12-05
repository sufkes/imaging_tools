#!/usr/bin/env python

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path, include_only, exclude_only):
    if (not include_only is None) and (not exclude_only is None):
        msg = 'Use either -i/--include_only or -e/--exclude_only, not both.'
        raise Exception(msg)

    ## Load the original label image.
    img_nii = nib.load(in_path)
    img_arr = np.round(img_nii.get_fdata()).astype(int)

    ## Set values to zero.
    if (not include_only is None):
        img_arr[~np.isin(img_arr, include_only)] = 0

    elif (not exclude_only is None):
        img_arr[np.isin(img_arr, exclude_only)] = 0

    ## Save the modified image.
    out_img_nii = nib.nifti1.Nifti1Image(img_arr, img_nii.affine, img_nii.header)
    nib.save(out_img_nii, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Set specified values in a NIFTI file to zero (e.g. to remove specific labels from an atlas). Input image should contain integer labels. Use either the -i/--include_only or -e/--exclude_only flag; not both."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="input NIFTI path")
    parser.add_argument("out_path", help="output NIFTI path")
    
    # Define optional arguments.
    parser.add_argument("-i", "--include_only", help="set of integers which will be kept; all others will be set to zero", nargs="+", type=int)
    parser.add_argument("-e", "--exclude_only", help="set of integers which will be set to zero; all others will be kept", nargs="+", type=int)

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

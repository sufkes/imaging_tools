#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_paths):
    ## Load the 4D image.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    ## Check that the number of volumes matches the number of output paths.
    if len(img_arr.shape) != 4:
        raise Exception("Input image must be 4D.")
    if img_arr.shape[3] != len(out_paths):
        raise Exception(f"Number of volumes in input image ({img_arr.shape[3]}) must match number of output paths ({len(out_paths)}).")

    ## Save each volume to a separate file.
    for index, out_path in enumerate(out_paths):
        vol_arr = img_arr[:, :, :, index]
        
        vol_nii = nib.nifti1.Nifti1Image(vol_arr, img_nii.affine, header=img_nii.header)

        ## Save the padded image.
        nib.save(vol_nii, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Save each volume in a 4D NIFTI image to a separate file.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file')
    parser.add_argument('out_paths', help='path to output NIFTI file', nargs='+')
    
    # Define optional arguments.

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

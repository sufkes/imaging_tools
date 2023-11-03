#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import nibabel as nib

def main(in_paths, out_path, threshold):

    first_iteration = True
    for in_path in in_paths:
        vis_nii = nib.load(in_path) # load visitation map NIFTI file
        vis_arr = vis_nii.get_fdata()  # load visition map array data

        vis_arr[vis_arr > 0] = 1 # binarize

        if first_iteration:
            # Generate the disconnectome array.
            dis_arr = vis_arr.copy()
            dis_affine = vis_nii.affine
            dis_header = vis_nii.header
            first_iteration = False
        else:
            dis_arr += vis_arr
        

    # Convert the sum of visitations to the frequency of visitations.
    num_controls = len(in_paths)
    dis_arr = dis_arr/num_controls

    # Threshold the disconnectome map.
    dis_arr[dis_arr < threshold] = 0
    
    # Convert disconnectome array to a NIFTI object.
    dis_nii = nib.nifti1.Nifti1Image(dis_arr, dis_affine, header=dis_header)
    
    # Save the disconnectome map.
    nib.save(dis_nii, out_path)

    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Generate a disconnnectome map from a list of visitation maps. Visition maps are binarized, averaged, then the output is thresholded to remove spurious disconnections."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_path", help="disconnectome output path, in NIFTI format", required=True, type=str)
    parser.add_argument("-i", "--in_paths", help="input visitation maps, in NIFTI format", nargs='+', required=True, type=str)
    parser.add_argument("-t", "--threshold", help="threshold probability for disconnectome map: zero regions in the disconnectome map with visitation frequency across control subjects below this number", type=float, default=0.5)

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

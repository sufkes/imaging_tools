#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np
import pandas as pd

def main(in_paths, label_path, out_path):

    # Initialize the output CSV.
    df = pd.DataFrame()
    df.index.name = 'image_path'
    
    # Open the label image.
    label_nii = nib.load(label_path)
    label_arr = label_nii.get_fdata().astype(int)
    label_values = np.unique(label_arr)
    label_values.sort()
    
    # Loop over input images.
    for in_path in in_paths:
        image_nii = nib.load(in_path)
        image_arr = image_nii.get_fdata()

        # Check that image is in the same space as the label image.
        if image_arr.shape != label_arr.shape:
            msg = f'Input image shape does not match label image.'
            raise Exception(msg)
        if (image_nii.affine != label_nii.affine).any():
            msg = f'Warning: Input image affine does not match label image affine.'
            print(msg)
        
        for label_val in label_values:
            mean = image_arr[label_arr==label_val].mean()
            df.loc[in_path, f'mean_label_{label_val}'] = mean

    # Save or print result
    if out_path is None:
        print(df.to_string())
    else:
        df.to_csv(out_path, index=True)
    
if (__name__ == '__main__'):
    # Create argument parser.
    description = """Input NIFTI images and compute their mean values within ROIs defined in a separate label image."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    
    
    # Define optional arguments.
    parser.add_argument('-i', '--in_paths', help='paths to NIFTI files', nargs='+')
    parser.add_argument('-l', '--label_path', help='paths to NIFTI label image; must be 3D image with integer labels in the same space as the input images.', required=True)
    parser.add_argument('-o', '--out_path', help='path to output CSV file storing the mean values; If not specified, print values to console and do not save.', required=False)

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

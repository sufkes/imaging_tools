#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path, x_axis, y_axis, z_axis):
    ## Load the original image and extra the data array.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    ## Store the header and affine from the original file.
    affine = img_nii.affine
    header = img_nii.header

    ## Perform the reflections.
    if x_axis:
        img_arr = np.flip(img_arr, axis=0)
    if y_axis:
        img_arr = np.flip(img_arr, axis=1)
    if z_axis:
        img_arr = np.flip(img_arr, axis=2)

    img_nii = nib.Nifti1Image(img_arr, affine, header)
        
    ## Save the padded image.
    nib.save(img_nii, out_path)
    return img_nii

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Reflect NIFTI image in specified axis. Input image must be 3D.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file')
    parser.add_argument('out_path', help='path to output NIFTI file')
    
    # Define optional arguments.
    parser.add_argument('-x', '--x_axis', action='store_true', help='Reverse the order of elements in the x-axis.')
    parser.add_argument('-y', '--y_axis', action='store_true', help='Reverse the order of elements in the y-axis.')
    parser.add_argument('-z', '--z_axis', action='store_true', help='Reverse the order of elements in the z-axis.')


    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path, x_unpad, y_unpad, z_unpad, t_unpad):
    ## Load the original image and extra the data array.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    ## Unpad the data array.
    if len(img_arr.shape) == 3:
        img_arr_unpad = img_arr[x_unpad[0]:img_arr.shape[0]-x_unpad[1], y_unpad[0]:img_arr.shape[1]-y_unpad[1], z_unpad[0]:img_arr.shape[2]-z_unpad[1]]
    elif len(img_arr.shape) == 4:
        img_arr_unpad = img_arr[x_unpad[0]:img_arr.shape[0]-x_unpad[1], y_unpad[0]:img_arr.shape[1]-y_unpad[1], z_unpad[0]:img_arr.shape[2]-z_unpad[1], t_unpad[0]:img_arr.shape[3]-t_unpad[1]] 
    else:
        raise Exception('Input image must be either 3D or 4D.')
    
    ## Generate a new NIFTI file.
    affine_unpad = img_nii.affine
    header_unpad = img_nii.header
    img_nii_unpad = nib.nifti1.Nifti1Image(img_arr_unpad, affine_unpad, header=header_unpad)

    ## Save the unpadded image.
    nib.save(img_nii_unpad, out_path)
    return img_nii_unpad

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Remove padding slices from NIFTI image. Input image must be 3D or 4D. The image affine will be unchanged, which is different behaviour from fslroi.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file')
    parser.add_argument('out_path', help='path to output NIFTI file')
    
    # Define optional arguments.
    parser.add_argument('-x', '--x_unpad', nargs=2, type=int, metavar=('UNPAD_X_START', 'UNPAD_X_END'), default=(0,0), help='number of slices to remove from the start and end of the first (x) dimension')
    parser.add_argument('-y', '--y_unpad', nargs=2, type=int, metavar=('UNPAD_Y_START', 'UNPAD_Y_END'), default=(0,0), help='number of slices to remove from the start and end of the second (y) dimension')
    parser.add_argument('-z', '--z_unpad', nargs=2, type=int, metavar=('UNPAD_Z_START', 'UNPAD_Z_END'), default=(0,0), help='number of slices to remove from the start and end of the third (z) dimension')
    parser.add_argument('-t', '--t_unpad', nargs=2, type=int, metavar=('UNPAD_T_START', 'UNPAD_T_END'), default=(0,0), help='number of volumes to remove from the start and end of the fourth (time/volume) dimension')

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

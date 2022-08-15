#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path, x_pad, y_pad, z_pad, t_pad):
    ## Load the original image and extra the data array.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    ## Pad the data array.
    # If image is 3D, and time/volume padding is not requested, leave it 3D.
    if len(img_arr.shape) == 3:
        if t_pad == (0,0):
            padding = (x_pad, y_pad, z_pad)
        else:
            img_arr = np.expand_dims(img_arr, axis=3)
            padding = (x_pad, y_pad, z_pad, t_pad)
    elif len(img_arr.shape) == 4:
        padding = (x_pad, y_pad, z_pad, t_pad)
    else:
        raise Exception('Input image must be either 3D or 4D.')
        
    img_arr_pad = np.pad(img_arr, padding)
    
    ## Generate a new NIFTI file.
    affine_pad = img_nii.affine
    header_pad = img_nii.header
    img_nii_pad = nib.nifti1.Nifti1Image(img_arr_pad, affine_pad, header=header_pad)

    ## Save the padded image.
    nib.save(img_nii_pad, out_path)
    return img_nii_pad

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Pad NIFTI image with zeros. Input image must be 3D or 4D.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file')
    parser.add_argument('out_path', help='path to output NIFTI file')
    
    # Define optional arguments.
    parser.add_argument('-x', '--x_pad', nargs=2, type=int, metavar=('PAD_X_START', 'PAD_X_END'), default=(0,0), help='number of slices to add the start and end of the first (x) dimension')
    parser.add_argument('-y', '--y_pad', nargs=2, type=int, metavar=('PAD_Y_START', 'PAD_Y_END'), default=(0,0), help='number of slices to add the start and end of the second (y) dimension')
    parser.add_argument('-z', '--z_pad', nargs=2, type=int, metavar=('PAD_Z_START', 'PAD_Z_END'), default=(0,0), help='number of slices to add the start and end of the third (z) dimension')
    parser.add_argument('-t', '--t_pad', nargs=2, type=int, metavar=('PAD_T_START', 'PAD_T_END'), default=(0,0), help='number of volumes to add the start and end of the fourth (time/volume) dimension')

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

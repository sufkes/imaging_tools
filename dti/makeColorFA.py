#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(FA_path, V1_path, out_path):
    fa_nifti = nib.load(FA_path)
    v1_nifti = nib.load(V1_path)

    fa_array = fa_nifti.get_fdata()
    v1_array = v1_nifti.get_fdata()

    ## Generate FA color image.
    # The FA color image should have the shape [dim1, dim2, dim3, color]. The color dimension should store values between 0 and 255.

    # For some reason, FSL's dtifit generates FA images with values > 1. Try simply setting them to 1.
    fa_array[fa_array>1] = 1
    fa_color_array = np.abs(v1_array) # use absolute value of principal eigenvector components to get left/right, posterior/anterior, superior/inferior component
    
    fa_color_array[:,:,:,0] = fa_color_array[:,:,:,0] * fa_array # multiple the direction components by FA
    fa_color_array[:,:,:,1] = fa_color_array[:,:,:,1] * fa_array
    fa_color_array[:,:,:,2] = fa_color_array[:,:,:,2] * fa_array

    fa_color_array *= 255 # convert to number between 0 and 255
    fa_color_array = np.rint(fa_color_array).astype(np.int8) # convert to integers
    
    ## Save the color FA map.
    # ras_pos is a 4-d numpy array, with the last dim holding RGB
    shape_3d = fa_color_array.shape[0:3]
    rgb_dtype = np.dtype([('R', 'u1'), ('G', 'u1'), ('B', 'u1')])
    fa_color_array = fa_color_array.copy().view(dtype=rgb_dtype).reshape(shape_3d)  # copy used to force fresh internal structure
    fa_color_nifti = nib.Nifti1Image(fa_color_array, fa_nifti.affine)
    nib.save(fa_color_nifti, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Generate a color fractional anisotropy (FA) image from an input FA image and principal eigenvector (V1) image. Both input images must be NIFTI files in the same image space."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("FA_path", type=str, help="path to fractional anisotropy (FA) NIFTI file")
    parser.add_argument("V1_path", type=str, help="path to principal eigenvector (V1) NIFTI file")
    parser.add_argument("out_path", type=str, help="path to output color FA map")
    
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

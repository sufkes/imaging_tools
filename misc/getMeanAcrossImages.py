#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_paths, out_path):
    first_image=True

    for index, in_path in enumerate(in_paths):
        ## Load the 3D image.
        img_nii = nib.load(in_path)
        img_arr = img_nii.get_fdata()

        if first_image:
            # Initialize the 4D array which will store all of the input images.
            all_images_shape = img_arr.shape + (len(in_paths),)
            all_images_arr = np.zeros(shape=all_images_shape, dtype=img_arr.dtype)

            # Copy header and affine to be used in mean image.
            mean_header = img_nii.header
            mean_affine = img_nii.affine

            first_image = False
        else:
            if img_arr.shape != all_images_shape[:3]:
                raise Exception("All input images must have the same dimensions")

        all_images_arr[:, :, :, index] = img_arr

    # Compute the mean image and save.
    mean_arr = all_images_arr.mean(axis=3)
    mean_nii = nib.nifti1.Nifti1Image(mean_arr, mean_affine, header=mean_header)
    nib.save(mean_nii, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Compute the mean across a series of 3D NIFTI images. The result is a 3D image. All input images must have the same dimensions.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('out_path', help='path to output NIFTI file')
    parser.add_argument('in_paths', help='path to input NIFTI file', nargs='+')
    
    # Define optional arguments.

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

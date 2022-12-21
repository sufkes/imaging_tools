#!/usr/bin/env python

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(image_path_1, image_path_2, mask_path):
    array_1 = nib.load(image_path_1).get_fdata().flatten()
    array_2 = nib.load(image_path_2).get_fdata().flatten()
    

    if mask_path is None:
        r = np.corrcoef(array_1, array_2)
    else:
        array_mask = nib.load(mask_path).get_fdata().flatten()
        array_1_masked = array_1[array_mask>0]
        array_2_masked = array_2[array_mask>0]

        r = np.corrcoef(array_1_masked, array_2_masked)[0, 1]

    print(r)
    return r

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Calculate the Pearson correlation coeffecient of voxel intensities between two images."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("image_path_1", type=str, help="path to first image")
    parser.add_argument("image_path_2", type=str, help="path to second image")
    
    # Define optional arguments.
    parser.add_argument("-m", "--mask_path", type=str, help="path to mask to which calculation will be restricted; if none specified, include all voxels in the input images")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

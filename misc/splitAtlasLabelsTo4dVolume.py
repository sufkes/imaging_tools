#!/usr/bin/env python

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path):
    ## Load the original label image.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    # Determine the number of labels in the image.
    unique_values = np.sort(np.unique(img_arr))
    unique_labels = unique_values[unique_values != 0]

    num_labels=len(unique_labels)
    
    print(f'Found {num_labels} unique labels. Labels will be mapped to the 4th dimension of the output image in the following order: {list(unique_labels)}')

    ## Generate the 4D array.
    shape_4d = img_arr.shape + (num_labels, )
    img_arr_4d = np.zeros(shape=shape_4d)

    # For loop method.
#    for i in range(num_labels):
#        old_label_value = unique_labels[i]
#        img_arr_4d[:, :, :, i][img_arr==old_label_value] = 1

    # Broadcasting method (faster?). 
    indices = np.arange(num_labels)
    img_arr_4d = (img_arr[..., np.newaxis] == unique_labels).astype(int)
        
    ## Generate the new NIFTI file.
    affine_4d = img_nii.affine
    header_4d = img_nii.header
    img_nii_4d = nib.nifti1.Nifti1Image(img_arr_4d, affine_4d, header=header_4d)

    ## Save the padded image.
    nib.save(img_nii_4d, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Convert a 3D NIFTI file with integer labels to a 4D NIFTI file with each label in a separate dimension. Assumes the labels are ordered 1, 2, ..., <number of labels>. Labels are assigned to the dimensions of the output image based on their value. For example, label=1 is mapped to the first volume, label=2 is mapped to the second volume, etc."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="path to input 3D NIFTI label file with integer labels")
    parser.add_argument("out_path", help="path to output 4D NIFTI label file with each label in a separate dimension")
    
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

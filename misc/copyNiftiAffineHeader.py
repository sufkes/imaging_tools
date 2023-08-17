#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib

def main(from_image_path, to_image_path, out_path, affine, header):
    from_image_nii = nib.load(from_image_path)
    to_image_nii = nib.load(from_image_path)

    # Select output image data, affine, and header.
    new_data = to_image.get_fdata() # We do not want to change the image data itself.
    
    if affine:
        new_affine = from_image_nii.affine
    else:
        new_affine = to_image_nii.affine

    if header:
        new_header = from_image_nii.header
    else:
        new_header = to_image_nii.header

    # Create the new NIFTI file using the chosen header and affine.
    new_nii = nib.nifti1.Nifti1Image(new_data, new_affine, header=new_header)

    # Save the modified image.
    if out_path is None:
        out_path = to_image_path
    nib.save(new_nii, out_path)
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Copy NIFTI affine or header from one image to another. By default, the image is modified in place."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("from_image_path", type=str, help='path to NIFTI file whose affine or header will be copied')
    parser.add_argument("to_image_path", type=str, help='path to NIFTI file whose affine or headers will be overwritten.')
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_path", help="new path of image whose affine or header was overwritten. If unspecified, the image is modified in place.", type=str)
    parser.add_argument("-a", "--affine", type='store_true', help='copy affine')
    parser.add_argument("-d", "--header", type='store_true', help='copy header')    

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib

def main(from_image_path, to_image_path, out_path):
    from_image_nii = nib.load(from_image_path)
    to_image_nii = nib.load(to_image_path)

    # Copy the qform and sform.
    to_image_nii.set_sform(from_image_nii.get_sform())
    to_image_nii.set_qform(from_image_nii.get_qform())

    # Save.
    nib.save(to_image_nii, out_path)
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Copy NIFTI sform and qform from one image onto another."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("from_image_path", type=str, help='path to NIFTI file whose qform and sform will be copied')
    parser.add_argument("to_image_path", type=str, help='path to NIFTI file whose qform and sform will be overwritten.')
    parser.add_argument("out_path", help="output path of NIFTI file with modified sform and qform", type=str)

    # Define optional arguments.
    

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

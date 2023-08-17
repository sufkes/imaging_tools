#!/usr/bin/env python

import os
import sys
import argparse

def main():
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Copy NIFTI affine or header from one image to another."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("from_image_path", type=str, help='path to NIFTI file whose affine or header will be copied')
    parser.add_argument("to_image_path", type=str, help='path to NIFTI file whose affine or headers will be overwritten.')
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_path", help="new path of image whose affine or header was overwritten. If unspecified, the image is modified in place.", type=str)

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

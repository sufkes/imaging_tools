#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import nibabel as nib

def main(in_path, out_path):
    in_nii = nib.load(in_path)
    data = in_nii.get_fdata()
    
    w_in, h_in, z_in, t_in = data.shape
    w = 256 # output width
    h = 256 # output height

    ## Fourier transform the input.
    a = np.fft.fft2(data, axes=(0,1))
    a = np.fft.fftshift(a, axes=(0,1))

    ## Zero-fill the high-frequencies to reach the desired shape.
    b = np.zeros((h, w, z_in, t_in), dtype=complex)
    b[int((h-h_in)/2):int((h+h_in)/2), int((w-w_in)/2):int((w+w_in)/2), :, :] = a 

    ## Inverse Fourier transform.
    b = np.fft.ifftshift(b, axes=(0,1))
    b = np.fft.ifft2(b, axes=(0,1))
    data_resamp = np.abs(b)

    # Copy the header from the input, and manually change the pixel width and height
    out_header = in_nii.header.copy()
    out_header['pixdim'][1:3] = [0.859400, 0.859400]

    out_nii = nib.Nifti1Image(data_resamp, in_nii.affine, header=out_header) # The q-form is apparently automatically set to "good" values during creation of the NIFTI object in the previous step, but the s-form does not change from its value which was copied from the input NIFTI. In the images resampled using AFNI, the sform and qform matched, so here, I will just copy the new q-form onto the s-form.
    out_nii.set_sform(out_nii.get_qform())

    # Save the resampled image.
    nib.save(out_nii, out_path)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Resample NIFTI file from shape (100, 100, H, T) to (256, 256, H, T) using zero-fill interpolation (ZIP). Basically, Fourier transform the slices, pad the high frequency terms with zeros until there are 256 x 256, and inverse Fourier transform. As far as I can tell, this script replicates the action of the GE scanner. This script could be improved to handle arbitrary input and output shapes."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("in_path", help="")
    parser.add_argument("out_path", help="")
    
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

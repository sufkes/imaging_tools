#!/usr/bin/env python3

import os
import sys
import argparse

import nibabel as nib
import numpy as np

def main(in_path, out_path, yz_line, above):
    ## Load the original image and extra the data array.
    img_nii = nib.load(in_path)
    img_arr = img_nii.get_fdata()

    # Zero points beyond the specified plane.
    img_arr_modified = img_arr

    if yz_line is not None:
        # Get equation of line in form z = m*y + b
        m = (yz_line[3] - yz_line[1]) / (yz_line[2]-yz_line[0]) # m = slope = "rise over run"

        b = yz_line[1] - m*yz_line[0] # b = z - m*y

        print(f'Zeroing points {("above" if above else "below")} line: z = {m}*y + {b}')

#        zvalues_2d = np.zeros(shape=img_arr_modified.shape[1:3])
 #       zvalues_2d[:,0] = range(img_arr_modified.shape[2]-1, -1, -1) # e.g. [3, 2, 1, 0] if z dim = 4
        zvalues_2d = np.array([range(img_arr_modified.shape[2])] * img_arr_modified.shape[1]) # e.g. [3, 2, 1, 0] if z dim = 4 # I should probably use something like np.repeat(a[:, :, np.newaxis], 3, axis=2) instead
#        zvalues_2d = np.rot90(zvalues_2d, k=3) # 2d array of z values with shape (dim y, dim z)

        zvalues_line_2d = np.array([m*np.array(range(img_arr_modified.shape[1])) + b] * img_arr_modified.shape[1]) # z value of line at each y value
        zvalues_line_2d = np.rot90(zvalues_line_2d, k=3)
        
        #print(zvalues_2d)
        #print(zvalues_line_2d)

        mask_2d = np.zeros(shape=zvalues_2d.shape, dtype=np.int64)
        if above: # if zeroing points above the plane
            mask_2d[zvalues_2d < zvalues_line_2d] = 1
        else:
            mask_2d[zvalues_2d > zvalues_line_2d] = 1

        #print(mask_2d.shape)
        mask_3d = np.repeat(mask_2d[np.newaxis, :, :], img_arr_modified.shape[0], axis=0) # now has shape (dim y, dim z, dim x)
        #print(mask_3d[0,0,:])
        #print(mask_3d[0,50,:])
        #print(mask_3d[20,0,:])
        #print(mask_3d[20,50,:])
        #print(mask_3d.shape)
        #print(img_arr_modified.shape)
        #print(mask_3d.shape)
        img_arr_modified[mask_3d==0] = 0
        #print(img_arr_modified.shape)

    ## Generate a new NIFTI file.
    affine_modified = img_nii.affine
    header_modified = img_nii.header
    img_nii_modified = nib.nifti1.Nifti1Image(img_arr_modified, affine_modified, header=header_modified)

    ## Save the modified image.
    nib.save(img_nii_modified, out_path)
    return img_nii_modified

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Zero 3D NIFTI image above or below a specified plane. Use to aid registration of images which cut off the brain at different planes.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file')
    parser.add_argument('out_path', help='path to output NIFTI file')
    
    # Define optional arguments.
    #parser.add_argument('-t', '--t_pad', nargs=2, type=int, metavar=('PAD_T_START', 'PAD_T_END'), default=(0,0), help='number of volumes to add the start and end of the fourth (time/volume) dimension')
    parser.add_argument('--yz_line', nargs=4, type=int, metavar=('Y1', 'Z1', 'Y2', 'Z2'), help='plane is constant in the x direction, specified by a line in the y-z plane. Enter two (y, z) pairs along the line.')
    parser.add_argument('-a', '--above', action='store_true', help='zero points above the specified plane. By default points below the specified plane are zeroed.')
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

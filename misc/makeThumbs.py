#!/usr/bin/env python3
import os
import sys
import argparse

import nibabel as nib
import numpy as np
from matplotlib import pyplot as plt

def main(in_path, out_path, volumes, mask_path, mask_cmap, mask_alpha, mask_binarize, mask_vmin):
    # Load image as numpy array.
    img_nii = nib.load(in_path)
    img_array = img_nii.get_fdata()
    if mask_path:
        mask_nii = nib.load(mask_path)
        mask_array= mask_nii.get_fdata()

    if mask_path and (img_array.shape != mask_array.shape):
        raise Exception(f'Main image shape: {img_array.shape} does not match mask image shape: {mask_array.shape}.')
        
    # Get the requested volume numbers only.
    if not volumes is None:
        img_array = img_array[:, :, :, volumes]
        if mask_path:
            mask_array = mask_array[:, :, :, volumes]

    # Binarize the mask image.
    if mask_binarize:
        mask_array[mask_array>0] = 1
        mask_array[mask_array<=0] = 0
            
            
    zooms = img_nii.header.get_zooms() # size of x, y, z[, TR]

    img_shape_real = np.array(img_array.shape)*np.array(zooms)
    
    ## Quick and dirty ortho mode plot, possibly with multiple volumes, each on a separate row.
    if len(img_array.shape) > 3:
        num_vols = img_array.shape[3]
    else:
        num_vols = 1

    #### Plot slices.
    ## Subplot sizing
    # Each subplot will have an aspect ratio which represents true lengths equally in all direction.
    # To ensure that lengths are represented the same in every subplot, set the width ratios using the real dimensions of the image frame.
    # Sagittal view has width of y dimension; height of z
    # Coronal view has width of x dimension; height of z
    # Axial view has width of x dimension; height of y
    width_ratios = [img_shape_real[1],
                    img_shape_real[0],
                    img_shape_real[0]]

    # Determine required size of figure. 
    fig_width_real_mm  = img_shape_real[1] + 2*img_shape_real[0] # width of figure if drawn to scale and without padding
    fig_height_real_mm = num_vols * max(img_shape_real[1], img_shape_real[2]) # height of figure if drawn to scale and without padding

    mm_per_inch = 25.4
    fig_width_real_inch = fig_width_real_mm/mm_per_inch
    fig_height_real_inch = fig_height_real_mm/mm_per_inch
    
    figsize = np.array([fig_width_real_inch, fig_height_real_inch])

    dpi = 92

    
    while True:
        fig_res = figsize * dpi
        if (fig_res > 10000).any():
            figsize = figsize/2
        else:
            break

    #print(f'Final figure size (inches) = {figsize}')
    
    cmap = 'gray'
    fig, axes = plt.subplots(nrows=num_vols, ncols=3, squeeze=False, figsize=figsize, gridspec_kw={'width_ratios': width_ratios})
    
    aspect_ratios = np.zeros((3,))
    aspect_ratios[0] = zooms[2]/zooms[1]
    aspect_ratios[1] = zooms[2]/zooms[0]
    aspect_ratios[2] = zooms[1]/zooms[0]
    
    for vol_index in range(num_vols):
        if num_vols > 1:
            plot_array = img_array[:,:,:,vol_index]
            if mask_path:
                mask_plot_array = mask_array[:,:,:,vol_index]
        else:
            plot_array = img_array
            if mask_path:
                mask_plot_array = mask_array

        midpoint = np.floor(np.array(plot_array.shape)/2).astype(int)

        axes[vol_index,0].imshow(np.rot90(plot_array[midpoint[0],:,:]), cmap=cmap, aspect=aspect_ratios[0])
        axes[vol_index,1].imshow(np.rot90(plot_array[:,midpoint[1],:]), cmap=cmap, aspect=aspect_ratios[1])
        axes[vol_index,2].imshow(np.rot90(plot_array[:,:,midpoint[2]]), cmap=cmap, aspect=aspect_ratios[2])

        if mask_path:
            from matplotlib import cm
            mask_cmap_modified = cm.get_cmap(mask_cmap).copy()
            mask_cmap_modified.set_under('k', alpha=0)
            axes[vol_index,0].imshow(np.rot90(mask_plot_array[midpoint[0],:,:]), cmap=mask_cmap_modified, aspect=aspect_ratios[0], alpha=mask_alpha, vmin=mask_vmin)
            axes[vol_index,1].imshow(np.rot90(mask_plot_array[:,midpoint[1],:]), cmap=mask_cmap_modified, aspect=aspect_ratios[1], alpha=mask_alpha, vmin=mask_vmin)
            axes[vol_index,2].imshow(np.rot90(mask_plot_array[:,:,midpoint[2]]), cmap=mask_cmap_modified, aspect=aspect_ratios[2], alpha=mask_alpha, vmin=mask_vmin)

    plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        
    for axis in fig.axes:
        axis.xaxis.set_visible(False)
        axis.yaxis.set_visible(False)
    
    #out_name = os.path.basename(in_path).split('.nii')[0] + '.png'
    #out_path = os.path.join(out_dir, out_name)
    #fig.set_tight_layout(True)
    fig.savefig(out_path, dpi=dpi)
    plt.close()
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Generate thumbnail of 3D NIFTI image (with optional overlay image) for rapid viewing.'''
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('in_path', help='path to input NIFTI file', type=str)
    parser.add_argument('out_path', help='path of output file; must have a compatible image file extension, e.g. ".png".', type=str)
    
    # Define optional arguments.
    #parser.add_argument('-m', '--display_mode', type=str, default='ortho', help='[NOT IMPLEMENTED] display mode (x: sagittal, y: coronal, z: axial, xyz: multiple slices in all three orientation, ortho: one slice in each direction)')
    #parser.add_argument('-n', '--num_cuts', help='[NOT IMPLEMENTED] number of slices to plot; ignored if display_mode is ortho')
    parser.add_argument('-v', '--volumes', type=int, nargs='+', help='volume numbers to include in figure (for 4D images); volumes should be separated by spaces; first volume has number 0.')
    parser.add_argument('--mask_path', type=str, help='path to mask NIFTI file; must have same dimensions as IN_PATH; mask is assumed to have same (tranformation matrix?) as main image.')
    parser.add_argument('--mask_cmap', type=str, default='Greens', help='mask overlay color map. Must be the name of a Matplotlib Colormap. See https://matplotlib.org/stable/tutorials/colors/colormaps.html for options.')
    parser.add_argument('--mask_alpha', type=float, default=0.5, help='mask overlay alpha (opacity). Must be value between 0 and 1.')
    parser.add_argument('--mask_binarize', action='store_true', help='binarize mask (set values > 0 to 1; values < 0 to 0.')
    parser.add_argument('--mask_vmin', type=float, help='minimum value to display for mask; values below with will have zero opacity.', default=0.1)
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))


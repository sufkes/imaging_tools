#!/usr/bin/env python3
import os
import sys
import argparse

import nibabel as nib
import numpy as np
from matplotlib import pyplot as plt

def main(in_path, out_path, mask_path, mask_cmap, mask_alpha, mask_binarize, mask_vmin, num_views, num_major_columns, view_spacing):
    # Load image as numpy array.
    img_nii = nib.load(in_path)
    img_array = img_nii.get_fdata()
    if mask_path:
        mask_nii = nib.load(mask_path)
        mask_array= mask_nii.get_fdata()

    if mask_path and (img_array.shape != mask_array.shape):
        raise Exception(f'Main image shape: {img_array.shape} does not match mask image shape: {mask_array.shape}.')
        
    # Binarize the mask image.
    if mask_binarize:
        mask_array[mask_array>0] = 1
        mask_array[mask_array<=0] = 0
            
            
    zooms = img_nii.header.get_zooms() # size of x, y, z[, TR]

    img_shape_real = np.array(img_array.shape)*np.array(zooms)
    
    ## Quick and dirty ortho mode plot, possibly with multiple volumes, each on a separate row.
    num_vols = 1

    # Determine number of rows and columns in figure
    ncols = 3 * num_major_columns
    nrows = np.ceil(num_views/num_major_columns).astype(int)
    
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

    fig_width_real_inch *= num_major_columns
    fig_height_real_inch *= nrows
    
    figsize = np.array([fig_width_real_inch, fig_height_real_inch])

    dpi = 92
    
    while True:
        fig_res = figsize * dpi
        if (fig_res > 10000).any():
            figsize = figsize/2
        else:
            break

    # Set view spacing if not specified.
    if view_spacing is None:
        view_spacing = 80//num_views # must be integer
    
    cmap = 'gray'


    width_ratios = width_ratios * num_major_columns
    
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, figsize=figsize, gridspec_kw={'width_ratios': width_ratios}, facecolor='black')
    
    aspect_ratios = np.zeros((3,))
    aspect_ratios[0] = zooms[2]/zooms[1]
    aspect_ratios[1] = zooms[2]/zooms[0]
    aspect_ratios[2] = zooms[1]/zooms[0]

    midpoint = np.floor(np.array(img_array.shape)/2).astype(int) # center of image

    # Want to set a common color map for the whole image.
    img_vmin = img_array.min()
    img_vmax = img_array.max()
    
    for view_index in range(num_views):
        plot_array = img_array
        if mask_path:
            mask_plot_array = mask_array

            
        viewpoint = midpoint - view_spacing*(num_views - 1)/2 + view_spacing*view_index # center of current orhogonal view
        viewpoint = viewpoint.astype(int)

        row_index = view_index % nrows
        start_col = view_index // nrows * 3
        try: # viewpoint could possibly extend outside array; if so, simply plot nothing here.
            axes[row_index, start_col + 0].imshow(np.rot90(plot_array[viewpoint[0],:,:]), cmap=cmap, aspect=aspect_ratios[0], vmin=img_vmin, vmax=img_vmax)
            axes[row_index, start_col + 1].imshow(np.rot90(plot_array[:,viewpoint[1],:]), cmap=cmap, aspect=aspect_ratios[1], vmin=img_vmin, vmax=img_vmax)
            axes[row_index, start_col + 2].imshow(np.rot90(plot_array[:,:,viewpoint[2]]), cmap=cmap, aspect=aspect_ratios[2], vmin=img_vmin, vmax=img_vmax)
        except IndexError:
            pass
        

        if mask_path:
            from matplotlib import cm
            mask_cmap_modified = cm.get_cmap(mask_cmap).copy()
            mask_cmap_modified.set_under('k', alpha=0)
            #axes[row_index, start_col + 0].imshow(np.rot90(mask_plot_array[viewpoint[0],:,:]), cmap=mask_cmap_modified, aspect=aspect_ratios[0], alpha=mask_alpha, vmin=mask_vmin)
            #axes[row_index, start_col + 1].imshow(np.rot90(mask_plot_array[:,viewpoint[1],:]), cmap=mask_cmap_modified, aspect=aspect_ratios[1], alpha=mask_alpha, vmin=mask_vmin)
            #axes[row_index, start_col + 2].imshow(np.rot90(mask_plot_array[:,:,viewpoint[2]]), cmap=mask_cmap_modified, aspect=aspect_ratios[2], alpha=mask_alpha, vmin=mask_vmin)
            try:
                axes[row_index, start_col + 0].imshow(np.rot90(mask_plot_array[viewpoint[0],:,:]), cmap=mask_cmap_modified, aspect=aspect_ratios[0], alpha=mask_alpha, vmin=mask_vmin)
                axes[row_index, start_col + 1].imshow(np.rot90(mask_plot_array[:,viewpoint[1],:]), cmap=mask_cmap_modified, aspect=aspect_ratios[1], alpha=mask_alpha, vmin=mask_vmin)
                axes[row_index, start_col + 2].imshow(np.rot90(mask_plot_array[:,:,viewpoint[2]]), cmap=mask_cmap_modified, aspect=aspect_ratios[2], alpha=mask_alpha, vmin=mask_vmin)
            except IndexError: # Viewpoint could possibly extend outside array; if so, simply plot nothing here.
                pass
                
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
    parser.add_argument('--mask_path', type=str, help='path to mask NIFTI file; must have same dimensions as IN_PATH; mask is assumed to have same (tranformation matrix?) as main image.')
    parser.add_argument('--mask_cmap', type=str, default='Greens', help='mask overlay color map. Must be the name of a Matplotlib Colormap. See https://matplotlib.org/stable/tutorials/colors/colormaps.html for options.')
    parser.add_argument('--mask_alpha', type=float, default=0.5, help='mask overlay alpha (opacity). Must be value between 0 and 1.')
    parser.add_argument('--mask_binarize', action='store_true', help='binarize mask (set values > 0 to 1; values < 0 to 0.')
    parser.add_argument('--mask_vmin', type=float, help='minimum value to display for mask; values below with will have zero opacity.', default=0.1)

    parser.add_argument('--num_views', type=int, help='number of points along a diagonal line on which orthogonal (3-axis) views will be centered', default=6)
    parser.add_argument('--num_major_columns', type=int, help='number of columns of orthogonal views', default=2)
    parser.add_argument('--view_spacing', type=int, help='number of voxels (in all 3 planes) separating points on which orthogonal (3-axis) views will be centered); if none specified, use 80/num_views')

    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))


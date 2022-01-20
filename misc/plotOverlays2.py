#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import nibabel as nib
from matplotlib import pyplot as plt

class Image(object):
    def __init__(self, image_path, image_type):
        self.image_path = image_path
        self.image_type = image_type
        self.data = nib.load(image_path).get_fdata()
            
        

def plotOverlays(main_path, labels_paths, display_mode, num_cuts, out_dir, out_path, alpha):
    ## Load main image.
    main = Image(main_path, 'main')
    
    #### Initialize figure.
    ## Set subplot grid dimensions.
    if display_mode == 'ortho':
        subplot_shape = (1,3)
    elif display_mode in ['x', 'y', 'z']:
        if (not num_cuts is None):
            subplot_shape = (1, num_cuts)
        elif display_mode == 'x':
            subplot_shape = (1, main.data.shape[0])
        elif display_mode == 'y':
            subplot_shape = (1, main.data.shape[1])
        elif display_mode == 'z':
            subplot_shape = (1, main.data.shape[2])
    elif display_mode == 'xyz':
        if (not num_cuts is None):
            subplot_shape = (3, num_cuts)
        else:
            subplot_shape = (3, max(main.data.shape))
    else:
        raise Exception('Invalid display_mode')
        
    ## Set the specific slices to show.
    slices = [[], [], []] # (x slices, y slices, z slices)
    if display_mode == 'ortho':
        slices[0] = np.round(main.data.shape[0]/2)
        slices[1] = np.round(main.data.shape[1]/2)
        slices[2] = np.round(main.data.shape[2]/2)
    elif display_mode in ['x', 'y', 'z']:
        if display_mode == 'x':
            if num_cuts is None:
                slices[0] = range(main.data.shape[2])
        elif display_mode == 'y':
            if num_cuts is None:
                slices[1] = range(main.data.shape[2])
        elif display_mode == 'z':
            if num_cuts is None:
                slices[2] = range(main.data.shape[2])
            else:
                slices[2] = np.linspace(0, main.data.shape[2]-1, num_cuts, dtype=int)
    elif display_mode == 'xyz':
        pass
    
    fig, axes = plt.subplots(subplot_shape[0], subplot_shape[1])#, sharex=True, sharey=True)
    if not hasattr(axes, 'shape'):
        axes = np.array([[axes]])
    elif len(axes.shape) == 1:
        if subplot_shape == (1,1):
            axes = np.array([[axes]])
        if subplot_shape[0] == 1:
            axes = np.expand_dims(axes, axis=0)
        elif subplot_shape[1] == 1:
            axes = np.expand_dims(axes, axis=1)

    if display_mode == 'ortho':
        pass
    elif display_mode in ['x', 'y', 'z']:
        if display_mode == 'x':
            pass
        elif display_mode == 'y':
            pass
        elif display_mode == 'z':
            for row in range(axes.shape[0]):
                for col in range(axes.shape[1]):
                    ax = axes[row, col]
                    ax.imshow(main.data[:,:,slices[2][col]])
    elif display_mode == 'xyz':
        pass
    
    for row in range(axes.shape[0]):
        for col in range(axes.shape[1]):            
            ax = axes[row, col]
            
            ## Remove axis ticks and labels.
            ax.set_xticks([])
            ax.set_yticks([])
    out_name = os.path.basename(main_path).split('.nii')[0]+"-with_overlay.png"
    out_path = os.path.join(out_dir, out_name)
    fig.savefig(out_path)
    #import code; code.interact(local=locals())
    sys.exit()
    for overlay_path in overlay_paths:
        display = plt.plot_img(target_path, display_mode=display_mode, cmap='gray', cut_coords=num_cuts)
        #display = plt.plot_img(target_path, cmap='gray', display_mode='z', cut_coords=10)
        #display.add_overlay(overlay_path, cmap='YlOrRd_r', alpha=alpha, colorbar=True, threshold=1.0)
        display.add_overlay(overlay_path, cmap='winter', alpha=alpha, threshold=0.5)#cmap=plt.cm.purple_green)
        
        if out_path is None:
            out_name = os.path.basename(overlay_path).split('.nii')[0]+"-with_overlay.png"
            if out_dir is None:
                out_dir = os.path.dirname(target_path)
            out_path = os.path.join(out_dir, out_name)
            
        display.savefig(out_path)
        display.close()


if __name__ == "__main__":
    description = """Create figures of images overlaid on a main image. For each input overlay image, create a separate figure."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Positional arguments.
    parser.add_argument("main_path", help='path to main image.', type=str)
    
    # Optional arguments
    parser.add_argument("-l", "--labels_paths", help='paths to overlay images.', nargs="*", type=str)
    parser.add_argument("-m", "--display_mode", help="display_mode argument passed to nilearn.plotting.plot_img. Can use x, y, z, xyz, and ortho", type=str, default="xyz")
    parser.add_argument("-n", "--num_cuts", help="number of cuts. Ignored unless display_mode is 'x', 'y', or 'z'", type=int, default=None)
    parser.add_argument("-d", "--out_dir", help="output directory", default=None, type=str)
    parser.add_argument("-o", "--out_path", help="output path", default=None, type=str)
    parser.add_argument("-a", "--alpha", help="Overlay opacity (alpha argument passed to add_overlay function", type=float, default=1)
    args = parser.parse_args()
    
    # Generate figures
    plotOverlays(**vars(args))

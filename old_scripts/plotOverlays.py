#!/usr/bin/env python

import argparse
import os, sys
from nilearn import plotting as plt
from matplotlib import pyplot

def plotOverlays(target_path, overlay_paths, display_mode, num_cuts, out_dir=None, out_path=None, alpha=0.4):

    if (not display_mode in ['x', 'y', 'z']):
        num_cuts = None
    
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
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("target_path", help='path to main image.', type=str)
    parser.add_argument("overlay_paths", help='paths to overlay images.', nargs="+", type=str)
    
    # Parse arguments.
    parser.add_argument("-m", "--display_mode", help="display_mode argument passed to nilearn.plotting.plot_img. Can use x, y, z, ortho, and more", type=str, default="ortho")
    parser.add_argument("-n", "--num_cuts", help="number of cuts. Ignored unless display_mode is 'x', 'y', or 'z'", type=int, default=10)
    parser.add_argument("-d", "--out_dir", help="output directory", default=None, type=str)
    parser.add_argument("-o", "--out_path", help="output path", default=None, type=str)
    parser.add_argument("-a", "--alpha", help="Overlay opacity (alpha argument passed to add_overlay function", type=float, default=1)
    args = parser.parse_args()
    
    # Do stuff.
    #plotRegistrations(args.target_path, args.overlay_paths, args.out_dir, args.display_mode, args.num_cuts)
    plotOverlays(**vars(args))

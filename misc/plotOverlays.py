#!/usr/bin/env python

import argparse
import os, sys
from nilearn import plotting as plt

def plotRegistrations(target_path, overlay_paths, out_dir, display_mode, num_cuts):

    if (not display_mode in ['x', 'y', 'z']):
        num_cuts = None
    
    for overlay_path in overlay_paths:
        display = plt.plot_img(target_path, display_mode=display_mode, cmap='gray', cut_coords=num_cuts)
        #display = plt.plot_img(target_path, cmap='gray', display_mode='z', cut_coords=10)
        display.add_overlay(overlay_path, cmap='YlOrRd_r', alpha=0.40)

        out_name = os.path.basename(overlay_path).split('.nii')[0]+"_on_main.png"
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
    parser.add_argument("-o", "--out_dir", help="output directory", default=os.getcwd(), type=str)
    parser.add_argument("-d", "--display_mode", help="display_mode argument passed to nilearn.plotting.plot_img. Can use x, y, z, ortho, and more. Default: 'ortho'", type=str, default="ortho")
    parser.add_argument("-n", "--num_cuts", help="number of cuts. Ignored unless display_mode is 'x', 'y', or 'z'. Default=10", type=int, default=10)
    args = parser.parse_args()

    # Do stuff.
    plotRegistrations(args.target_path, args.overlay_paths, args.out_dir, args.display_mode, args.num_cuts)

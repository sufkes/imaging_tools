#!/usr/bin/env python

import argparse
import os, sys
from nilearn import plotting as plt

def plotRegistrations(target_path, overlay_paths, out_dir):

    for overlay_path in overlay_paths:
        display = plt.plot_img(target_path, display_mode='ortho', cmap='gray') # automatically generate cuts (number of cuts = cut_coods)
        display.add_overlay(overlay_path, cmap='YlOrRd_r', alpha=0.40)

        out_name = os.path.basename(overlay_path).split('.nii.gz')[0]+"_on_main.png"
        out_path = os.path.join(out_dir, out_name)
        display.savefig(out_path)
        display.close()


if __name__ == "__main__":
    description = """Create figures of images registered to a target. For each input overlay image, create a separate figure of the overlay image overlaid on the target iamge."""
    epilog = '' # """Text to follow argument explantion """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("target_path", help='path to main image.', type=str)
    parser.add_argument("overlay_paths", help='paths to overlay images.', nargs="+", type=str, metavar=("OVERLAY_PATH_1", "OVERLAY_PATH_2"))
        
    # Parse arguments.
    parser.add_argument("-o", "--out_dir", help="output directory", default=os.getcwd(), type=str)
    args = parser.parse_args()

    # Do stuff.
    plotRegistrations(args.target_path, args.overlay_paths, args.out_dir)

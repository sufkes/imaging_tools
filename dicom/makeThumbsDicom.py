#!/usr/bin/env python

import os, sys
import glob
import argparse
import matplotlib.pyplot as plt
import pydicom

from pprint import pprint


def makeThumbs(in_dirs, out_dir):
    for study_dir in in_dirs:
        series_dirs = sorted(glob.glob(study_dir+"/*"))
        study_name = os.path.basename(study_dir)

#        for series_dir in series_dirs:
        columns = 10
        width = 2*columns
        height = 2*len(series_dirs)/columns
        fig = plt.figure(figsize=(width, height))
        out_name = "{study_name}".format(study_name=study_name)
        out_path = os.path.join(out_dir, out_name)
        for ii, series_dir in enumerate(series_dirs):
            files = sorted(glob.glob(series_dir+"/*"))

            # Pick a file near the "middle" - hopefull will show some brain (or other body part) instead of empty space.
            thumb_file = files[len(files)/2] # floored integer division
            series_name = os.path.basename(series_dir)

            try: # might not have PixelData
                dataset = pydicom.dcmread(thumb_file)
#                out_name = "{study_name}-Series-{series_name}".format(study_name=study_name, series_name=series_name)
#                out_path = os.path.join(out_dir, out_name)
                
#                fig = plt.imshow(dataset.pixel_array, cmap=plt.cm.bone)
#                plt.axis("off")
#                plt.savefig(out_path)
#                plt.close()


                

                sp = plt.subplot(len(series_dirs) / columns + 1, columns, ii + 1)
#                sp.title.set_text(series_name)
                sp.set_title(series_name[:30], fontsize=8)
                plt.imshow(dataset.pixel_array, cmap=plt.cm.bone)
                plt.axis("off")
                
            except:
#                print "Something went wrong for file: '{thumb_file}'".format(thumb_file=thumb_file)
                pass
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close()

if (__name__ == '__main__'):
    # Create argument parser
    description = """Input a directory containing DICOM Studies, generate thumbnails for each Series."""
    parser = argparse.ArgumentParser(description=description)
    
    # Define positional arguments.
    parser.add_argument("in_dirs", help="paths to directories containing DICOM Studies. DICOM Series should be stored in subfolders within DICOM Study folders.", nargs="+", type=str)
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_dir", help="path to directory to which thumbnails will be saved", type=str, default=os.getcwd())

    # Parse arguments.
    args = parser.parse_args()

    # Run the function.
    makeThumbs(args.in_dirs, args.out_dir)

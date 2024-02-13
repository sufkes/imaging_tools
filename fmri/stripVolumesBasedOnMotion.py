#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import nibabel as nib

def main(rel_rms_motion_path, in_image_path, out_image_path, out_volume_report_path, min_vols, max_motion_for_extension, verbose):
    # Read motion file.
    motion = np.loadtxt(rel_rms_motion_path)

    num_edges = len(motion) # edge = between two volumes
    num_vols = len(motion) + 1 # number of motion values between consecutive volumes = number of volumes - 1 

    min_edges = min_vols - 1

    ## Determine which volumes to strip.
    # Determine first volume in motion-minimizing window of minimum size.
    min_subseries_motion_mean = np.inf
    for edge_start in range(0, num_edges - min_edges + 1):
        edge_stop = edge_start + min_edges
        subseries_motion = motion[edge_start:edge_stop]
        subseries_motion_mean = subseries_motion.mean()
        if subseries_motion_mean < min_subseries_motion_mean:
            min_subseries_motion_mean = subseries_motion_mean
            best_edge_start = edge_start
            best_edge_stop = edge_stop
    window_edge_start = best_edge_start
    window_edge_stop = best_edge_stop
    window_vol_start = window_edge_start
    window_vol_stop = window_edge_stop + 1 

    if verbose: print(f'Motion-minimizing window of minimum size: [{window_vol_start}, {window_vol_stop-1}] ({window_vol_stop-window_vol_start} volumes)')

    # Extend backward in time.
    final_edge_start = window_edge_start
    for edge_start in range(window_edge_start - 1, -1, -1):
        if motion[edge_start] <= max_motion_for_extension:
            final_edge_start = edge_start
        else:
            break

    # Extend forward in time.
    final_edge_stop = window_edge_stop
    for edge_stop in range(window_edge_stop + 1, num_edges + 1):
        last_edge = edge_stop - 1
        if motion[last_edge] <= max_motion_for_extension:
            final_edge_stop = edge_stop
        else:
            break

    final_vol_start = final_edge_start
    final_vol_stop = final_edge_stop + 1

    if verbose: print(f'Final volumes after stripping: [{final_vol_start}, {final_vol_stop-1}] ({final_vol_stop-final_vol_start} volumes)')

    ## Strip volumes from NIFTI and save. 
    if not in_image_path is None:
        img_nii = nib.load(in_image_path)
        img_arr = img_nii.get_fdata()
        in_dtype = img_arr.dtype
        #import code; code.interact(local=locals())

        # Strip volumes.
        stripped_img_arr = img_arr[:, :, :, final_vol_start:final_vol_stop]
        if verbose: print(f'Final NIFTI shape: {stripped_img_arr.shape}')

        # Save stripped NIFTI.
        stripped_img_nii = nib.nifti1.Nifti1Image(stripped_img_arr, img_nii.affine)#, img_nii.header)
        nib.save(stripped_img_nii, out_image_path)
        if verbose: print(f'Saved volume-stripped image: {out_image_path}')

        # Save a record of which volumes were included in the striped image.
        if out_volume_report_path is None:
            if out_image_path.endswith('.nii'):
                out_volume_report_path = out_image_path[:-4] + '-included_volumes.txt'
            elif out_image_path.endswith('.nii.gz'):
                out_volume_report_path = out_image_path[:-7] + '-included_volumes.txt'
        report_text = f'{final_vol_start},{final_vol_stop-1}\n'
        with open(out_volume_report_path, 'w') as handle:
            handle.write(report_text)
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Remove consecutive volumes from the start and end of an fMRI series that have excessive motion.  User specifies a minimum number of volumes that must be present in the final stripped series. First, the window of this minimum number of volumes within which the mean relative RMS motion between frames is minimized is determined. This motion-minimizing window is then expanded backward in forward in time, adding new volumes provided that the relative RMS motion between the new volume and its already-included neighbor does not exceed a user-specified value. 

# Example: Input fMRI series has 120 volumes. Run with --min_volumes 100 --max_motion_for_extension 0.25
- The window of 100 volumes within which the mean relative RMS motion between consecutive frames is minimized is [5, 104].
- This 100-volume series is extended backward in time: frames 3 and 4 are added. 
- The relative motion between frames 2 and 3 is 0.5 mm (greater than --max_motion_for_extension 0.25), so frame 2 is not added. The extension backward in time is stopped.
- This 102-volume series is extended forward in time: frames 105-110 are added. 
- The relative motion between frames 110 and 111 is 0.4 mm (greater than --max_motion_for_extension 0.25), so frame 111 is not added. The extension forward in time is stopped.
- The final volume stripped series includes volumes [3, 110]. 
"""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
#    parser.add_argument("", help="")
    
    # Define optional arguments.
    parser.add_argument("-m", "--rel_rms_motion_path", help="path to report of wise relative RMS motion between consecutive frames output by the FSL 'mcflirt' function (the file called <prefix>_rel.rms)", required=True, type=str)
    parser.add_argument("-i", "--in_image_path", help="path to input fMRI NIFTI file to volume-strip. If not specified, the script will print the range of volumes to include based on the motion data.", type=str)
    parser.add_argument("-o", "--out_image_path", help="path to output volume-stripped NIFTI file.", type=str)
    parser.add_argument("-r", "--out_volume_report_path", help="path to output text file which will record the included volumes in the stripped series.", type=str)
    parser.add_argument("-n", "--min_vols", help="minimum number of consecutive volumes to include in volume-stripped fMRI series", type=int, default=100)
    parser.add_argument("-d", "--max_motion_for_extension", help="the maximum relative RMS motion allowed between a new volume and its already-included neighbor for it to be included in the volume-stripped series, in mm", type=float, default=0.25)
    parser.add_argument("-v", "--verbose", action="store_true", help="print information")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

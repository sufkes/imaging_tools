#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import nibabel as nib

from dipy.core.gradients import gradient_table
from dipy.data import get_fnames, default_sphere
from dipy.direction import peaks_from_model
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti, load_nifti_data, save_nifti
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import save_trk
from dipy.reconst.csdeconv import auto_response_ssst
from dipy.reconst.shm import CsaOdfModel
#from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.streamline import Streamlines
import dipy.reconst.dti as dti
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
#from dipy.tracking.stopping_criterion import ActStoppingCriterion
from dipy.tracking._utils import _mapping_to_voxel, _to_voxel_coordinates # for custom exclude_mask and terminal_mask

def getStreamlineLabels(streamline, label_image, affine):
    """Input a streamline (shape Nx3), label image, and affine; return an array (shape N) with the corresponding value in the label image at each point in the streamline"""
    # Determine mapping to voxel coordinates (to determine whether streamline point lies in a given region; only used if certain options are used).
    mapping_to_voxel_linear_transformation, mapping_to_voxel_offset = _mapping_to_voxel(affine)

    # Get the label values at each point in the streamline.
    streamline_voxel_coordinates = _to_voxel_coordinates(streamline, mapping_to_voxel_linear_transformation, mapping_to_voxel_offset)
    i, j, k = streamline_voxel_coordinates.T
    streamline_labels = label_image[i, j, k]
    return streamline_labels

def main(dmri_path, bval_path, bvec_path, out_trk_path, seed_mask_path, tracking_mask_path, exclude_mask_path, terminal_mask_path, wm_trim_mask_path, roi_path, out_network_matrix_path, seed_density, seed_mask_min, tracking_mask_min, exclude_mask_min, terminal_mask_min, wm_trim_mask_min):
    # Load dMRI image and bvecs/bvals.
    data, affine, dmri_img = load_nifti(dmri_path, return_img=True)
    bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
    gtab = gradient_table(bvals, bvecs)
    
    # Load masks.
    seed_mask_data = load_nifti_data(seed_mask_path)
    seed_mask = seed_mask_data > seed_mask_min

    tracking_mask_data = load_nifti_data(tracking_mask_path)
    tracking_mask = tracking_mask_data > tracking_mask_min

    # Define stopping criterion.
    #stopping_criterion = ActStoppingCriterion(include_map, exclude_map)
    stopping_criterion = BinaryStoppingCriterion(tracking_mask)

    ## Step 1: Getting directions from a diffusion dataset
    tenmodel = dti.TensorModel(gtab)
    csa_peaks = peaks_from_model(tenmodel, data, default_sphere, relative_peak_threshold=0.8, min_separation_angle=45, mask=tracking_mask)

    ## Step 3: Defining a set of seeds from which to begin tracking
    seeds = utils.seeds_from_mask(seed_mask, affine, density=seed_density)
    
    # Initialization of LocalTracking. The computation happens in the next step.
    print('Calculating streamlines')
    streamlines_generator = LocalTracking(csa_peaks, stopping_criterion, seeds, affine=affine, step_size=0.5)

    # Generate streamlines object.
    streamlines = Streamlines(streamlines_generator)

    ##### Remove and modify streamlines based on options specified.
    ## Remove streamlines of length 1 or less, since the connectivity matrix calculator will throw an error.
    streamlines = [streamline for streamline in streamlines if len(streamline)>1]

    #save_trk(StatefulTractogram(streamlines, dmri_img, Space.RASMM), out_trk_path.replace('.trk', '-before_exclusions.trk'), streamlines)

    save_trk(StatefulTractogram(streamlines, dmri_img, Space.RASMM), out_trk_path.replace('.trk','-before_exclusions.trk'), streamlines)
    
    # Store excluded streamlines for testing.
    streamlines_excluded = []
    
    ## Remove streamlines which enter exclude mask.
    if not exclude_mask_path is None:
        # Load exclude mask
        exclude_mask_data = load_nifti_data(exclude_mask_path)
        exclude_mask = exclude_mask_data > exclude_mask_min

        # Check whether streamline enters exlude mask, and if it does, exclude it.
        print(f'Numer of streamlines before removing streamlines that enter the exclude mask: {len(streamlines)}')
        streamlines_after_exclusion = []
        for streamline in streamlines:
            streamline_exclude_mask_values = getStreamlineLabels(streamline, exclude_mask, affine)
            if (streamline_exclude_mask_values == False).all():
                streamlines_after_exclusion.append(streamline)
            else:
                streamlines_excluded.append(streamline)
        streamlines = streamlines_after_exclusion
        print(f'Numer of streamlines after removing streamlines that enter the exclude mask:  {len(streamlines)}')
        #save_trk(StatefulTractogram(streamlines, dmri_img, Space.RASMM), out_trk_path.replace('.trk', '-after_exclude_mask.trk'), streamlines)

    ## Remove streamlines which never leave the WM trim mask.
    if not wm_trim_mask_path is None:
        # Load WM trim mask.
        wm_trim_mask_data = load_nifti_data(wm_trim_mask_path)
        wm_trim_mask = wm_trim_mask_data > wm_trim_mask_min

        # Remove streamlines which are only in the WM trim mask.
        print(f'Numer of streamlines before removing WM-only streamlines: {len(streamlines)}')
        streamlines_after_exclusion = []
        for streamline in streamlines:
            streamline_wm_trim_mask_values = getStreamlineLabels(streamline, wm_trim_mask, affine)
            if (streamline_wm_trim_mask_values == False).any():
                streamlines_after_exclusion.append(streamline)
            else:
                streamlines_excluded.append(streamline)
        streamlines = streamlines_after_exclusion
        print(f'Numer of streamlines after removing WM-only streamlines:  {len(streamlines)}')
        #save_trk(StatefulTractogram(streamlines, dmri_img, Space.RASMM), out_trk_path.replace('.trk', '-after_wm_trim_mask.trk'), streamlines)

    # Remove streamlines which do not terminate in terminal_mask.
    if not terminal_mask_path is None:
        # Load terminal mask.
        terminal_mask_data = load_nifti_data(terminal_mask_path)
        terminal_mask = terminal_mask_data > terminal_mask_min

        # Remove streamlines which do not terminate in the terminal mask
        print(f'Numer of streamlines before removing streamlines not terminating in terminal mask: {len(streamlines)}')
        streamlines_after_exclusion = []
        for streamline in streamlines:
            streamline_terminal_mask_values = getStreamlineLabels(streamline, terminal_mask, affine)
            if (streamline_terminal_mask_values[0] == True) and (streamline_terminal_mask_values[-1] == True):
                streamlines_after_exclusion.append(streamline)
            else:
                streamlines_excluded.append(streamline)
        streamlines = streamlines_after_exclusion
        print(f'Numer of streamlines after removing streamlines not terminating in terminal mask:  {len(streamlines)}')
        #save_trk(StatefulTractogram(streamlines, dmri_img, Space.RASMM), out_trk_path.replace('.trk', '-after_terminal_mask.trk'), streamlines)

    # Save TRK file for current ROI.
    print('Saving tractogram')
    sft = StatefulTractogram(streamlines, dmri_img, Space.RASMM)
    save_trk(sft, out_trk_path, streamlines)

    # Save excluded streamlines to TRK file for testing.
    save_trk(StatefulTractogram(streamlines_excluded, dmri_img, Space.RASMM), out_trk_path.replace('.trk','-excluded_streamlines.trk'), streamlines)
    
    # Generate ROI x ROI connectivity matrix.
    if not roi_path is None:
        roi_img = load_nifti_data(roi_path).astype(int)
        #network_matrix, network_mapping = utils.connectivity_matrix(streamlines, affine, roi_img, inclusive=True, symmetric=True, return_mapping=True, mapping_as_streamlines=False)
        network_matrix = utils.connectivity_matrix(streamlines, affine, roi_img, inclusive=False, symmetric=True, return_mapping=False, mapping_as_streamlines=False)
        network_matrix = network_matrix[1:, 1:] # remove connectivity to ROI label 0
        np.savetxt(out_network_matrix_path, network_matrix)
        #np.savetxt(out_network_mapping_path, network_mapping) # not working and I don't care.
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Run deterministic tractography, with seeds and streamlines restricted to a single binary mask. Output a .trk file, and an ROI x ROI connectivity matrix (optional)."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)

    # Define required arguments.
    parser.add_argument("-i", "--dmri_path", help="path to 4D diffusion MRI NIFTI file", required=True)
    parser.add_argument("-b", "--bval_path", help="path to BVAL file, in format output by dcm2niix", required=True)
    parser.add_argument("-g", "--bvec_path", help="path to BVEC file, in format output by dcm2niix", required=True)
    parser.add_argument("-o", "--out_trk_path", help="path to output .trk file", required=True)
    parser.add_argument("-s", "--seed_mask_path", help="path to NIFTI file which defines the seed voxels (e.g. a WM mask or partial volume estimate)", required=True)
    parser.add_argument("-t", "--tracking_mask_path", help="path to NIFTI file which defines the voxels in which tracking can occur (e.g. a WM mask or partial volume estimate)", required=True)
    
    # Define optional arguments.
    parser.add_argument("--exclude_mask_path", help="if specified, streamlines which enter this mask will be excluded (e.g. a CSF mask or partial volume estimate)")
    parser.add_argument("--terminal_mask_path", help="if specified, retain only streamlines which terminate within this mask (e.g. a mask or partial volume estimate including the grey matter and background (for brainstem tracts)")
    #parser.add_argument("--wm_trim_mask_path", help="if specified, the ends of streamlines lying in this white matter mask will be trimmed from the streamline (e.g. WM-GM-WM-GM would be trimmed to GM-WM-GM); then, streamlines which never enter this white matter mask will be removed (e.g. WM-only streamlines would be removed; WM-GM-WM would be trimmed and then removed). The purpose of this is to remove streamlines which never enter the WM, and to trim streamlines such that they terminate outside of the WM. Can input a binary mask or partial volume estimate. Note: streamlines that enter and leave the WM mask multiple times will still be allowed by this option (e.g. GM-WM-GM-WM-GM would not be trimmed or removed).")
    parser.add_argument("--wm_trim_mask_path", help="if specified, streamlines lying solely within this mask will be removed. This is intended to remove fragmentary streamlines which never enter GM regions.")
    parser.add_argument("-r", "--roi_path", help="path to NIFTI file with integer labels defining regions of interest (ROIs) to be used to generate an ROI x ROI connectivity matrix")
    parser.add_argument("-n", "--out_network_matrix_path", help="path to output connectivity matrix")
    parser.add_argument("--seed_density", help="number of tractography seeds per seed voxel specify either a single integer, or three integers to specify a grid of seed points to be used in each seed voxel (e.g. for a 2x2x2 grid of seeds in each voxel, use '--seed_density 2 2 2').", default=1, nargs="+", type=int)
    parser.add_argument("--seed_mask_min", help='values greater than SEED_MASK_MIN in SEED_MASK_PATH will be used to define the seed mask', default=0.0, type=float)
    parser.add_argument("--tracking_mask_min", help='values greater than TRACKING_MASK_MIN in TRACKING_MASK_PATH will be used to define the tracking mask', default=0.0, type=float)
    parser.add_argument("--exclude_mask_min", help="values greater than EXCLUDE_MASK_MIN in EXCLUDE_MASK_PATH will be used to define the exclusion mask", default=0.0, type=float)
    parser.add_argument("--terminal_mask_min", help="values greater than TERMINAL_MASK_MIN in TERMINAL_MASK_PATH will be used to define the terminal mask", default=0.0, type=float)
    parser.add_argument("--wm_trim_mask_min", help="values greater than WM_TRIM_MASK_MIN in WM_TRIM_MASK_PATH will be used to define the white matter trim mask", default=0.0, type=float)
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

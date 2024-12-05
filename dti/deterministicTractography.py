#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np

from dipy.core.gradients import gradient_table
from dipy.data import get_fnames, default_sphere
from dipy.direction import peaks_from_model
from dipy.io.gradients import read_bvals_bvecs
from dipy.io.image import load_nifti, load_nifti_data
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.io.streamline import save_trk
from dipy.reconst.csdeconv import auto_response_ssst
from dipy.reconst.shm import CsaOdfModel
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.streamline import Streamlines
import dipy.reconst.dti as dti
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion

def main(dmri_path, bval_path, bvec_path, out_trk_path, seed_path, tracking_boundary_path, roi_path, out_network_matrix_path, out_network_mapping_path, seed_lower_threshold, tracking_boundary_lower_threshold):
    # Load dMRI image and bvecs/bvals.
    data, affine, dmri_img = load_nifti(dmri_path, return_img=True)
    bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
    gtab = gradient_table(bvals, bvecs)

    # Load seed and tracking boundary images.
    seed_img = load_nifti_data(seed_path)
    tracking_boundary_img = load_nifti_data(tracking_boundary_path)

    # Apply thresholds to get masks.
    seed_mask = seed_img >= seed_lower_threshold
    tracking_boundary_mask = tracking_boundary_img >= tracking_boundary_lower_threshold

    ## Step 1: Getting directions from a diffusion dataset
    tenmodel = dti.TensorModel(gtab)
    csa_peaks = peaks_from_model(tenmodel, data, default_sphere, relative_peak_threshold=.8, min_separation_angle=45, mask=tracking_boundary_mask)

    ## Step 2: Identifying when the tracking must stop
    stopping_criterion = BinaryStoppingCriterion(tracking_boundary_mask)

    ## Step 3: Defining a set of seeds from which to begin tracking
    seeds = utils.seeds_from_mask(seed_mask, affine, density=[1, 1, 2])
    
    # Initialization of LocalTracking. The computation happens in the next step.
    streamlines_generator = LocalTracking(csa_peaks, stopping_criterion, seeds, affine=affine, step_size=.5)

    # Generate streamlines object
    streamlines = Streamlines(streamlines_generator)

    # Save TRK file.
    sft = StatefulTractogram(streamlines, dmri_img, Space.RASMM)
    save_trk(sft, out_trk_path, streamlines)

    # Generate ROI x ROI connectivity matrix.
    if not roi_path is None:
        roi_img = load_nifti_data(roi_path).astype(int)
        network_matrix, network_mapping = utils.connectivity_matrix(streamlines, affine, roi_img, inclusive=True, symmetric=True, return_mapping=True, mapping_as_streamlines=False)
        network_matrix=network_matrix[1:, 1:]
        np.savetxt(out_network_matrix_path, network_matrix)
        np.savetxt(out_network_mapping_path, network_mapping)
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = """Run deterministic tractography. Output a .trk file, and a ROI x ROI connectivity matrix (optional)."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("-i", "--dmri_path", help="path to 4D diffusion MRI NIFTI file", required=True)
    parser.add_argument("-b", "--bval_path", help="path to BVAL file, in format output by dcm2niix", required=True)
    parser.add_argument("-g", "--bvec_path", help="path to BVEC file, in format output by dcm2niix", required=True)

    parser.add_argument("-o", "--out_trk_path", help="path to output .trk file", required=True)

    parser.add_argument("-s", "--seed_path", help="path to NIFTI file which defines the seed voxels", required=True)
    parser.add_argument("-m", "--tracking_boundary_path", help="path to NIFTI file which defines the region within which fibre tracking occurs; tracking stops when it leads outside of this mask", required=True)

    # Define optional arguments.
    parser.add_argument("-r", "--roi_path", help="path to NIFTI file with integer labels defining regions of interest (ROIs) to be used to generate an ROI x ROI connectivity matrix")
    parser.add_argument("-n", "--out_network_matrix_path", help="path to output connectivity matrix")
    parser.add_argument("-p", "--out_network_mapping_path", help="path to network mapping (MATRIX?) which maps network matrix maps matrix indices to streamlines")
    parser.add_argument("--seed_lower_threshold", type=float, default=0.1, help="voxels above this value in the the seed image will be used as seeds")
    parser.add_argument("--tracking_boundary_lower_threshold", type=float, default=0.1, help="tracking will occur in voxels above this value in the the tracking image")

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

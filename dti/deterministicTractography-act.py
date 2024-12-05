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
from dipy.tracking.stopping_criterion import ThresholdStoppingCriterion
from dipy.tracking import utils
from dipy.tracking.local_tracking import LocalTracking
from dipy.tracking.streamline import Streamlines
import dipy.reconst.dti as dti
from dipy.tracking.stopping_criterion import BinaryStoppingCriterion
from dipy.tracking.stopping_criterion import ActStoppingCriterion

def main(dmri_path, bval_path, bvec_path, out_trk_path, roi_path, out_network_matrix_path, gm_path, wm_path, csf_path, seed_density):
    # Load dMRI image and bvecs/bvals.
    data, affine, dmri_img = load_nifti(dmri_path, return_img=True)
    bvals, bvecs = read_bvals_bvecs(bval_path, bvec_path)
    gtab = gradient_table(bvals, bvecs)

    # Load GM, WM, CSF probability maps.
    pve_csf_data = load_nifti_data(csf_path)
    pve_gm_data = load_nifti_data(gm_path)
    pve_wm_data = load_nifti_data(wm_path)

    # Generate a backgroud map (where P(WM or GM or CSF) = 0)
    background = np.ones(pve_gm_data.shape)
    background[(pve_gm_data + pve_wm_data + pve_csf_data) > 0] = 0

    # Generate the include_map; streamlines terminating in these regions will be considered valid.
    include_map = pve_gm_data + background
    exclude_map = pve_csf_data

    # Define stopping criterion.
    stopping_criterion = ActStoppingCriterion(include_map, exclude_map)

    # Define seeds.
    seed_mask = pve_wm_data >= 0.5

    ## Step 1: Getting directions from a diffusion dataset
    tenmodel = dti.TensorModel(gtab)
    print('Calculating csa_peaks in WM + GM + CSF (not sure if this is ideal)')
    csa_peaks = peaks_from_model(tenmodel, data, default_sphere, relative_peak_threshold=0.8, min_separation_angle=45, mask=(pve_gm_data + pve_wm_data + pve_csf_data))

    ## Step 3: Defining a set of seeds from which to begin tracking
    seeds = utils.seeds_from_mask(seed_mask, affine, density=seed_density)
    
    # Initialization of LocalTracking. The computation happens in the next step.
    print('Calculating streamlines')
    streamlines_generator = LocalTracking(csa_peaks, stopping_criterion, seeds, affine=affine, step_size=0.5)

    # Generate streamlines object.
    streamlines = Streamlines(streamlines_generator)

    # Remove streamlines of length 1 or less, since the connectivity matrix calculator will throw an error.
    streamlines = [streamline for streamline in streamlines if len(streamline)>1]

    # Save TRK file for current ROI.
    print('Saving tractogram')
    sft = StatefulTractogram(streamlines, dmri_img, Space.RASMM)
    save_trk(sft, out_trk_path, streamlines)
    
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
    description = """Run deterministic tractography using Anatomically Constrained Tractography (ACT) method. Output a .trk file, and an ROI x ROI connectivity matrix (optional)."""
    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument("-i", "--dmri_path", help="path to 4D diffusion MRI NIFTI file", required=True)
    parser.add_argument("-b", "--bval_path", help="path to BVAL file, in format output by dcm2niix", required=True)
    parser.add_argument("-g", "--bvec_path", help="path to BVEC file, in format output by dcm2niix", required=True)
    parser.add_argument("-o", "--out_trk_path", help="path to output .trk file", required=True)

    # Define optional arguments.
    parser.add_argument("-r", "--roi_path", help="path to NIFTI file with integer labels defining regions of interest (ROIs) to be used to generate an ROI x ROI connectivity matrix")
    parser.add_argument("-n", "--out_network_matrix_path", help="path to output .trk file")
    #parser.add_argument("-p", "--out_network_mapping_path", help="path to network mapping (MATRIX?) which maps network matrix maps matrix indices to streamlines")
    parser.add_argument("--gm_path", help='path to grey matter probability map')
    parser.add_argument("--wm_path", help='path to white matter probability map')
    parser.add_argument("--csf_path", help='path to cerebrospinal fluid probability map')
    parser.add_argument("--seed_density", help="number of tractography seeds per seed voxel specify either a single integer, or three integers to specify a grid of seed points to be used in each seed voxel (e.g. for a 2x2x2 grid of seeds in each voxel, use '--seed_density 2 2 2').", default=1, nargs="+", type=int)
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

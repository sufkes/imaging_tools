#!/usr/bin/env python

import os, sys
import argparse

# Hide annoying pydicom warning
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import nibabel as nib
warnings.filterwarnings("default", category=UserWarning)

import numpy as np

def buildConnectivityMatrix(bold_path, mask_path, mask_thres=0.0, corr_dist_min=10.0, no_fisher_r_to_z=False, no_absolute=False):
    
    bold_nii = nib.load(bold_path)
    mask_nii = nib.load(mask_path)
    bold = np.asarray(bold_nii.get_fdata()).astype(np.single)
    mask = np.asarray(mask_nii.get_fdata()).astype(np.single)

    #print("BOLD:", bold.shape, bold.dtype)
    #print("Mask:", mask.shape, mask.dtype)

    zooms = bold_nii.header.get_zooms() # "resolution" (x, y, z, TR)
    #print(zooms)
    num_voxels = bold.shape[0]*bold.shape[1]*bold.shape[2]

    #### Generate an NxN matrix storing the Euclidean distance between each pair of voxels.
    #coords = np.zeros((num_voxels, 3), dtype=int) # Store the coordinate of each voxel.
    #coords[:,0] = range(bold.shape[0] * bold.shape[1]*bold.shape[2]) // bold.shape[1]*bold.shape[2]
    #coords[:,1] = range(np.arange(bold.shape[1] * bold.shape[0]*bold.shape[2]) // bold.shape[0]*bold.shape[2]
    #coords[:,2] = np.arange(bold.shape[2] * bold.shape[0]*bold.shape[1]) // bold.shape[0]*bold.shape[1]

    ## Generate an Nx3 matrix which stores the coordinates for each voxel.
    # This method works fine but is not clean.
    #coords = np.zeros((num_voxels, 3), dtype=int)
    #coords[:, 0] = np.tile(np.repeat(np.arange(bold.shape[0]), bold.shape[1]*bold.shape[2]), 1                          )
    #coords[:, 1] = np.tile(np.repeat(np.arange(bold.shape[1]), bold.shape[2]              ), bold.shape[0]              )
    #coords[:, 2] = np.tile(np.repeat(np.arange(bold.shape[2]), 1                          ), bold.shape[0]*bold.shape[1])
    #print coords.shape

    coords = np.indices((bold.shape[0], bold.shape[1], bold.shape[2])) # the indices for voxel [x, y, z] are coords[:, x, y, z] 
    coords = np.moveaxis(coords, 0, -1) # the indices for voxel [x, y, z] are coords[x, y, z, :] (analogous to storage of bold data; do this so that bold and coordinate map can be reshaped the same way).
    #print coords.shape
    
    #euclid = np.zeros((num_voxels, num_voxels))
    #euclid = np.zeros((bold.shape[0], bold.shape[1], bold.shape[2], bold.shape[0], bold.shape[1], bold.shape[2]))
    #print euclid.shape

#    H, W, D = (3,4,5) # DEBUG: test values only
#    x1, y1, z1, x2, y2, z2 = np.ogrid[:H, :W, :D, :H, :W, :D] # for method 1, 3    
#    x, y, z = np.ogrid[:H, :W, :D] # for method 2

    
    #x1, y1, z1, x2, y2, z2 = np.ogrid[:bold.shape[0], :bold.shape[1], :bold.shape[2], :bold.shape[0], :bold.shape[1], :bold.shape[2]]
    #euclid1 = np.sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )
    
    #x, y, z = np.ogrid[:bold.shape[0], :bold.shape[1], :bold.shape[2]]
    #euclid2 = np.sqrt( (x[:,:,:,None,None,None]-x)**2 + (y[:,:,:,None,None,None]-y)**2 + (z[:,:,:,None,None,None]-z)**2)

    #import numexpr as ne
    #euclid3 = ne.evaluate('sqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )')
    
    #print euclid1.shape
    #print euclid2.shape
    #print euclid3.shape
    #print (euclid1 == euclid2).all()
    #print (euclid1 == euclid3).all()
    

    #print(bold.shape, mask.shape)
    
    # Flatten the BOLD image and mask
    bold_flat = bold.reshape(num_voxels, -1) # automatically set the size along the time-axis
    coords_flat = coords.reshape(num_voxels, -1)
    mask_flat = mask.reshape(num_voxels)
    #print bold_flat.shape
    #print coords_flat.shape
    
    #print(bold_flat.shape, mask_flat.shape)

    # Mask the flattened BOLD image
    nonzero_bold = (bold_flat != 0.0).all(axis=1)
    #print (nonzero_bold.shape)
    print "Number of nonzero BOLD voxels:", np.count_nonzero(nonzero_bold)
    in_mask = mask_flat > mask_thres
    #print (in_mask.shape)
    print "Number of voxels in thresholded mask:", np.count_nonzero(in_mask)
    nonzero_bold_and_in_mask = nonzero_bold & in_mask
    #print (nonzero_bold_and_in_mask.shape)
    print "Number of nonzero BOLD voxels in thresholded mask:", np.count_nonzero(nonzero_bold_and_in_mask)

    #bold_mask = bold_flat[(mask_flat>mask_thres), :] # flattened; only including mask voxels.
    bold_mask = bold_flat[nonzero_bold_and_in_mask, :] # flattened; only include BOLD voxels that are nonzero and inside the thresholded mask.
    coords_mask = coords_flat[nonzero_bold_and_in_mask, :]
    #print bold_mask.shape
    #print coords_mask.shape


    #### Compute the Pearson coorelation coefficient for each ROI pair.
    
    #corr = np.zeros((num_nodes, num_nodes))
    num_nodes = bold_mask.shape[0]
    ### Cao:
    # - A gray matter mask (number of voxels = 7101) was predefined through thresholding the combing cortex and deep gray matter probability templates.
    # - Pearson's correlation between the BOLD time series of each pair of voxels within the gray matter mask was calculated. 
    # - Fisher's r-to-z transformation was applied to improve the normality of the correlation coefficients
    # - Absolute values of all correlations were used to obtain the correlation matrix for each subject.
    # - Connectivity terminating within 10 mm of each source voxel center was set to zero to avoid potential shared signals between nearby voxels.

    corr = np.corrcoef(bold_mask)
    #print corr[:10,:10]
    #print corr.shape

    ## Get the positions in milimeters but multiplying the voxel size along each dimension. Handle non-cubic voxels accordingly
    position = np.zeros(coords_mask.shape, dtype=np.single)
    position[:, 0] = coords_mask[:, 0]*zooms[0]
    position[:, 1] = coords_mask[:, 1]*zooms[1]
    position[:, 2] = coords_mask[:, 2]*zooms[2]

    ## Compute the matrix of Euclidean distances between each point using Einstein summation notation. I don't know exactly how this works. But it's fast!
    # https://stackoverflow.com/questions/43367001/how-to-calculate-euclidean-distance-between-pair-of-rows-of-a-numpy-array/43368088
    position_new_view = position.reshape(position.shape[0], 1, position.shape[1])
    dist = np.sqrt(np.einsum('ijk, ijk->ij', position-position_new_view, position-position_new_view))
    print "\nPairwise distances between voxels:"
    print dist[:10,:10]
    print dist.shape

    ## Set correlations between voxels to zero if they lie too close to each other.
    corr[dist < corr_dist_min] = 0.0
    print "\nCorrelation matrix after zeroing correlations between voxels within "+str(corr_dist_min)+"mm of each other:"
    print corr[:10,:10]

    ## Apply Fisher's r-to-z transformation.
    if not no_fisher_r_to_z:
        corr = np.arctanh(corr)
        print "\nCorrelation matrix after Fisher's r-to-z transformation:"
        print corr[:10,:10]

    ## Take the absolute values of the correlations.
    if not no_absolute:
        corr = np.abs(corr)
        print "\nCorrelation matrix after taking absolute values:"
        print corr[:10,:10]
    
    ## Test modifying the data. Here I want to check that I can get from the list of flattened, masked voxels back to the image space. I test on the bold image.
    # Reconstruct the bold series from the masked image.
    new_bold_flat = np.zeros(bold_flat.shape)
    #new_bold_flat[mask_flat>mask_thres] = bold_mask # add the BOLD voxels that were included in the mask
    new_bold_flat[nonzero_bold_and_in_mask] = bold_mask # add the BOLD voxels that were included in the mask

    #new_mask = mask_flat.reshape(mask.shape)
    new_bold = new_bold_flat.reshape(bold.shape)
    #new_mask_nii = nib.Nifti1Image(new_mask, mask_nii.affine, mask_nii.header) # create NIFTI image with the new data, old header, and old affine.
    new_bold_nii = nib.Nifti1Image(new_bold, bold_nii.affine, bold_nii.header) # create NIFTI image with the new data, old header, and old affine.    
    
    # Test saving a modified NIFTI file.
    #new_mask_path = os.path.join(os.path.dirname(mask_path), "modified-"+os.path.basename(mask_path))
    #print("Saving to", new_mask_path)
    #nib.save(new_mask_nii, new_mask_path)

    new_bold_path = os.path.join(os.path.dirname(bold_path), "modified-"+os.path.basename(bold_path))
    print("Saving to", new_bold_path)
    nib.save(new_bold_nii, new_bold_path)
    
    


if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate a voxel-wise functional connectivity matrix from fMRI data."""
    epilog = '' # """Text to follow argument explantion """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("bold_path", type=str, help="path to BOLD time series in NIFTI format.")
    parser.add_argument("mask_path", type=str, help="path to mask in NIFTI format.")
    
    # Define optional arguments.
    parser.add_argument("-t", "--mask_thres", type=float, help="include only voxels for which the mask has value greater than MASK_THRES. Default: 0", default=0.0)
    parser.add_argument("-d", "--corr_dist_min", type=float, help="set the correlation between pairs of voxels to zero if their centroids are within CORR_DIST_MIN of each other; enter in mm. Default=10", default=10.0)
    parser.add_argument("-z", "--no_fisher_r_to_z", action="store_true", help="do not apply Fisher's r-to-z transformation.")
    parser.add_argument("-a", "--no_absolute", action="store_true", help="do not take the absolute values of the correlations.")

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Do stuff.
    buildConnectivityMatrix(**vars(args))
#    buildConnectivityMatrix(bold_path=args.bold_path, mask_path=args.mask_path, mask_thres=args.mask_thres, corr_dist_min=args.corr_dist_min, no_fisher_r_to_z=args.no_fisher_r_to_z, no_absolute=args.no_absolute)

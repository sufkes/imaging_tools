#!/usr/bin/env python

import os, sys
import argparse

# Hide annoying pydicom warning
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import nibabel as nib
warnings.filterwarnings("default", category=UserWarning)

import numpy as np

def flatMaskedToImageSpace(flat_masked_array, flat_mask, spatial_image_shape, measure_dim=0):
    if (measure_dim == 0): # if the measure being mapped back to the image space is a scalar (i.e. a single number)
        flat_image_shape = (spatial_image_shape[0]*spatial_image_shape[1]*spatial_image_shape[2], )
        final_image_shape = spatial_image_shape
    elif (measure_dim == 1): # if the measure being mapped back to the image space is a 1D vector.
        flat_image_shape = (spatial_image_shape[0]*spatial_image_shape[1]*spatial_image_shape[2], flat_masked_array.shape[1])
        final_image_shape = spatial_image_shape + (flat_masked_array.shape[1],)
    else:
        message = "Cannot map measure back to image space - measure must be 0-dimensional (scalar) or 1-dimensional (vector)."
        raise Exception(message)

    # Initialize the flat matrix which has a row for every voxel in the image space.
    measure_image_space = np.zeros((flat_image_shape), dtype=np.float32)
        
    # Populate the voxels in the flat array with the values taken from the flat, masked array.
    #print measure_image_space.shape, flat_mask.shape, flat_masked_array.shape
    
    measure_image_space[flat_mask] = flat_masked_array
    
    # Reshape to the original image space.
    measure_image_space = measure_image_space.reshape(final_image_shape)

    return measure_image_space

def buildConnectivityMatrix(bold_path, mask_path=None, out_dir=None, mask_thres=0.0, corr_dist_min=10.0, fisher_r_to_z=False, no_absolute=False, network_density=0.05, debug=False, voxel_map_only=False, fcs_only=False, **kwargs):

    if (out_dir is None):
        out_dir = os.path.dirname(bold_path) # save it wherever the BOLD image is
    
    bold_nii = nib.load(bold_path)
    mask_nii = nib.load(mask_path)
    bold = np.asarray(bold_nii.get_fdata()).astype(np.float32)
    if (mask_path is None):
        mask = np.ones(bold.shape, dtype=np.float32)
    else:
        mask = np.asarray(mask_nii.get_fdata()).astype(np.float32)

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

    #coords = np.indices((bold.shape[0], bold.shape[1], bold.shape[2])) 
    coords = np.indices(bold.shape[:3]) # the indices for voxel [x, y, z] are coords[:, x, y, z]
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

    # Special option to fix a problem.
    if voxel_map_only:
        # "Voxel" map, which is required to map from the voxels in the map back to the voxels in the image. Depends on the BOLD image, mask, and mask threshold.
        voxel_map_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_voxelmap.npy")
        print "Saving voxel map to:", voxel_map_path
        np.save(voxel_map_path, nonzero_bold_and_in_mask)
        return
    
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

    corr = np.corrcoef(bold_mask).astype(np.float32) # numpy.corrcoef returns array of 64-bit floats; apparently no way to use 32-bit, but doing so might cause problems anyway.
    if debug:
        print "Raw correlation matrix:"
        print corr[:10,:10]
        print corr.shape
        print

    ## Get the positions in milimeters but multiplying the voxel size along each dimension. Handle non-cubic voxels accordingly
    position = np.zeros(coords_mask.shape, dtype=np.single)
    position[:, 0] = coords_mask[:, 0]*zooms[0]
    position[:, 1] = coords_mask[:, 1]*zooms[1]
    position[:, 2] = coords_mask[:, 2]*zooms[2]

    ## Compute the matrix of Euclidean distances between each point using Einstein summation notation. I don't know exactly how this works. But it's fast!
    # https://stackoverflow.com/questions/43367001/how-to-calculate-euclidean-distance-between-pair-of-rows-of-a-numpy-array/43368088
    position_new_view = position.reshape(position.shape[0], 1, position.shape[1])
    dist = np.sqrt(np.einsum('ijk, ijk->ij', position-position_new_view, position-position_new_view))
    if debug:
        print "Pairwise distances between voxels:"
        print dist[:10,:10]
        print dist.shape
        print

    ## Set correlations between voxels to zero if they lie too close to each other.
    corr[dist < corr_dist_min] = 0.0
    if debug:
        print "Correlation matrix after zeroing correlations between voxels within "+str(corr_dist_min)+"mm of each other:"
        print corr[:10,:10]
        print corr.shape
        print

    ## Take the absolute values of the correlations.
    # From Wang 2015 - GRETNA:
    # """Previous R-fMRI studies have found that certain functional systems are anti-correlated (i.e., have a negative correlation) in their spontaneous brain activity (Greicius et al., 2003; Fox et al., 2005). However, negative correlations may also be introduced by global signal removal, a preprocessing step that is currently controversial (Fox et al., 2009; Murphy et al., 2009; Weissenbacher et al., 2009; Scholvinck et al., 2010). For network topology, negative correlations may have detrimental effects on TRT reliability (Wang et al., 2011) and exhibit organizations different from positive correlations (Schwarz and McGonigle, 2011). Accordingly, GRETNA provides options for researchers to determine the network connectivity members, based on which subsequent graph analyses are implemented: positive network (composed of only positive correlations), negative network (composed of only absolute negative correlations) or full network (composed of both positive correlations and the absolute values of the negative correlations)."""
    if (not no_absolute):
        corr = np.abs(corr)
        if debug:
            print "Correlation matrix after taking absolute values:"
            print corr[:10,:10]
            print corr.shape
            print

    ## Apply Fisher's r-to-z transformation.
    corr_r = corr
    corr_z = np.arctanh(corr)
    if debug:
        print "Correlation matrix after Fisher's r-to-z transformation:"
        print corr_z[:10,:10]
        print corr_z.shape
        print
    if fisher_r_to_z: # If fisher_r_to_z is True, the transformed correlation matrix will be used construct the adjacency matrix, otherwise the z-transformed result will just be saved separately.
        corr = corr_z
    
    #### Calculate the network adjaceny matrices by thresholding the correlation matrix.
    ## Matrix requirements for the Brain Connectivity Toolbox (BCT), from their website:
    # """ 
    # Which network matrices can I use with the Brain Connectivity Toolbox? 
    # - Most functions do not explicitly check the validity of the input network matrices; it is crucial to manually ensure that these matrices are suitable for their intended use.
    # - The network matrices should be square; rows and columns in these matrices should represent network nodes, matrix entries should represent network links.
    # - The network matrices should not be too small. As a rule of thumb, the toolbox is designed to be used with networks of greater than 20 nodes.
    # - The network matrices should preferably be in double-precision and non-sparse formats. Sparse, single-precision or logical formats may sometimes cause errors.
    # - The network matrices may be binary or weighted, directed or undirected. Each function specifies the network type for which it is suitable.
    # - The network matrices should not contain self-self connections. In other words, all values on the main diagonal of these matrices should be set to 0.
    # - In most cases, the network matrices should not contain negative weights. However, a substantial number of functions can process matrices with positive and negative weights. These functions typically end with sign.m (for signed networks).
    # - In general, randomization functions are designed for non-dense matrices; many randomization functions will be too slow and/or ineffective in dense matrices. However, some randomization functions are specifically designed for dense and weighted matrices.

    ## Generate a binary undirected network matrix using a density threshold.
    #
    # network density = (number of connections)/(number of possible connections)
    #
    #   number of connections = (sum of adjacency matrix 1s)/2 = (sum of node degrees)/2
    #   number of possible connections = n*(n-1)/2
    # 
    # The adjacency matrix is nxn, so we need to find the correlation threshold such that the adjacency matrix filled with a ratio of (n-1)/n*(density).
    # I.e. The adjacency matrix has n zeros which do not correspond to possible connections, so the proportion of nonzero elements in the adjacency matrix will be slightly less than the network density (a factor of (n-1)/2 less).

    percentile = ( 1 - network_density * float(num_nodes-1)/float(num_nodes) ) * 100.0
    correlation_thres = np.percentile(corr, percentile)

    weighted_adjacency= corr.copy()
    weighted_adjacency[corr < correlation_thres] = 0 # correlations below threshold are not considered connections; weights of connection are set to the correlation coefficient magnitude.
    
    binary_adjacency = np.zeros(corr.shape, dtype=np.float64)
    print "Warning: Saving adjacency matrix for binary network as 64-bit float for compatibilty with BCT. Saving as (int32/bool?) would save space."
    binary_adjacency[corr >= correlation_thres] = 1
    print "Thresholding correlation matrix to create an adjacency matrix for a network with density "+str(network_density)
    print "Actual network density:", float(binary_adjacency.sum())/ ( float(num_nodes)*float(num_nodes-1) )
    print "Corresponding correlation coefficient (either r or z depending on options): "+str(correlation_thres)
    if debug:
        print "Adjacency matrix for binary network:"
        print binary_adjacency[:10, :10]
        print

    # Calculate the node degree to compare with BNA
    node_degree = binary_adjacency.sum(axis=1)
    if debug:
        print "Node degree matrix:"
        print node_degree
        print node_degree.shape
        print

    # Calculate the node degree to compare with BNA
    node_degree_weighted = weighted_adjacency.sum(axis=1)
    if debug:
        print "Weighted node degree matrix:"
        print node_degree_weighted
        print node_degree_weighted.shape
        print
        
    #### Calculate measures from the correlation matrix.
    ## Calcuate the functional connectivity strength (FCS) as in Cao et al. (2017):
    # "Specifically, for each voxel, the FCS was calculated as the average of the correlations between this voxel and all other voxels in the brain."
    fcs_r = corr_r.mean(axis=1)
    if debug:
        print "FCS matrix:"
        print fcs_r
        print fcs_r.shape
        print
    fcs_z = corr_z.mean(axis=1)
    if debug:
        print "FCS matrix:"
        print fcs_z
        print fcs_z.shape
        print

    #### Map measures masked, flatten voxel arrays back to image space.
        
    ## Test mapping with the BOLD data. Here I want to check that I can get from the list of flattened, masked voxels back to the image space. I test on the bold image.
    # Reconstruct the bold series from the masked image.
    # Old way (works):
    #new_bold_flat = np.zeros(bold_flat.shape)
    #new_bold_flat[nonzero_bold_and_in_mask] = bold_mask # add the BOLD voxels that were included in the mask
    #new_bold = new_bold_flat.reshape(bold.shape)
    # New way:
    #new_bold = flatMaskedToImageSpace(bold_mask, nonzero_bold_and_in_mask, bold.shape[:3], measure_dim=1) # Here, the "measure" is a intensity time series at a single voxel, which is a 1D vector.
    #new_bold_nii = nib.Nifti1Image(new_bold, bold_nii.affine, bold_nii.header) # create NIFTI image with the new data, old header, and old affine.
    
    #new_bold_path = os.path.join(os.path.dirname(bold_path), "modified-"+os.path.basename(bold_path))
    #print "Saving test BOLD image to:", new_bold_path
    #nib.save(new_bold_nii, new_bold_path)
    
    # Node degree
    # save flat, masked numpy array for a  test.
#    node_degree_npy_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_bdeg.npy")
#    np.save(node_degree_npy_path, node_degree)    
#    node_degree_image_space = flatMaskedToImageSpace(node_degree, nonzero_bold_and_in_mask, bold.shape[:3])
#    node_degree_nii = nib.Nifti1Image(node_degree_image_space, bold_nii.affine, bold_nii.header)
#    node_degree_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_bdeg.nii.gz")
#    print "Saving node degree to:", node_degree_path
#    nib.save(node_degree_nii, node_degree_path)

    # FCS
    fcs_r_image_space = flatMaskedToImageSpace(fcs_r, nonzero_bold_and_in_mask, bold.shape[:3])
    fcs_r_nii = nib.Nifti1Image(fcs_r_image_space, bold_nii.affine, bold_nii.header)
    fcs_r_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_fcs_r.nii.gz")
    print "Saving FCS to:", fcs_r_path
    nib.save(fcs_r_nii, fcs_r_path)

    fcs_z_image_space = flatMaskedToImageSpace(fcs_z, nonzero_bold_and_in_mask, bold.shape[:3])
    fcs_z_nii = nib.Nifti1Image(fcs_z_image_space, bold_nii.affine, bold_nii.header)
    fcs_z_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_fcs_z.nii.gz")
    print "Saving FCS to:", fcs_z_path
    nib.save(fcs_z_nii, fcs_z_path)
    # If only FCS NIFTI file is desired, return from function now.
    if fcs_only: 
        return
    
    # Weighted node degree
    node_degree_weighted_image_space = flatMaskedToImageSpace(node_degree_weighted, nonzero_bold_and_in_mask, bold.shape[:3])
    node_degree_weighted_nii = nib.Nifti1Image(node_degree_weighted_image_space, bold_nii.affine, bold_nii.header)
    node_degree_weighted_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_wdeg.nii.gz")
    print "Saving weighted node degree to:", node_degree_weighted_path
    nib.save(node_degree_weighted_nii, node_degree_weighted_path)    

    # "Voxel" map, which is required to map from the voxels in the map back to the voxels in the image. Depends on the BOLD image, mask, and mask threshold.
    voxel_map_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_voxelmap.npy")
    print "Saving voxel map to:", voxel_map_path
    np.save(voxel_map_path, nonzero_bold_and_in_mask)
    
    ## Correlation matrix.
    corr_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_corr.npy")
    print "Saving correlation matrix to:", corr_path
    np.save(corr_path, corr)

    # Weighted adjacency matrix
    weighted_adjacency_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_wadj.npy")
    print "Saving weighted adjacency matrix to:", weighted_adjacency_path
    np.save(weighted_adjacency_path, weighted_adjacency)
    
    # Binary adjacency matrix
    binary_adjacency_path = os.path.join(out_dir, os.path.basename(bold_path).rstrip(".gz").rstrip(".nii")+"_badj.npy")
    print "Saving binary adjacency matrix to:", binary_adjacency_path
    np.save(binary_adjacency_path, binary_adjacency)

def mapMeasureToVoxels(bold_path, voxel_map_path, measure_path, out_dir):
    bold_nii = nib.load(bold_path)
    bold = np.asarray(bold_nii.get_fdata()).astype(np.float32)

    voxel_map = np.load(voxel_map_path)
    measure = np.load(measure_path)

    if (out_dir is None):
        out_dir = os.path.dirname(bold_path) # save it wherever the BOLD image is
    
    measure_image_space = flatMaskedToImageSpace(measure, voxel_map, bold.shape[:3])
    measure_nii = nib.Nifti1Image(measure_image_space, bold_nii.affine, bold_nii.header)
    measure_path = os.path.join(out_dir, os.path.basename(measure_path).rstrip(".npy")+".nii.gz")
    if os.path.exists(measure_path):
        message = "File exists; not overwriting: "+str(measure_path)
        #raise Exception(message)
        warnings.warn(message)
    else:
        print "Saving measure to:", measure_path
        nib.save(measure_nii, measure_path)    
    
if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate a voxel-wise functional connectivity matrix from fMRI data."""
    epilog = '' # """Text to follow argument explantion """
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("bold_path", type=str, help="path to BOLD time series in NIFTI format.")
    parser.add_argument("mask_path", type=str, help="path to mask in NIFTI format. If none provided, all nonzero voxels in the BOLD image will be used as nodes.")
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_dir", type=str, help="directory to save results to. Default: same as BOLD directory.")

    parser.add_argument("-t", "--mask_thres", type=float, help="include only voxels for which the mask has value greater than MASK_THRES. Default: 0", default=0.0)
    parser.add_argument("--corr_dist_min", type=float, help="set the correlation between pairs of voxels to zero if their centroids are within CORR_DIST_MIN of each other; enter in mm. Default=10", default=10.0)
    parser.add_argument("-d","--network_density", type=float, help="density of network for binary network. Density = (number of connections)/(number of potential connections). Default: 0.05", default=0.05)
    parser.add_argument("-z", "--fisher_r_to_z", action="store_true", help="apply Fisher's r-to-z transformation.")
    parser.add_argument("-a", "--no_absolute", action="store_true", help="do not take the absolute values of the correlations.")
    parser.add_argument("-f", "--fcs_only", action="store_true", help="Save only the functional connectivity strength NIFTI file.")
    
    parser.add_argument("--debug", action="store_true", help="debug mode - print lots of stuff.")
    parser.add_argument("--voxel_map_only", action="store_true", help="special option for debug")
    
    parser.add_argument("-v", "--measure_to_voxels", type=str, help="map a measure back to image space. Must include paths to the voxel map which was used when the connectivity matrix was constructed, and the measure (e.g. node degree), stored as NumPy arrays. The first axis of the measure array must have size equal to the number of nodes.", nargs=2, metavar=("VOXEL_MAP", "MEASURE"))    
    
    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    if args.measure_to_voxels:
        # Map measure (e.g. node degree) back to image space).
        mapMeasureToVoxels(bold_path=args.bold_path, voxel_map_path=args.measure_to_voxels[0], measure_path=args.measure_to_voxels[1], out_dir=args.out_dir)
    else:
        # Build connectivity matrix.
        buildConnectivityMatrix(**vars(args))

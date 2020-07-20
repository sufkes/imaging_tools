#!/usr/bin/env python

import os, sys
import argparse
import nibabel
import numpy as np

def demarcateCluster(volume, coord, include_func=None): # this function should only be entered once per cluster
    if (include_func == None):
        include_func = lambda x: x > 0.0
        
    dim_x = np.shape(volume)[0]
    dim_y = np.shape(volume)[1]
    dim_z = np.shape(volume)[2]

    coord_list = [coord] # list containing triplets of cluster coordinates
    seed_list = [coord]
    checked_list = [] # list of voxels which have been checked (>1) and used as a seed.
    while (len(seed_list) > 0):
        seed_coord = seed_list[0]
        neighbors = [] # list of 26 neighboring voxel coordinates.
        seed_ii = seed_coord[0]
        seed_jj = seed_coord[1]
        seed_kk = seed_coord[2]
        for d_ii in range(-1,2):
            for d_jj in range(-1,2):
                for d_kk in range(-1,2):
                    if (0 <= seed_ii+d_ii < dim_x) and (0 <= seed_jj+d_jj < dim_y) and (0 <= seed_kk+d_kk < dim_z): # if neighbor is inside volume
                        neighbors.append((seed_ii+d_ii, seed_jj+d_jj, seed_kk+d_kk))
        neighbors.remove(seed_coord) # the seed is not its own neighbor
        
        for neighbor in neighbors:
            voxel = volume[neighbor]
            if include_func(voxel): # if voxel is in a cluster
                if (not neighbor in checked_list): # if voxel has not been checked
                    if (not neighbor in coord_list):
                        coord_list.append(neighbor)
                    if (not neighbor in seed_list):
                        seed_list.append(neighbor)
        
        # Remove current voxel from seed list, and add it to the checked list.
        seed_list.remove(seed_coord)
        checked_list.append(seed_coord)        
    return coord_list

def findClusters(volume, include_func=None):
    clusters = [] # list containing clusters sets. Each cluster is a list of coordinate triplets.
    clusters_found = False
    for ii in range(np.shape(volume)[0]):
        for jj in range(np.shape(volume)[1]):
            for kk in range(np.shape(volume)[2]):
                coord = (ii,jj,kk)
                voxel = volume[ii,jj,kk]
                if include_func(voxel): # if voxel is in a cluster
                    if (clusters_found): # if a cluster has already been identified
                        # Look for the current voxel in known clusters.
                        cluster_known = False
                        for cluster in clusters:
                            if coord in cluster: # if voxel is in a known cluster (i.e. if the cluster has already been identified)
                                cluster_known = True
                                break
                    else: # if no clusters have been identified
                        cluster_known = False 
                        clusters_found = True

                    # Demarcate cluster if it is currently unknown.
                    if (not cluster_known): # if voxel is in an unknown cluster
                        cluster_coord_list = demarcateCluster(volume, coord, include_func) # coordinates of all voxels in cluster.
                        clusters.append(cluster_coord_list)
    return clusters

def getTotalClusterVol(clusters, vol_vox):
    """Return total volume of all clusters in clusters"""
    vol = 0
    for cluster in clusters:
        vol += len(cluster)
    vol = vol*vol_vox
    return vol

def getTotalBigClusterVol(clusters, vol_vox, vol_thres=1600):
    """Return total volume of all clusters in clusters which have a volume greater than vol_thres (in mm^3)"""
    vol = 0
    for cluster in clusters:
        if len(cluster)*vol_vox >= vol_thres:
            vol += len(cluster)
    vol = vol*vol_vox
    return vol
    
def clusterStats(path_unmasked, path_masked):
    # Open NIFTI image.
    img_unmasked = nibabel.load(path_unmasked)
    img_masked = nibabel.load(path_masked)
    
    zooms = img_unmasked.header.get_zooms()
    vol_vox = zooms[0]*zooms[1]*zooms[2] # volume of a voxel in mm; assume the same with and without mask.

    volume_unmasked = img_unmasked.get_fdata() # pixel data in a numpy array
    volume_masked = img_masked.get_fdata()

    # Find clusters of significant "activation" (z>0) and "deactivation" (z<0).
    def pos(x):
        return x>0.0
    def neg(x):
        return x<0.0
    clusters_activation_unmasked = findClusters(volume_unmasked, include_func=pos) # list of lists of (x,y,z) coordinate tuples
    clusters_activation_masked = findClusters(volume_masked, include_func=pos) # list of lists of (x,y,z) coordinate tuples
    clusters_deactivation_unmasked = findClusters(volume_unmasked, include_func=neg)
    clusters_deactivation_masked = findClusters(volume_masked, include_func=neg)

    # Get volumes
    vol_activation_unmasked = getTotalClusterVol(clusters_activation_unmasked, vol_vox)
    vol_activation_masked = getTotalClusterVol(clusters_activation_masked, vol_vox)
    vol_deactivation_unmasked = getTotalClusterVol(clusters_deactivation_unmasked, vol_vox)
    vol_deactivation_masked = getTotalClusterVol(clusters_deactivation_masked, vol_vox)

    vol_significant_unmasked = vol_activation_unmasked + vol_deactivation_unmasked # total volume of all significant voxels.
    vol_significant_masked = vol_activation_masked + vol_deactivation_masked # total volume of all significant voxels.
    
    vol_big_cluster_activation_masked = getTotalBigClusterVol(clusters_activation_masked, vol_vox)
    vol_big_cluster_deactivation_masked = getTotalBigClusterVol(clusters_deactivation_masked, vol_vox)
    vol_big_cluster_significant_masked = vol_big_cluster_activation_masked + vol_big_cluster_deactivation_masked # volume of significant voxels lying in "big" clusters within regions of interest (not in WM or CSF).

    ratio_non_wm_csf = vol_significant_masked/vol_significant_unmasked
    ratio_relevant = vol_big_cluster_significant_masked/vol_significant_unmasked # signal if >10%, noise if <10%
    
    print "Volume of significant voxels                      :", vol_significant_unmasked
    print "Volume of significant voxels outside of WM and CSF:", vol_significant_masked
    print "Volume of 'relevant' voxels                       :", vol_big_cluster_significant_masked    
    print "Ratio V(significant non-WM/CSF)/V(significant)    :", ratio_non_wm_csf
    print "Ratio V(relevant)/V(significant)                  :", ratio_relevant
    return ratio_relevant


if (__name__ == '__main__'):
    # Create argument parser
    description = """Input a NIFTI image of a thresholded z statistic output from FSL MELODIC, and the same image after masking to exclude voxels in CSF and WM. Report information about significant clusters as described in R. Kelly et al. 'Visual inspection of independent components: Defining a procedure for artifact removal from fMRI data' Journal of Neuroscience Methods 189 (2010) 233-245 """
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("path_unmasked", help="path to NIFTI file of z statistic output by MELODIC")
    parser.add_argument("path_masked", help="path to NIFTI file of z statistic output by MELODIC with CSF and WM voxels excluded")
    
    # Define optional arguments.

    # Parse arguments.
    args = parser.parse_args()

    # Do stuff.
    clusterStats(args.path_unmasked, args.path_masked)

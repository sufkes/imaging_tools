#!/usr/bin/env python3

import os
import sys
import argparse
from collections import OrderedDict
from pprint import pprint

import numpy as np
import pandas as pd
import pyminc.volumes.factory as pyminc
import nibabel as nib
from skimage import measure

#class Cluster(object):
#    """Class storing data for a single cluster of labelled voxels"""
#    def __init__(self, clusters_labelled, cluster_val):
#        self.setClusterData(clusters_labelled, cluster_val)
#
#    def setClusterData(self, clusters_labelled, cluster_val):
#        cluster_data = clusters_labelled.copy()
#        cluster_data[cluster_data != cluster_val] = 0 # Remove data corresponding to other clusters.
#        cluster_data[cluster_data == cluster_val] = 1 # Set values for the current cluster to 1.
#        self.cluster_data = cluster_data

class LabelVolume(object):
    """Class storing data corresponding to a single unique label value"""
    def __init__(self, data, label_val):
        self.setLabelData(data, label_val)
        self.setClusters()
        
    def setLabelData(self, data, label_val):
        """Define the array corresponding to a single unique label"""
        label_data = data.copy()
        label_data[label_data != label_val] = 0 # Remove data corresponding to a different label.
        label_data[label_data == label_val] = 1 # Set values for the current label 1 for this LabelVolume.
        self.label_data = label_data
        
    def setClusters(self):
        # Convert the single-valued label data to separate values for each cluster.
        clusters_labelled = measure.label(self.label_data, connectivity=3)
        clusters = []
        for cluster_val in np.unique(clusters_labelled):
            if cluster_val == 0:
                continue
            #cluster = Cluster(clusters_labelled, cluster_val)
            cluster = clusters_labelled.copy()
            cluster[cluster != cluster_val] = 0 # Remove data corresponding to other clusters.
            cluster[cluster == cluster_val] = 1 # Set values for the current cluster to 1.
            
            clusters.append(cluster)
        self.clusters = clusters
        
class LabelImage(object):
    """Class storing all relevant data from an input NIFTI or MINC file"""
    def __init__(self, file_path, combine_labels=False):
        self.file_path = file_path
        self.combine_labels = combine_labels
        self.setFileType()
        self.readFile() # get the array data and the size of each pixel along each direction
        self.setPixVol()
        self.setLabelVolumes()
        self.setStats()
        
    def setFileType(self):
        """Determine the type of the input file"""
        if self.file_path.endswith('.mnc'):
            self.file_type = 'minc'
        elif self.file_path.endswith('.nii') or self.file_path.endswith('.nii.gz'):
            self.file_type = 'nifti'
        else:
            raise Exception('Unknown file extension')

    def readFile(self):
        """Read the pixel data and dimensions from file."""
        if self.file_type == 'minc':
            self.minc = pyminc.volumeFromFile(self.file_path, labels=True)
            self.nifti = None
            self.pixdim = np.abs(np.array(self.minc.separations), dtype=float)
            self.data = np.array(self.minc.data, dtype=int)

            # MINC files apparently sometimes store the x, y, z dimensions in strange orders. Try to fix this such that the data and pixdim arrays are ordered (x, y, z)
            desired_dimnames = ['xspace', 'yspace', 'zspace']
            transpose_axes = [0, 0, 0] # generate 'axes' argument for numpy.transpose(), to rearrange the data.
            if self.minc.dimnames != desired_dimnames:
                # Reorder the pixdim and data arrays.
                new_pixdim = np.zeros(shape=(3,), dtype=float)
                for desired_index, desired_dimname in enumerate(desired_dimnames):
                    current_index = self.minc.dimnames.index(desired_dimname)
                    new_pixdim[desired_index] = self.pixdim[current_index]
                    transpose_axes[desired_index] = current_index
                self.pixdim = new_pixdim
                self.data = np.transpose(self.data, axes=transpose_axes)

        elif self.file_type == 'nifti':
            self.minc = None
            self.nifti = nib.load(self.file_path)
            self.pixdim = np.array(self.nifti.header['pixdim'][1:4], dtype=float)
            self.data = np.array(self.nifti.get_fdata(), dtype=int)

        # Combine unique labels if requested.
        if self.combine_labels:
            self.data[self.data != 0] = 1 # Set all nonzero values to 1.

    def setPixVol(self):
        """Get the volume of a pixel in mm^3."""
        self.pixvol = np.product(self.pixdim)

    def setLabelVolumes(self):
        self.label_volumes = OrderedDict()
        for label_val in np.unique(self.data):
            if label_val == 0:
                continue
            self.label_volumes[label_val] = LabelVolume(self.data, label_val)

    def calculateExtent(self, cluster):
        # Calculate the greatest distance between a pair of voxels in the cluster. Use brute force since I can't find a slick way that is easy to understand.
        roi_coords = np.argwhere(cluster > 0) # Nx3 array with all coordinates lying in cluster.
        num_voxels = roi_coords.shape[0]
        extent = 0
        for v1 in range(num_voxels):
            for v2 in range(v1+1, num_voxels):
                c1 = roi_coords[v1] # coordinates of first voxel
                c2 = roi_coords[v2] # coordinates of second voxel
                d = np.abs((c2 - c1) * self.pixdim) + self.pixdim # Vector whose length is equal to the longest distance between the two voxels, accounting for the size of the voxels themselves.
                length = np.sqrt(d.dot(d))
                #print('c1:',c1)
                #print('c2:',c2)
                #print('pixdim:', self.pixdim)
                #print('d:',d)
                #print('length:',length)
                if length > extent:
                    extent = length
        return extent
            
    def setStats(self):
        """Compute stats for each label in the image."""
        image_stats = OrderedDict() # key is label_val, value is stats dictionary

        for label_val, label_volume in self.label_volumes.items():
            label_stats = {}
            num_clusters = len(label_volume.clusters)
            
            volumes_px = np.zeros(num_clusters, dtype=float)
            extents_mm = np.zeros(num_clusters, dtype=float)
            
            for cluster_num, cluster in enumerate(label_volume.clusters):
                volume_px = np.count_nonzero(cluster)
                extent_mm = self.calculateExtent(cluster)
                #region_props = measure.regionprops(cluster) # this nicely calculates a bunch of different stats for the cluster, non of which I really need.
                volumes_px[cluster_num] = volume_px
                extents_mm[cluster_num] = extent_mm

            label_stats['num_clusters'] = len(label_volume.clusters)

            label_stats['volume_tot_px'] = volumes_px.sum()
            label_stats['volume_mean_px'] = volumes_px.mean()
            label_stats['volume_max_px'] = volumes_px.max()
            label_stats['volume_tot_mm3'] = label_stats['volume_tot_px']*self.pixvol
            label_stats['volume_mean_mm3'] = label_stats['volume_mean_px']*self.pixvol
            label_stats['volume_max_mm3'] = label_stats['volume_max_px']*self.pixvol

            label_stats['extent_mean_mm'] = extents_mm.mean()
            label_stats['extent_max_mm'] = extents_mm.max()
            
            image_stats[label_val] = label_stats
            
        self.image_stats = image_stats

def reportStats(file_path, img, datasheet_path):
    ## Add stats for current image to dataframe.
    if (not datasheet_path is None) and os.path.isfile(datasheet_path):
        df = pd.read_csv(datasheet_path)
    else:
        df = pd.DataFrame()        
    
    file_name = os.path.basename(file_path)
    for label_val, label_stats in img.image_stats.items():
        ## Add row for current (file_name, label_val):
        new_index = len(df)
        df.loc[new_index, 'file_name'] = file_name
        df.loc[new_index, 'label'] = label_val
        for stat_name, stat_val in label_stats.items():
            df.loc[new_index, stat_name] = stat_val

    if not datasheet_path is None:
        df.to_csv(datasheet_path, index=False)
    print(df.to_string())
    return df

def getRoiStats(file_path, datasheet_path, combine_labels):
    """Input a 3D NIFTI or MINC file containing regions of interest (ROIs) labelled with integers. For each integer label, report statistics on size and position. Pairs of labelled voxels lying adjacent to each other along a square or cubic diagonal are assumed to belong to the same region."""

    ## Open image and compute statistics on the labelled ROIs.
    img = LabelImage(file_path, combine_labels=combine_labels)

    ## Report the statistics.
    reportStats(file_path, img, datasheet_path)
    
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = getRoiStats.__doc__
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("file_path", type=str, help="path to input NIFTI or MINC file")
    
    # Define optional arguments.
    parser.add_argument("-d", "--datasheet_path", type=str, help='Path to CSV containing statistics for input images. If file exists, data will be appended to it.')
    parser.add_argument("-c", "--combine_labels", help="Combine all nonzero labels. By default, statistics are computed for each nonzero label separately. If this flag is set, all nonzero labels will be treated as one. This may be useful, for example, if lesions in the left and right hemispheres are labelled differently, but you want statistics across all lesions.", action='store_true')

    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    getRoiStats(**vars(args))

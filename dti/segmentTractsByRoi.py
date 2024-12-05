#!/usr/bin/env python3

import os
import sys
import argparse

from pprint import pprint

from dipy.io.streamline import load_trk, save_trk
from dipy.tracking._utils import _mapping_to_voxel, _to_voxel_coordinates
from dipy.tracking.utils import target, density_map
from dipy.io.stateful_tractogram import StatefulTractogram, Space

import nibabel as nib

import numpy as np

class Metric(object):
    def __init__(self, in_metric_path, out_metric_matrix_path, num_rois):
        metric_nii = nib.load(in_metric_path)
        self.in_metric_path = in_metric_path
        self.out_metric_matrix_path = out_metric_matrix_path
        self.nii = metric_nii
        self.data = metric_nii.get_fdata()
        self.affine = metric_nii.affine
        self.metric_matrix = np.zeros(shape=(num_rois, num_rois)) # connectivity matrix based on this metric

class TractManager(object):
    def __init__(self, trk_path, label_path, skip_segmentation, include_terminal_rois, exclude_terminal_rois, exclude_terminal_roi_interiors, out_segmented_trk_path,  out_nos_matrix_path, in_metric_paths, out_metric_matrix_paths, save_roi_pairs, out_roi_pair_trk_paths, out_roi_pair_density_paths, verbose):
        self.verbose = verbose
        if self.verbose:
            print('Parsing inputs and loading files.')
            
        ## Parse arguments, checking for errors.
        self.trk_path = trk_path
        self.label_path = label_path

        self.skip_segmentation = skip_segmentation
        
        # ROI pair arguments (mutually exclusive).
        if (not include_terminal_rois is None) and (not exclude_terminal_rois is None):
            msg = 'Use "--include_terminal_rois", "--exclude_terminal_rois", or neither, but not both'
            raise Exception(msg)
        self.include_terminal_rois = include_terminal_rois
        self.exclude_terminal_rois = exclude_terminal_rois

        self.exclude_terminal_roi_interiors = exclude_terminal_roi_interiors
        self.out_segmented_trk_path = out_segmented_trk_path
        self.out_nos_matrix_path = out_nos_matrix_path

        # Metric arguments (number of outputs must equal number of inputs if outputs requested).
        if (not out_metric_matrix_paths is None):
            if (in_metric_paths is None) or (len(in_metric_paths) != len(out_metric_matrix_paths)):
                msg = 'Number of --out_metric_matrix_paths arguments must match number of --in_metric_paths arguments.'
                raise Exception(msg)
        self.in_metric_paths = in_metric_paths
        self.out_metric_matrix_paths = out_metric_matrix_paths
        
        # Specified ROI pairs (number of outputs must match number of specified ROI pairs).
        if (not save_roi_pairs is None):
            roi_pairs = []
            try:
                for pair_string in save_roi_pairs:
                    roi_pair = tuple(int(label.strip()) for label in pair_string.split(','))
                    assert len(roi_pair) == 2
                    roi_pairs.append(roi_pair)
            except AssertionError:
                msg = 'Each pair of ROIs in "--save_roi_pairs" must be specified by a comma-separated pair of integers, and each pair of ROIs must be separated by a space. For example, to specify the ROI pairs 1<->2 and 1<->3, use "--save_roi_pairs 1,2 1,3'
                raise Exception(msg)
            self.save_roi_pairs = roi_pairs
        
        if (not out_roi_pair_trk_paths is None):
            if (save_roi_pairs is None) or (len(save_roi_pairs) != len(out_roi_pair_trk_paths)):
                msg = 'Number of "--out_roi_pair_trk_paths" arguments must match number of "--save_roi_pairs" arguments.'
                raise Exception(msg)
        if (not out_roi_pair_density_paths is None):
            if (save_roi_pairs is None) or (len(save_roi_pairs) != len(out_roi_pair_density_paths)):
                msg = 'Number of "--out_roi_pair_density_paths" arguments must match number of "--save_roi_pairs" arguments.'
                raise Exception(msg)

        self.out_roi_pair_trk_paths = out_roi_pair_trk_paths
        self.out_roi_pair_density_paths = out_roi_pair_density_paths

        ## Load mandatory inputs.
        self.loadStreamlines()
        self.loadLabels()

    def loadStreamlines(self):
        tractogram = load_trk(self.trk_path, reference='same')
        self.tractogram_original = tractogram
        self.streamlines_original = tractogram.streamlines
        
    def loadLabels(self):
        label_nii = nib.load(self.label_path)
        self.label_nii = label_nii
        label_data = label_nii.get_fdata()
        try:
            self.label_data = label_data.astype(int)
        except:
            msg = 'Failed to convert label image data type integer.'
            raise Exception(msg)
        self.label_affine = label_nii.affine

    def setMappingToVoxelCoordinates(self):
        self.mapping_to_voxel_linear_transformation, self.mapping_to_voxel_offset = _mapping_to_voxel(self.label_affine)
        
    def computeStreamlineLabels(self, streamline):
        """Input streamline of shape (3, N), return array of length N indicating the ROI label at each point in the streamline"""

        # Logic from https://github.com/dipy/dipy/blob/master/dipy/tracking/utils.py
        streamline_voxel_coordinates = _to_voxel_coordinates(streamline, self.mapping_to_voxel_linear_transformation, self.mapping_to_voxel_offset)
        i, j, k = streamline_voxel_coordinates.T
        streamline_labels = self.label_data[i, j, k]
        return streamline_labels

    def setAllStreamlineLabels(self):
        if self.verbose:
            print('Getting labels for streamlines.')
        self.setMappingToVoxelCoordinates()
        
        #all_streamline_labels = nib.streamlines.array_sequence.ArraySequence() # way slower than a Python list for some reason.
        all_streamline_labels = []
        for streamline in self.streamlines_original:
            streamline_labels = self.computeStreamlineLabels(streamline)
            all_streamline_labels.append(streamline_labels)
        self.all_streamline_labels = all_streamline_labels

    def computeStreamlineSegments(self, streamline, streamline_labels):
        """Split streamline into segments that terminate at each pair of valid ROIs."""
        boundary_list = []
        # Add background label to ends of the streamline so that they are automatically treated as boundaries.
        padded_streamline_labels = np.zeros(shape=streamline_labels.shape[0]+2, dtype=streamline_labels.dtype)
        padded_streamline_labels[1:-1] = streamline_labels
        padded_streamline = np.zeros(shape=(streamline.shape[0]+2, 3), dtype=streamline.dtype)
        padded_streamline[1:-1, :] = streamline
        streamline_labels = padded_streamline_labels
        streamline = padded_streamline

        # Find boundaries between labels along the streamline.
        for index in range(len(streamline_labels)-1):
            if streamline_labels[index] != streamline_labels[index+1]:
                boundary_list.append( (index, streamline_labels[index], streamline_labels[index+1]) )

        #print(streamline)
        #print('streamline_labels:')
        #print(streamline_labels)
        #print('')
        #print('unique labels:')
        #print(np.unique(np.array(streamline_labels)))
        #print('')
        #print([x[1:] for x in boundary_list])
        #print('')
                
        streamline_segments = []
        streamline_segments_labels = []

        if not self.exclude_terminal_roi_interiors:
            for boundary_list_index, boundary_tuple in enumerate(boundary_list):
                origin = boundary_tuple[2]
                # If entering the background, continue, because we want to start streamline segments where regions of interest begin.
                if boundary_tuple[2] == 0:
                    continue

                # Labels of other ROIs to which the current departure ROI connects.
                destinations = []

                # Loop over boundaries further down the streamline.
                for next_boundary_tuple in boundary_list[boundary_list_index+2:]: # by definition, the next boundary is the edge of the starting ROI, so we never want to consider it; therefore start at second-next boundary.
                    destination = next_boundary_tuple[1]
                    # If destination is the same region from which we started (i.e. if looping back to the same ROI), stop search.
                    if origin == destination:
                        break

                    # If destination is the background, do not record.
                    if destination == 0:
                        continue

                    # If destination was already recorded, continue
                    if destination in destinations:
                        continue

                    # Record segment leading to new destination.
                    start_index = boundary_tuple[0]+1
                    stop_index = next_boundary_tuple[0]+1
                    segment_labels = streamline_labels[start_index: stop_index]
                    segment = streamline[start_index: stop_index, :]
                    streamline_segments_labels.append(segment_labels)
                    streamline_segments.append(segment)
            
                    destinations.append(destination)
        else:        
            for boundary_list_index, boundary_tuple in enumerate(boundary_list):
                origin = boundary_tuple[1]
                # If entering a region from the background, continue, because we want to start streamline segments when we leave a region of interest.
                if origin == 0:
                    continue

                # Labels of other ROIs to which the current departure ROI connects.
                destinations = []

                # If this is a boundary between two ROIs, record two-point streamline.
                if boundary_tuple[2] != 0:
                    destination = boundary_tuple[2]
                
                    start_index = boundary_tuple[0]
                    stop_index = boundary_tuple[0] + 2
                    segment_labels = streamline_labels[start_index: stop_index]
                    segment = streamline[start_index: stop_index, :]
                    streamline_segments_labels.append(segment_labels)
                    streamline_segments.append(segment)

                    destinations.append(destination)

                # Loop over boundaries further down the streamline.
                for next_boundary_tuple in boundary_list[boundary_list_index+1:]:
                    destination = next_boundary_tuple[2]
                    # If entering the same region from which we started (i.e. if looping back to the same ROI), stop search.
                    if origin == destination:
                        break

                    # If destination is the background, do not record.
                    if destination == 0:
                        continue

                    # If destination was already recorded, do not record.
                    if destination in destinations:
                        continue

                    # Record segment leading to new destination.
                    start_index = boundary_tuple[0]
                    stop_index = next_boundary_tuple[0]+2 # include the entire next boundary in the streamline (+1), and Python stops at following index (+1)
                    segment_labels = streamline_labels[start_index: stop_index]
                    segment = streamline[start_index: stop_index, :]
                    streamline_segments_labels.append(segment_labels)
                    streamline_segments.append(segment)
            
                    destinations.append(destination)

        streamline_segments = nib.streamlines.array_sequence.ArraySequence(streamline_segments) # Convert back to ArraySequence for consistency with dipy/nibabel streamlines.
        #print(streamline_segments_labels)
        #print('-'*50)
        return streamline_segments, streamline_segments_labels

    def setValidTerminalRois(self):
        if not self.include_terminal_rois is None:
            valid_terminal_rois = np.unique(self.include_terminal_rois)
        elif not self.exclude_terminal_rois is None:
            all_labels = [n for n in np.unique(self.label_data) if n != 0]
            valid_terminal_rois = np.array([n for n in all_labels if not n in self.exclude_terminal_rois])
        else:
            valid_terminal_rois = np.array([n for n in np.unique(self.label_data) if n != 0])
        valid_terminal_rois.sort()
        if 0 in valid_terminal_rois:
            msg = 'Terminal ROIs may not include 0.'
            raise Exception(msg)
        self.valid_terminal_rois = valid_terminal_rois
        self.num_valid_terminal_rois = len(valid_terminal_rois)

    def streamlineTerminatesAtRoiPair(self, streamline_labels, roi_pair):
        streamline_terminates_at_roi_pair = (streamline_labels[0] in roi_pair) and (streamline_labels[-1] in roi_pair)
        return streamline_terminates_at_roi_pair

    def findStreamlinesThatTerminateAtRoiPair(self, streamlines, streamlines_labels, roi_pair, return_labels=True):
        streamlines_that_terminate_at_roi_pair = [] # faster to append to Python list than nib.streamlines.array_sequence.ArraySequence
        if return_labels:
            streamlines_labels_that_terminate_at_roi_pair = []
        if return_labels: # do separate loops depending on outcome because maybe num_rois^2*num_streamlines extra if statements will take a noticeable amount of time.
            for streamline, streamline_labels in zip(streamlines, streamlines_labels):
                if self.streamlineTerminatesAtRoiPair(streamline_labels, roi_pair):
                    streamlines_that_terminate_at_roi_pair.append(streamline)
                    streamlines_labels_that_terminate_at_roi_pair.append(streamline_labels)
        else:
            for streamline, streamline_labels in zip(streamlines, streamlines_labels):
                if self.streamlineTerminatesAtRoiPair(streamline_labels, roi_pair):
                    streamlines_that_terminate_at_roi_pair.append(streamline)
        streamlines_that_terminate_at_roi_pair = nib.streamlines.array_sequence.ArraySequence(streamlines_that_terminate_at_roi_pair) # convert to nib.streamlines.array_sequence.ArraySequence for consistency with nibabel/dipy
        if return_labels:
            return_val = (streamlines_that_terminate_at_roi_pair, streamlines_labels_that_terminate_at_roi_pair)
        else:
            return_val = streamlines_that_terminate_at_roi_pair
        return return_val

    def streamlineTerminatesAtRois(self, streamline_labels, rois):
        streamline_terminates_at_rois = (streamline_labels[0] in rois) and (streamline_labels[-1] in rois)
        return streamline_terminates_at_rois
    
    def findStreamlinesThatTerminateAtRois(self, streamlines, streamlines_labels, rois):
        streamlines_that_terminate_at_rois = []
        streamlines_labels_that_terminate_at_rois = [] # faster to append to Python list than nib.streamlines.array_sequence.ArraySequence
        for streamline, streamline_labels in zip(streamlines, streamlines_labels):
            if self.streamlineTerminatesAtRois(streamline_labels, rois):
                streamlines_that_terminate_at_rois.append(streamline)
                streamlines_labels_that_terminate_at_rois.append(streamline_labels)
        streamlines_that_terminate_at_rois = nib.streamlines.array_sequence.ArraySequence(streamlines_that_terminate_at_rois) # convert to nib.streamlines.array_sequence.ArraySequence for consistency with nibabel/dipy
        return streamlines_that_terminate_at_rois, streamlines_labels_that_terminate_at_rois
    
    def removeInvalidStreamlineSegments(self):
        self.setValidTerminalRois()
        
        streamlines_segmented_valid, streamlines_segmented_labels_valid = self.findStreamlinesThatTerminateAtRois(self.streamlines_segmented, self.streamlines_segmented_labels, self.valid_terminal_rois)
        self.streamlines_segmented = streamlines_segmented_valid
        self.streamlines_segmented_labels = streamlines_segmented_labels_valid
    
    def segmentStreamlines(self):        
        if self.skip_segmentation:
            self.streamlines_segmented = self.streamlines_original
        else:
            if self.verbose:
                print('Segmenting streamlines.')
            streamlines_segmented = [] # Convert to nib.streamlines.array_sequence.ArraySequence later for consistency with nibabel, but appending to Python list is way faster.
            streamlines_segmented_labels = []

            for streamline, streamline_labels in zip(self.streamlines_original, self.all_streamline_labels):
                streamline_segments, streamline_segments_labels = self.computeStreamlineSegments(streamline, streamline_labels)
                streamlines_segmented.extend(streamline_segments)
                streamlines_segmented_labels.extend(streamline_segments_labels)

            streamlines_segmented = nib.streamlines.array_sequence.ArraySequence(streamlines_segmented)
                    
            self.streamlines_segmented = streamlines_segmented
            self.streamlines_segmented_labels = streamlines_segmented_labels

            # Remove segments terminating in invalid regions, as specified by user.
            self.removeInvalidStreamlineSegments()

    def saveTractogram(self, streamlines, out_path):
        tractogram = StatefulTractogram(streamlines, self.label_nii, Space.RASMM)
        save_trk(tractogram, out_path)
        
    def saveSegmentedTractogram(self):
        if (not self.out_segmented_trk_path is None):
            if self.verbose:
                print('Saving segmented streamlines to TRK file.')
            self.saveTractogram(self.streamlines_segmented, self.out_segmented_trk_path)
            
    def generateSpecifiedTracts(self):
        if (not self.out_roi_pair_trk_paths is None):
            if self.verbose:
                print('Generating TRK files for specified ROI pairs.')
            for roi_pair, out_roi_pair_trk_path in zip(self.save_roi_pairs, self.out_roi_pair_trk_paths):
                streamlines_segmented_filtered = self.findStreamlinesThatTerminateAtRoiPair(self.streamlines_segmented, self.streamlines_segmented_labels, roi_pair, return_labels=False)
                self.saveTractogram(streamlines_segmented_filtered, out_roi_pair_trk_path)

    def calculateStreamlineDensity(self, streamlines):
        streamline_density = density_map(streamlines, affine=self.label_affine, vol_dims=self.label_nii.shape).astype(float)
        return streamline_density
                
    def calculateStreamlineDensityConnectingRoiPair(self, streamlines, streamlines_labels, roi_pair):
        streamlines_that_terminate_at_roi_pair = self.findStreamlinesThatTerminateAtRoiPair(streamlines, streamlines_labels, roi_pair, return_labels=False)
        streamline_density_data = self.calculateStreamlineDensity(streamlines_that_terminate_at_roi_pair)
        return streamline_density_data
                
    def generateSpecifiedStreamlineDensityMaps(self):
        if (not self.out_roi_pair_density_paths is None):
            if self.verbose:
                print('Generating streamline density maps for specified ROI pairs.')
            for roi_pair, out_roi_pair_density_path in zip(self.save_roi_pairs, self.out_roi_pair_density_paths):
                streamline_density_data = self.calculateStreamlineDensityConnectingRoiPair(self.streamlines_segmented, self.streamlines_segmented_labels, roi_pair)
                streamline_density_nii = nib.Nifti1Image(streamline_density_data, self.label_affine)
                nib.save(streamline_density_nii, out_roi_pair_density_path)

    def setMetricsList(self):
        if (not self.out_metric_matrix_paths is None):
            metrics_list = []
            for in_metric_path, out_metric_matrix_path in zip(self.in_metric_paths, self.out_metric_matrix_paths):
                metric = Metric(in_metric_path, out_metric_matrix_path, self.num_valid_terminal_rois)
                metrics_list.append(metric)
            self.metrics_list = metrics_list
                
    def generateConnectivityMatrices(self):
        if (not self.out_metric_matrix_paths is None) or (not self.out_nos_matrix_path is None):
            if self.verbose:
                print('Generating connectivity matrices.')

            if (not self.out_nos_matrix_path is None):
                nos_matrix = np.zeros(shape=(self.num_valid_terminal_rois, self.num_valid_terminal_rois))
            if (not self.out_metric_matrix_paths is None):
                self.setMetricsList()
                
            # Loop over all valid ROI pairs. Since connectivity matrices will be symmetric, compute the upper-triangle part of the matrix, then add transpose later.
            num_connections = self.num_valid_terminal_rois*(self.num_valid_terminal_rois - 1)/2
            count = 0
            if self.verbose:
                print(f'  Computing strengths of {num_connections} connections...')
            for index_1, roi_1 in enumerate(self.valid_terminal_rois[:-2]):
                for index_2, roi_2 in enumerate(self.valid_terminal_rois[index_1 + 1:], index_1+1):
                    roi_pair = (roi_1, roi_2)
                    streamlines_connecting_roi_pair = self.findStreamlinesThatTerminateAtRoiPair(self.streamlines_segmented, self.streamlines_segmented_labels, roi_pair, return_labels=False)
                    num_streamlines_connecting_roi_pair = len(streamlines_connecting_roi_pair)

                    # Add element to number of streamlines matrix.
                    if (not self.out_nos_matrix_path is None):
                        nos_matrix[index_1, index_2] = num_streamlines_connecting_roi_pair

                    # Compute elements of metric-mean matrices.
                    if (not self.out_metric_matrix_paths is None):
                        if num_streamlines_connecting_roi_pair > 0:
                            streamline_density_connecting_roi_pair = self.calculateStreamlineDensity(streamlines_connecting_roi_pair)

                        for metric in self.metrics_list:
                            # Compute metric mean weighted by streamline density.
                            if num_streamlines_connecting_roi_pair > 0:                               
                                weighted_mean = np.multiply(metric.data, streamline_density_connecting_roi_pair).sum() / streamline_density_connecting_roi_pair.sum()
                            else:
                                weighted_mean = 0
                            metric.metric_matrix[index_1, index_2] = weighted_mean
                    count += 1
                    if self.verbose:
                        if count % 10 == 0:
                            print(f'\r    Finished {count}/{num_connections}', end='', flush=True)
            if self.verbose:
                print('')

            # Add matrix transposes in order to fill elements at [i,j] where i>j.
            if (not self.out_nos_matrix_path is None):
                nos_matrix += nos_matrix.T
            if (not self.out_metric_matrix_paths is None):
                for metric in self.metrics_list:
                    metric.metric_matrix += metric.metric_matrix

            # Save matrices.
            if (not self.out_nos_matrix_path is None):
                np.savetxt(self.out_nos_matrix_path, nos_matrix)
            if (not self.out_metric_matrix_paths is None):
                for metric in self.metrics_list:
                    np.savetxt(metric.out_metric_matrix_path, metric.metric_matrix)
    
    def debug(self):
        #print(self.label_nii)
        #print(self.label_data)
        #print(self.label_affine)
        #print(self.tractogram_original)
        #print(self.streamlines_original)
        #pprint(dir(self))
        #print(len(self.streamlines_original))
        #print(self.streamlines_original[70].shape)
        #print(self.all_streamline_labels[70].shape)
        #print(len(self.streamlines_segmented_labels))
        #print(len(self.streamlines_segmented))
        #print(len(self.streamlines_segmented_labels[50]))
        #print(len(self.streamlines_segmented[50]))
        #print(self.streamlines_segmented_labels)
        #pprint([(x[0],x[-1]) for x in self.streamlines_segmented_labels if (x[0] == 46) or (x[-1] == 46)])
        print('')
        print(f'rows and columns in connectivity matrices correspond to the following ROI labels: {self.valid_terminal_rois}')
        pass
            
def main(trk_path, label_path, skip_segmentation, include_terminal_rois, exclude_terminal_rois, exclude_terminal_roi_interiors, out_segmented_trk_path, out_nos_matrix_path, in_metric_paths, out_metric_matrix_paths, save_roi_pairs, out_roi_pair_trk_paths, out_roi_pair_density_paths, verbose):
    tract_manager = TractManager(trk_path, label_path, skip_segmentation, include_terminal_rois, exclude_terminal_rois, exclude_terminal_roi_interiors, out_segmented_trk_path, out_nos_matrix_path, in_metric_paths, out_metric_matrix_paths, save_roi_pairs, out_roi_pair_trk_paths, out_roi_pair_density_paths, verbose)

    # Get labels for all streamlines.
    tract_manager.setAllStreamlineLabels()

    # Segment the streamlines.
    tract_manager.segmentStreamlines()

    # Save segmented tractogram.
    #tract_manager.saveSegmentedTractogram()
    
    # Compute connectivity matrices based on number of streamlines, and mean values of metrics in segmeneted tracts weighted by streamline density.
    tract_manager.generateConnectivityMatrices()

    # Generate TRK files including only streamlines terminating at specified pairs of ROIs
    tract_manager.generateSpecifiedTracts()
    
    # Generate streamline density maps including only streamlines terminating at specified pairs of ROIs
    tract_manager.generateSpecifiedStreamlineDensityMaps()
    
    tract_manager.debug()
    return

if (__name__ == '__main__'):
    # Create argument parser.
    description = '''Segment tractogram (TRK file) into sets of streamlines terminating at each pair of ROIs defined in a label image (NIFTI file with integer labels). Streamlines passing through a series of ROIs will be converted into a set of streamlines which terminate at each pair of ROIs. Segments may not loop onto terminal ROIs. For example, if a streamline passes through ROIs 1->2->5->1->5, it will be converted into segments (1->2), (1->2->5), (2->5), (2->5->1), (5->1), (1->5).

Available outputs (all optional):
  (1) modified TRK file containing all segmented tracts;
  (2) an ROI x ROI connectivity matrix based on number of streamlines
  (3) ROI x ROI connectivity matrices based on the mean value of a metric (e.g. fractional anisotropy) defined in a separate input metric image, weighted by streamline density
  (4) TRK files including only tract segments terminating at a specified pairs of ROIs
  (5) streamline density maps (NIFTI image) indicating the number of streamlines at each voxel for sets of tract segments terminating at a specified pair of ROIs'''

    class formatter_class(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    parser = argparse.ArgumentParser(description=description, formatter_class=formatter_class)
    
    # Define positional arguments.
    parser.add_argument('trk_path', help='path to tractogram (TRK file).', type=str)
    parser.add_argument('label_path', help='path to label image (NIFTI file with integer labels)', type=str)
    
    # Define optional arguments.
    parser.add_argument('--skip_segmentation', help='do not segment the streamlines (assume input TRK file is already segmented). Useful if you have already segmented the streamlines in an earlier run.', action='store_true')
    parser.add_argument('--include_terminal_rois', help='retain only streamline segments which terminate at these ROIs. Use "--include_terminal_rois", "--exclude_terminal_rois", or neither, but not both; if neither option is used, segments terminating at any ROI will be retained.', type=int, nargs='+')
    parser.add_argument('--exclude_terminal_rois', help='retain all streamline segments except those which terminate at these ROIs; if this option is not used. Use "--include_terminal_rois", "--exclude_terminal_rois", or neither, but not both; if neither option is used, segments terminating at any ROI will be retained.', type=int, nargs='+')
    parser.add_argument('--exclude_terminal_roi_interiors', help='By default, when segmenting a streamline to terminate at a pair of ROIs, the portion of the streamline inside the terminal ROIs is included in the segment. Use this option to exclude portions of the streamline segment inside the terminal ROIs. Note: regardless of whether this option is used, streamline segments are always allowed to pass through ROIs other than the terminal ROIs.', action='store_true')
    parser.add_argument('--out_segmented_trk_path', help='path to output TRK file will all segmented streamlines', type=str)
    parser.add_argument('--out_nos_matrix_path', help='path to output connectivity matrix based on the number of streamlines terminating at each pair of ROIs', type=str)
    parser.add_argument('--in_metric_paths', help='paths to input metric NIFTI image (e.g. fractional anisotropy); the mean value of each metric weighted by streamline density will be computed in the segmented tracts and used to generate a separate connectivity matrix for each metric.', nargs='+', type=str)
    parser.add_argument('--out_metric_matrix_paths', help='paths to output connectivity matrices based on the mean value of metrics in each segmented tract. The number of --out_metric_matrix_paths entries must match the number of --in_metric_paths entries.', nargs='+')
    parser.add_argument('--save_roi_pairs', help='pairs of ROIs for which separate TRK and streamline density maps will be saved, depending on the --out_roi_pair_trk_paths and --out_roi_pair_density_paths options. Each pair of ROIs must be specified by a comma-separated pair of integers, and each pair of ROIs must be separated by a space. For example, to specify the ROI pairs 1<->2 and 1<->3, use "--save_roi_pairs 1,2 1,3".', nargs='+')
    parser.add_argument('--out_roi_pair_trk_paths', help='paths to output TRK files including only streamline segments terminating at a specific pair of ROIs; the number of output paths must match the number of "--save_roi_pairs" entries', nargs='+')
    parser.add_argument('--out_roi_pair_density_paths', help='paths to output streamline density maps (NIFTI files) including only streamline segments terminating at a specific pair of ROIs; the number of output paths must match the number of "--save_roi_pairs" entries', nargs='+')
    parser.add_argument('-v', '--verbose', help='print progress to console', action='store_true')
    
    # Print help if no arguments input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    main(**vars(args))

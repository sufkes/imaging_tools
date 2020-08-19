#!/usr/bin/env python

import os, sys
import argparse
import glob
import numpy as np
from collections import OrderedDict
from pprint import pprint
from distanceCorrelation import distanceCorrelation, normalizeTimeSeries

def getRoiName(path):
    # Get the ROI name from the time series path.
    # e.g.  /path/to/time_series/9283U039_12.txt -> 12
    return os.path.basename(path).split("_")[1].split(".")[0]
    
class Subject(object):
    def __init__(self, subject_id):
        self.subject_id = subject_id
        self.time_series_paths = []
        self.time_series = OrderedDict() # key is ROI name or number; value is numpy array with time series for each voxel in a separate row.
        
    def addPath(self, path):
        self.time_series_paths.append(path)

    def sortAndCountPaths(self):
        self.time_series_paths.sort(key=lambda x:int(getRoiName(x)))
        self.num_rois = len(self.time_series_paths)

    def loadTimeSeries(self):
        for path in self.time_series_paths:
            roi_name = getRoiName(path)
            time_series = np.loadtxt(path, skiprows=3)
            self.time_series[roi_name] = time_series

    def makeDistanceCorrelationMatrix(self):
        dcor_matrix = np.identity(self.num_rois)
        for roi_name_i, time_series_i in self.time_series.iteritems():
            roi_index_i = int(roi_name_i) - 1
            for roi_index_j in range(roi_index_i+1, self.num_rois):
                roi_name_j = str(roi_index_j + 1)
                time_series_i = self.time_series[roi_name_i]
                time_series_j = self.time_series[roi_name_j]
                print "Computing distance correlation between ROIs: ("+roi_name_i+", "+roi_name_j+")."
                distance_correlation = distanceCorrelation(time_series_i, time_series_j)
                print distance_correlation
                dcor_matrix[roi_index_i, roi_index_j] = distance_correlation
                dcor_matrix[roi_index_j, roi_index_i] = distance_correlation
            print roi_name_i
            
        
def makeDistanceCorrelationMatrices(in_dir):
    files = glob.glob(os.path.join(in_dir,"*"))
    files.sort()

    subjects = {}
    
    for path in files:
        subject_id = os.path.basename(path).split("_")[0]
#        roi_index = os.path.basename(path).split("_")[1].split(".")[0]
        roi_name = getRoiName(path)

        # Add instance of Subject object for current subject to dict of all subjects.
        if (not subject_id in subjects):
            subjects[subject_id] = Subject(subject_id)

        subjects[subject_id].addPath(path)
    for subject_id, subject in subjects.iteritems():
        num_rois = set()
        subject.sortAndCountPaths()
        num_rois.add(subject.num_rois)
        #pprint(subject.time_series_paths)

    # Verify that all subjects have time series for the same number of ROIs.
    if (len(num_rois) != 1):
        raise Exception("Subjects must all have time series for the same number of ROIs")
    else:
        num_rois = list(num_rois)[0]
        print "Number of ROIs:", num_rois    

    # Load the time series as numpy arrays.
    for subject_id, subject in subjects.iteritems():
        subject.loadTimeSeries()
        
    # Generate a distance correlation matrix for each subject.
    for subject_id, subject in subjects.iteritems():
        pass
    return

if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate distance correlation matrices for time series extracted from ROIs. 

Each time series will be separately variance normalized, as in L. Geerligs et al. (2016).

Time series should all be in same folder, with names "<subject_ID>_<ROI_number>.txt". This is the format generated using getVoxelTimeSeries.sh. Each time series text file should contain one column for each voxel. The first 3 rows give the (coordinates?) of the voxel, and the remaining rows give the intesity values. This is the format output by the FSL function fslmeants with the '--showall' option.."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_dir", help="path to directory containing time series text files.")
    
    # Define optional arguments.
#    parser.add_argument("-n", "--name", help="Subject ID to set PatientName tags to", type=str)

    # Parse arguments.
    args = parser.parse_args()

    makeDistanceCorrelationMatrices(args.in_dir)

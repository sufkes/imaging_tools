#!/usr/bin/env python

import os, sys
import argparse
import glob
import numpy as np
np.seterr(all='raise') # raise exceptions instead of warnings.
#from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
from pprint import pprint
from distanceCorrelation import distanceCorrelation

import scipy.stats # for pearsonr()

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import misc
#from misc.Timer import Timer
from misc.ProgressBar import ProgressBar

def getRoiName(path):
    # Get the ROI name from the time series path.
    # e.g.  /path/to/time_series/9283U039_12.txt -> 12
    return os.path.basename(path).split("_")[1].split(".")[0]
    
class Subject(object):
    roi_sizes = OrderedDict()
    
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
            if (not roi_name in Subject.roi_sizes):
                Subject.roi_sizes[roi_name] = [time_series.shape[1], 0, time_series.shape[1]] # store [original ROI size, number of voxels removed, final ROI size]

    def normalizeTimeSeries(self):
        for roi_name, time_series in self.time_series.iteritems():
            self.time_series[roi_name] = (time_series - time_series.mean(axis=0))/time_series.std(axis=0, ddof=1)
            #print "Error while normalizing time series for (subject, ROI):", self.subject_id, roi_name
            
    def makeDistanceCorrelationMatrix(self):
        pearson_matrix = np.identity(self.num_rois)
        dcor_matrix = np.identity(self.num_rois)
        for roi_name_i, time_series_i in self.time_series.iteritems():
            roi_index_i = int(roi_name_i) - 1
            for roi_index_j in range(roi_index_i+1, self.num_rois):
                roi_name_j = str(roi_index_j + 1)
                time_series_i = self.time_series[roi_name_i]
                time_series_j = self.time_series[roi_name_j]
                #print "Computing distance correlation between ROIs: ("+roi_name_i+", "+roi_name_j+")."
                distance_correlation_squared = distanceCorrelation(time_series_i, time_series_j)[0]
                geerligs_measure = np.sqrt(max(0.0, distance_correlation_squared))
                dcor_matrix[roi_index_i, roi_index_j] = geerligs_measure
                dcor_matrix[roi_index_j, roi_index_i] = geerligs_measure

                # Also compute the Pearson correlation between the mean time series.
                time_series_mean_i = time_series_i.mean(axis=1)
                time_series_mean_j = time_series_j.mean(axis=1)
                pearson_r = scipy.stats.pearsonr(time_series_mean_i,time_series_mean_j)[0]
                pearson_matrix[roi_index_i, roi_index_j] = pearson_r
                pearson_matrix[roi_index_j, roi_index_i] = pearson_r
        self.dcor_matrix = dcor_matrix
        self.pearson_matrix = pearson_matrix

    def saveDistanceCorrelationMatrix(self, out_dir):
        out_filename = self.subject_id + ".txt"
        out_path = os.path.join(out_dir, out_filename)
        np.savetxt(out_path, self.dcor_matrix)

        out_filename = self.subject_id + "_pearson.txt"
        out_path = os.path.join(out_dir, out_filename)
        np.savetxt(out_path, self.pearson_matrix)

    def plotDistanceCorrelationMatrix(self, out_dir):
        out_filename = self.subject_id + ".png"
        out_path = os.path.join(out_dir, out_filename)
        plt.imshow(self.dcor_matrix, cmap='hot', vmin=0.0, vmax=1.0)
        plt.colorbar()
        plt.savefig(out_path)
        plt.close()
        
def makeDistanceCorrelationMatrices(in_dir, out_dir, skip_normalization=False):
    # Get list of time series files.
    files = glob.glob(os.path.join(in_dir,"*"))
    files.sort()

    # Build dict of subjects
    subjects = OrderedDict()
    for path in files:
        subject_id = os.path.basename(path).split("_")[0]
        roi_name = getRoiName(path)

        # Add instance of Subject object for current subject to dict of all subjects.
        if (not subject_id in subjects):
            subjects[subject_id] = Subject(subject_id)
        subjects[subject_id].addPath(path)
        
    for subject_id, subject in subjects.iteritems():
        num_rois = set()
        subject.sortAndCountPaths()
        num_rois.add(subject.num_rois)

    # Count subjects
    num_subjects = len(subjects)
        
    # Verify that all subjects have time series for the same number of ROIs.
    if (len(num_rois) != 1):
        raise Exception("Subjects must all have time series for the same number of ROIs")
    else:
        num_rois = list(num_rois)[0]
        #print "Number of ROIs:", num_rois    

    # Load the time series as numpy arrays.
    p=ProgressBar("Loading time series into NumPy arrays.")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.loadTimeSeries()
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()

    # Restrict each ROI to a region in which voxels are non-zero (i.e. not masked out) for all subjects.
    valid_voxels = OrderedDict()
    #print Subject.roi_sizes
    for roi_name, roi_size_list in Subject.roi_sizes.iteritems():
        roi_size_orig = roi_size_list[0]
        valid_voxels[roi_name] = np.zeros((num_subjects, roi_size_orig), dtype=bool)
        
        # Find the nonzero voxels for each subject.
        i=0
        for subject_id, subject in subjects.iteritems():
            time_series = subject.time_series[roi_name]
            valid_voxels[roi_name][i, :] = ~(time_series == 0.0).all(axis=0)
            i+=1

        # Find the voxels which are nonzero for all subjects.
        valid_voxels[roi_name] = valid_voxels[roi_name].all(axis=0)

        # Record the number of voxels removed.
        roi_size_list[2] = np.count_nonzero(valid_voxels[roi_name])
        roi_size_list[1] = roi_size_list[0] - roi_size_list[2]

        # Remove the voxels from each subject's time series.
        for subject_id, subject in subjects.iteritems():
            time_series = subject.time_series[roi_name]
            nonzero_time_series = time_series[:, valid_voxels[roi_name]]
            subject.time_series[roi_name] = nonzero_time_series

    print "Voxels removed from ROIs [# original, # removed, # final]:"
    for roi_name, roi_size_list in Subject.roi_sizes.iteritems():
        print roi_name, roi_size_list
    
    # Variance normalize the time series.
    if (not skip_normalization):
        p=ProgressBar("Variance normalizing time series.")
        n=0
        for subject_id, subject in subjects.iteritems():
            subject.normalizeTimeSeries()
            n+=1
            p.update(float(n)/float(num_subjects))
        p.stop()
        
    # Generate a distance correlation matrix for each subject.
    p=ProgressBar("Generating distance correlation matrices.")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.makeDistanceCorrelationMatrix()
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()
        
    # Save distance correlation matrix for each subject.
    p=ProgressBar("Save distance correlation matrices to text.")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.saveDistanceCorrelationMatrix(out_dir)
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()
        
    # Save distance correlation matrix plots each subject.
    p=ProgressBar("Plot distance correlation matrices.")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.plotDistanceCorrelationMatrix(out_dir)
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()
    return


if (__name__ == '__main__'):
    # Create argument parser
    description = """Generate distance correlation matrices for time series extracted from ROIs. 

Each time series will be separately variance normalized, as in L. Geerligs et al. (2016).

Time series should all be in same folder, with names "<subject_ID>_<ROI_number>.txt". This is the format generated using getVoxelTimeSeries.sh. Each time series text file should contain one column for each voxel. The first 3 rows give the (coordinates?) of the voxel, and the remaining rows give the intesity values. This is the format output by the FSL function fslmeants with the '--showall' option.."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_dir", help="path to directory containing time series text files.")
    parser.add_argument("out_dir", help="path to directory to which correlation matrices will be saved.")
    
    # Define optional arguments.
    parser.add_argument("-n", "--skip_normalization", help="skip variance normalization of time series if already done", action="store_true")

    # Parse arguments.
    args = parser.parse_args()

    makeDistanceCorrelationMatrices(args.in_dir, args.out_dir, skip_normalization=args.skip_normalization)

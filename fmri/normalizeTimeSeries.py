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
from distanceCorrelation import distanceCorrelation, normalizeTimeSeries

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import misc
#from misc.Timer import Timer
from misc.ProgressBar import ProgressBar

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

    def normalizeTimeSeries(self):
        for roi_name, time_series in self.time_series.iteritems():
            try:
                self.time_series[roi_name] = (time_series - time_series.mean(axis=0))/time_series.std(axis=0, ddof=1)
            except:
                print "Error normalizing time series: (subject, ROI, voxels in ROI, zero voxels in ROI):", self.subject_id, roi_name, time_series.shape[1], list(time_series.std(axis=0, ddof=1) <= 0.0).count(True)
#                print "Number of voxels with zero variance:", list(time_series.std(axis=0, ddof=1) <= 0.0).count(True)
            
    def saveTimeSeries(self, out_dir):
        for roi_name, time_series in self.time_series.iteritems():
            out_filename = self.subject_id + "_" + roi_name + ".txt"
            out_path = os.path.join(out_dir, out_filename)
            np.savetxt(out_path, self.time_series)

def normalizeTimeSeries(in_dir, out_dir):
    # Get list of time series files.
    files = glob.glob(os.path.join(in_dir,"*"))
    files.sort()

    # Build dict of subjects
    subjects = {}
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
    p=ProgressBar("Loading time series into NumPy arrays and variance normalizing")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.loadTimeSeries()
        subject.normalizeTimeSeries()
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()

    # Variance normalize the time series.
#    p=ProgressBar("Variance normalizing time series")
#    n=0
#    for subject_id, subject in subjects.iteritems():
#        subject.normalizeTimeSeries()
#        n+=1
#        p.update(float(n)/float(num_subjects))
#    p.stop()
        
    # Save normalized time series for each subject.
    p=ProgressBar("Save time series to text")
    n=0
    for subject_id, subject in subjects.iteritems():
        subject.saveTimeSeries(out_dir)
        n+=1
        p.update(float(n)/float(num_subjects))
    p.stop()        
    return


if (__name__ == '__main__'):
    # Create argument parser
    description = """Normalize time series extracted from ROIs. 

Each time series will be separately variance normalized, as in L. Geerligs et al. (2016).

Time series should all be in same folder, with names "<subject_ID>_<ROI_number>.txt". This is the format generated using getVoxelTimeSeries.sh. Each time series text file should contain one column for each voxel. The first 3 rows give the (coordinates?) of the voxel, and the remaining rows give the intesity values. This is the format output by the FSL function fslmeants with the '--showall' option.."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("in_dir", help="path to directory containing time series text files.")
    parser.add_argument("out_dir", help="path to directory to which correlation matrices will be saved.")
    

    # Parse arguments.
    args = parser.parse_args()

    normalizeTimeSeries(args.in_dir, args.out_dir)

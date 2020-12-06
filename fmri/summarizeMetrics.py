#!/usr/bin/env python

import glob
import os, sys
import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from pprint import pprint

class MetricNifti(object):
    def __init__(self, path, subject):
        self.path = path
        self.name = "_".join(os.path.basename(path).split('.')[0].split("_")[1:])
        #self.subject_name = os.path.basename(path).split("_")[0]
        self.subject = subject
        self.nii = nib.load(path)
        self.arr = self.nii.get_fdata()
        
class MetricNumpy(object):
    def __init__(self, path, subject):
        self.path = path
        self.name = "_".join(os.path.basename(path).split('.')[0].split("_")[1:])
        #self.subject_name = os.path.basename(path).split("_")[0]
        self.subject = subject
        self.arr = np.load(path)

class RegionOfInterest(object):
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path).split('.')[0]
        self.arr = nib.load(path).get_fdata()

class Mask(object):
    def __init__(self, path, mask_thres=0):
        self.path = path
        arr = nib.load(path).get_fdata()
        #print arr.shape
        print "Number of voxels in mask before thresholding (val > "+str(mask_thres)+"):"
        print np.count_nonzero(arr)
        arr[arr<=mask_thres] = 0.0
        print "Number of voxels in mask after thresholding (val > "+str(mask_thres)+"):"
        print np.count_nonzero(arr)
        #print "Shape of mask before np.squeeze:"
        #print arr.shape
        #print arr.shape
        #arr = np.squeeze(arr)
        #print "Shape of mask after np.squeeze:"
        #print arr.shape
        self.arr = arr
        
def summarizeMetrics(metric_dir, roi_dir, mask_path, mask_thres, out_dir):    
    # Initialize the dataframe which will store all of the data.
    df = pd.DataFrame()

    # Get the mask file
    mask = Mask(mask_path, mask_thres=mask_thres)
        
    # Find the ROI NIFTI files.
    rois = []
    
    for root, directories, file_names in os.walk(roi_dir):
#        for directory in directories:
#            print os.path.join(root, directory)
        for file_name in file_names:
            if ".nii" in file_name:
                file_path = os.path.join(root, file_name)
                roi = RegionOfInterest(file_path)
                rois.append(roi)
    rois.sort(key=lambda x: x.path)
    #pprint([(roi.path, roi.name) for roi in rois])
    
    
    # Find subject dirs. Assume that the name of the directory is the subject name.
    metric_dir = metric_dir.rstrip('/') # trailing slash seems to trick os.path.basename()
    metrics_nifti_filter = ['betweenness_w', 'charpath_ecc_w', 'clustering_coef_w', 'degree_w', 'local_efficiency_w', 'fcs', 'reho']
    metrics_npy_filter = ['betweenness_w', 'charpath_diameter_w', 'charpath_ecc_w', 'charpath_efficiency_w', 'charpath_lambda_w', 'charpath_radius_w', 'clustering_coef_w', 'corr', 'degree_w', 'fcs', 'global_efficiency_w', 'local_efficiency_w', 'reho', 'transitivity_w']
    
    subject_dirs = glob.glob(metric_dir+'/*')
    subject_dirs.sort()
    
    sane = False # Not sure if sane yet. Let's assume not. 

    volume_df = pd.DataFrame(index=[roi.name for roi in rois], columns=["num_voxels"]) # to store the volumes of the ROI/mask overlap.
    for roi in rois:
        volume_df.loc[roi.name, "num_voxels"] = len(roi.arr[(mask.arr>0.0) & (roi.arr>0.0)])

    out_name = os.path.basename(metric_dir) + "_roi_mask_overlap_volumes.csv"
    out_path = os.path.join(out_dir, out_name)
    volume_df.to_csv(out_path)
            
    
    # Add row for each subject.
    for subject_dir in subject_dirs:
        subject = os.path.basename(subject_dir)
        print subject
        subject_dir_abs = os.path.abspath(subject_dir)
        df = df.append(pd.Series(name=subject)) # set index to subject ID.
        
        df.loc[subject, "metric_dir"] = subject_dir_abs

        ## NumPy array metrics
        # Find the metric NumPy array files. Some of these are voxel-wise, some are not.
        metric_paths_npy = glob.glob(subject_dir+"/*.npy")
        #pprint(metric_paths_npy)
        
        # Exclude any npy files which are not in the list of metrics specified above.
        metrics_npy = []
        for metric_path in metric_paths_npy:
            for string in metrics_npy_filter:
                if string in metric_path:
                    metric = MetricNumpy(metric_path, subject)
                    metrics_npy.append(metric)
                    break

        metrics_npy.sort(key=lambda x: x.path)
        #pprint(metrics_npy)
        
        for metric in metrics_npy:
            if (len(metric.arr.shape) > 0):
                col_prefix = metric.name
                col_mean = col_prefix+"_mean"
                col_median = col_prefix+"_median"
                col_std = col_prefix+"_std"

                if not col_mean in df:
                    df[col_mean] = 0.0
                if not col_median in df:
                    df[col_median] = 0.0
                if not col_std in df:
                    df[col_std] = 0.0
            
                df.loc[subject, col_mean] = metric.arr.mean()
                df.loc[subject, col_median] = np.median(metric.arr)
                df.loc[subject, col_std] = metric.arr.std(ddof=1)
            else:
                if not metric.name in df:
                    df[metric.name] = 0.0
                df.loc[subject, metric.name] = metric.arr
            
        ## NIFTI voxel-wise metrics
        # Find the metric NIFTI files.
        metric_paths_nifti = glob.glob(subject_dir+"/*.nii.gz")

        # Exclude any NIFTI files which are not in the list of metrics specified above.
        metrics_nifti = []
        for metric_path in metric_paths_nifti:
            for string in metrics_nifti_filter:
                if string in metric_path:
                    metric = MetricNifti(metric_path, subject)
                    metrics_nifti.append(metric)
                    break

        metrics_nifti.sort(key=lambda x: x.path)

        for metric in metrics_nifti:
            # Get the mean and standard deviation over the whole grey matter mask. This will be redundant for most measures, whose GM-averaged values are calculated above from the Numpy arrays. This is not redundant for FCS, which does not have a numpy array version.
            if (metric.name == 'fcs'):
                col_prefix = metric.name + "_gm_"
                col_mean = col_prefix+"_mean"
                col_median = col_prefix+"_median"
                col_std = col_prefix+"_std"
                if not col_mean in df:
                    df[col_mean] = 0.0
                if not col_median in df:
                    df[col_median] = 0.0
                if not col_std in df:
                    df[col_std] = 0.0

                # Get the portion of the metric lying within the mask.
                metric_in_mask = metric.arr[mask.arr > 0.0] # This is flat.

                df.loc[subject, col_mean] = metric_in_mask.mean()
                df.loc[subject, col_median] = np.median(metric_in_mask)
                df.loc[subject, col_std] = metric_in_mask.std(ddof=1)
                    
            for roi in rois:
                col_prefix = metric.name + "_" + roi.name
                col_mean = col_prefix+"_mean"
                col_median = col_prefix+"_median"
                col_std = col_prefix+"_std"
                if not col_mean in df:
                    df[col_mean] = 0.0
                if not col_median in df:
                    df[col_median] = 0.0
                if not col_std in df:
                    df[col_std] = 0.0

                # Get the portion of the metric lying within both the mask and the ROI.
                metric_in_roi_and_mask = metric.arr[((roi.arr > 0.0) & (mask.arr > 0.0))] # This is flat.

                df.loc[subject, col_mean] = metric_in_roi_and_mask.mean()
                df.loc[subject, col_median] = np.median(metric_in_roi_and_mask)
                df.loc[subject, col_std] = metric_in_roi_and_mask.std(ddof=1)

                if not sane:
                    # Save a NIFTI file of the metric limited to both the mask and the ROI, as a sanity check. Also save an NIFTI file of ones limited the same way as a check.
                    print "Saving sanity check images for", metric.subject, metric.name, roi.name
                    metric_lim_arr = metric.arr.copy()
                    metric_lim_arr[~((roi.arr > 0.0) & (mask.arr > 0.0))] = 0.0
                    metric_lim_nii = nib.Nifti1Image(metric_lim_arr, metric.nii.affine, metric.nii.header)
                    metric_lim_path = os.path.join(out_dir, "-".join([metric.subject, metric.name, roi.name]) + ".nii.gz")
                    nib.save(metric_lim_nii, metric_lim_path)
                    
                    test_ones_arr = np.ones(metric.arr.shape)
                    test_ones_arr[~((roi.arr > 0.0) & (mask.arr > 0.0))] = 0.0
                    test_ones_nii = nib.Nifti1Image(test_ones_arr, metric.nii.affine, metric.nii.header)
                    test_ones_path = os.path.join(out_dir, "-".join(['ones', roi.name]) + ".nii.gz")
                    nib.save(test_ones_nii, test_ones_path)
                    sane = True
    
    # Save the dataframe.
    out_name = os.path.basename(metric_dir) + "_data.csv"
    out_path = os.path.join(out_dir, out_name)
    df.to_csv(out_path, index=True, encoding='utf-8')
    return
    
if (__name__ == '__main__'):
    # Create argument parser
    description = """Calculate mean and standard deviation of voxel-wise metrics within regions of interest. Input metrics and ROIs must be in NIFTI format. Additionally, summarize numpy arrays (maybe)."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument("metric_dir", type=str, help="directory in which metric .nii.gz and .npy files live.")
    parser.add_argument("roi_dir", type=str,  help="directory in which region of interest (ROI) .nii.gz files live. Mean and std. dev. values of each metric will be calculated within these ROIs.")
    parser.add_argument("mask_path", type=str, help="path to NIFTI mask which defines a restricted region within which the metrics were calculated")
    
    # Define optional arguments.
    parser.add_argument("-o", "--out_dir", type=str, help="output directory. Default: current directory.", default=os.getcwd())
    parser.add_argument("--mask_thres", type=float, help="mask threshold -- only consider metric values within the portion of the mask whose values are above this threshold. By default, all positive values in the mask are included.", default=0)

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run the
    summarizeMetrics(**vars(args))

    

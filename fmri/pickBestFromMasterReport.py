#!/usr/bin/env python

import os, sys
from collections import OrderedDict
import numpy as np

import pandas as pd

with open(str(sys.argv[1])) as handle:
    lines = handle.readlines()

# Convert the report lines to a Pandas DataFrame
best_array = np.zeros(shape=(120*121,4),dtype='object')
for line_index in range(len(lines)):
    line = lines[line_index]
    subject = str(line.split()[0])
    refvol = int(line.split()[1])
    abs_mean_rms = float(line.split()[2])
    status = line.split()[3]
    if (status == 'Pass'):
        status = 1
    else:
        status = 0

    best_array[line_index, 0] = subject
    best_array[line_index, 1] = refvol
    best_array[line_index, 2] = abs_mean_rms
    best_array[line_index, 3] = status    
best = pd.DataFrame(best_array, columns=['subject', 'refvol', 'abs_mean_rms', 'status'])

# Pick the best scan for each subject.
for subject in best['subject'].unique():
    if (len(best.loc[(best['subject'] == subject) & (best['status'] == 1), :]) == 0): # If the motion threshold was exceeded for all reference volumes, find min abs_mean_rms over all refvols
        subject_rows = best.loc[(best['subject'] == subject), :]
        min_index = subject_rows.loc[subject_rows['abs_mean_rms'] == subject_rows['abs_mean_rms'].min(), :].index[0] # index of refvol with smallest abs_mean_rms for subject
        min_refvol = best.loc[min_index, 'refvol']
        abs_mean_rms_of_refvol = best.loc[min_index, 'abs_mean_rms']

        print subject, min_refvol, abs_mean_rms_of_refvol, 'Fail'
        
    else: # If the motion threshold is not exceeded for at least one reference volume, find min abs_mean_rms over all passing refvols
        subject_rows = best.loc[(best['subject'] == subject) & (best['status'] == 1), :]
        min_index = subject_rows.loc[subject_rows['abs_mean_rms'] == subject_rows['abs_mean_rms'].min(), :].index[0] # index of passing refvol with smallest abs_mean_rms for subject
        min_refvol = best.loc[min_index, 'refvol']
        abs_mean_rms_of_refvol = best.loc[min_index, 'abs_mean_rms']
        
        print subject, min_refvol, abs_mean_rms_of_refvol, 'Pass'

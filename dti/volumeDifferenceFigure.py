#!/usr/bin/env python

import os, sys 
from matplotlib import pyplot as plt
import numpy as np

## script must be run in /hpf/largeprojects/smiller/scratch/sufkes/tbss_fig3_volume_difference_distribution/

def makePlot(path1, path2, out_path):
    # Read text file storing (target_volume, warped_image_volume) for all subject-subject registrations (excluding self-self registrations).
    path_list = [path1, path2]
    
    lines_list = []
    for path in path_list:
        with open(path, 'rb') as handle:
            lines_list.append(handle.readlines())

    # Get volume differences
    deltas_list = [np.array([],dtype=np.float)]*2
#    for lines in lines_list:
    for ii in range(2):
        lines = lines_list[ii]
        for line in lines:
            tar_vol = np.float(line.split()[0])
            wrp_vol = np.float(line.split()[1])
            delta = wrp_vol - tar_vol

            # convert from mm^3 to cm^3
            delta = delta/1000

            # Add difference to array of differences
            deltas_list[ii] = np.append(deltas_list[ii], delta)

            # Convert differences to absolute differences
            deltas_list[ii] = np.abs(deltas_list[ii])


    x = deltas_list[0]
    y = deltas_list[1]

    mu_x = np.mean(x)
    sigma_x = np.std(x)
    median_x = np.median(x)
    x_25th = np.percentile(x, 25)
    x_75th = np.percentile(x, 75)
    mu_y = np.mean(y)
    sigma_y = np.std(y)
    median_y = np.median(y)
    y_25th = np.percentile(y, 25)
    y_75th = np.percentile(y, 75)

    n_y, bins_y, patches_y = plt.hist(y, bins='auto', alpha=0.5, label=r'6 DoF and 12 DoF', density=True, color='red')
    n_x, bins_x, patches_x = plt.hist(x, bins=bins_y, alpha=0.5, label=r'12 DoF', density=True, color='blue')
    plt.legend(loc='upper right')

    # the histogram of the data

    # add a 'best fit' line
    x_line = ((1 / (np.sqrt(2 * np.pi) * sigma_x)) * np.exp(-0.5 * (1 / sigma_x * (bins_x - mu_x))**2))
    plt.plot(bins_x, x_line, '--', color='blue')

    y_line = ((1 / (np.sqrt(2 * np.pi) * sigma_y)) * np.exp(-0.5 * (1 / sigma_y * (bins_y - mu_y))**2))
    plt.plot(bins_y, y_line, '-', color='red')

    plt.xlabel(r'Volume difference (cm$^3$)')
    plt.ylabel('Frequency')

    plt.text(median_x, np.max(x_line)*1.05, r'median: '+"{0:.1f}".format(median_x)+' cm$^3$, IQR: '+"{0:.1f}".format(x_25th)+r'$-$'+"{0:.1f}".format(x_75th)+' cm$^3$')
    plt.text(median_y, np.max(y_line)*1.05, r'median: '+"{0:.1f}".format(median_y)+' cm$^3$, IQR: '+"{0:.1f}".format(y_25th)+'$-$'+"{0:.1f}".format(y_75th)+' cm$^3$')

    plt.savefig(out_path)
    plt.close()

for group in ['group1', 'group2', 'group3', 'group4']:
#for group in ['group1']:
    vols_12_path = os.path.join('12', group+'_volumes.txt')
    vols_6_12_path = os.path.join('6_12', group+'_volumes.txt')

    makePlot(vols_12_path, vols_6_12_path, group+'_volume_differences.png')

#!/usr/bin/env python

import os, sys 
from matplotlib import pyplot as plt
import numpy as np

## script must be run in /hpf/largeprojects/smiller/scratch/sufkes/tbss_tests_20191011

# $7      $8     $9      $10    $11     $12    $13         $14            $15
# dex_yes dex_no inf_yes inf_no bpd_yes bpd_no age_atbirth days_intubated

def makePlot(data_path, out_path):
    # Read text file storing (target_volume, warped_image_volume) for all subject-subject registrations (excluding self-self registrations).

    lines = []
    with open(data_path, 'rb') as handle:
        lines = handle.readlines()

    bpd_din_list = []
    for line in lines:
        bpd_din_list.append([line.split()[ii] for ii in [10,13]])
    
    din_bpd_y = np.array([])
    din_bpd_n = np.array([])
    for row in bpd_din_list:
        if (row[0] == '0'):
            din_bpd_n = np.append(din_bpd_n, float(row[1]))
        else:
            din_bpd_y = np.append(din_bpd_y, float(row[1]))

    print din_bpd_y, din_bpd_n

    n_n, bins_n, patches_n = plt.hist(din_bpd_n, bins=range(0, int(din_bpd_n.max()), 1), alpha=0.5, label=r'no BPD', color='blue', range=(0, din_bpd_n.max()))
    n_y, bins_y, patches_y = plt.hist(din_bpd_y, bins=bins_n, alpha=0.5, label=r'BPD', color='red', range=(0, din_bpd_y.max()))
    plt.legend(loc='upper right')

    # the histogram of the data

    # add a 'best fit' line
#    x_line = ((1 / (np.sqrt(2 * np.pi) * sigma_x)) * np.exp(-0.5 * (1 / sigma_x * (bins_x - mu_x))**2))
#    plt.plot(bins_x, x_line, '--', color='blue')

#    y_line = ((1 / (np.sqrt(2 * np.pi) * sigma_y)) * np.exp(-0.5 * (1 / sigma_y * (bins_y - mu_y))**2))
#    plt.plot(bins_y, y_line, '-', color='red')

    plt.xlabel(r'Days intubated')
    plt.ylabel('Number of patients')

#    plt.text(median_x, np.max(x_line)*1.05, r'median: '+"{0:.1f}".format(median_x)+' cm$^3$, IQR: '+"{0:.1f}".format(x_25th)+r'$-$'+"{0:.1f}".format(x_75th)+' cm$^3$')
#    plt.text(median_y, np.max(y_line)*1.05, r'median: '+"{0:.1f}".format(median_y)+' cm$^3$, IQR: '+"{0:.1f}".format(y_25th)+'$-$'+"{0:.1f}".format(y_75th)+' cm$^3$')

    plt.savefig(out_path)
    plt.close()    

for group in ['group1', 'group2', 'group3', 'group4']:
    data_path = os.path.join(group, 'data_'+group+'.txt')
    out_path = os.path.join(os.getcwd(), group+'_days_intubated_histogram.png')

    makePlot(data_path, out_path)

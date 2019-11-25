#!/usr/bin/env python

import os, sys
from nilearn import plotting as plt

p_path = str(sys.argv[1]) # 1-p image
p_name = os.path.basename(p_path)
stats_dir = os.path.dirname(p_path) # Assume the directory containing the 1-p image also contains the mean_FA and mean_FA_skeleton_mask images.
mean_FA_path = os.path.join(stats_dir, 'mean_FA.nii.gz')
mean_FA_skeleton_mask_path = os.path.join(stats_dir, 'mean_FA_skeleton_mask.nii.gz')

#display = plt.plot_img(mean_FA_path, display_mode='z', cut_coords=8, cmap='gray') # automatically generate cuts (number of cuts = cut_coods)
display = plt.plot_img(mean_FA_path, display_mode='z', cut_coords=(10,), cmap='gray') # Generate cuts at specified points.
display.add_overlay(mean_FA_skeleton_mask_path, cmap='ocean')
display.add_overlay(p_path, cmap='YlOrRd_r',threshold=0.95,colorbar=True)

#out_path = os.path.join(stats_dir, p_name.split(".nii")[0]+'.png') # automatically generate name for image, and save to directory where the (1-p) is stored.
out_path = str(sys.argv[2])
display.savefig(out_path)

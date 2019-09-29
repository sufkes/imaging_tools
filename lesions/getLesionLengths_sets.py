#!/usr/bin/env python
"""This script takes a minc file which maps out lesions, groups the voxels 
into sets which corresponds to, e.g., a lesion, determines the maximum length
of the volume (i.e., the length of the volume measured along whichever 
direction, not necessarily along the coordinate axis, gives the longest length),
and returns the length of each lesion, the average coordinate of the lesion, 
and the length of the longest lesion.

Usage:
  python getLesionLengths.py <minc file>

Notes: 
- This code assumes that each voxel is a cube in the real world.
- The lesion minc file must be all zeros except for the lesion voxels, which 
must have an intensity value greater than zero.
- This code ignore the naming of the x, y, z axes, and refers to the 1st, 2nd, 
and 3rd dimensions in the exported volumne as x, y, and z respectively.
# THERE APPEARS TO BE A GENUINE ERROR IN MINCINFO WITH REGARD TO THE DIMENSION REPORT!
"""

# Standard modules
import os, sys
#import argparse

# Non-standard modules
import numpy as np
import pyminc.volumes.factory as pyminc

# My modules in other directories
sufkes_git_repo_dir = "/Users/steven ufkes/scripts" # change this to the path to which the sufkes Git repository was cloned.
sys.path.append(os.path.join(sufkes_git_repo_dir, "misc"))
from Color import Color


minc_file_path = str(sys.argv[1])
minc_vol = pyminc.volumeFromFile(minc_file_path, labels=True)
separations = minc_vol.separations
volume = minc_vol.data

def getVolumeSet(volume, val_list):
    """Separate a single volume into multiple volumes of the same size
    by sorting voxels into groups based on the number assigned to them.
    For example, suppose you labelled right-hemisphere lesions with '1' 
    and left-hemisphere lesions with '2', this function will return a
    list of 2 volumes, each with the original size. One will contain all
    zeros except for the voxels labelled '1'; the other will contain all 
    zeros except for the voxels labelled '2'. 
    This function will also convert all the labels to the int type by 
    rounding to the nearest integer."""
    volume_list = [] # list containing volumes extracted from the input volume according to the number assigned to the ROI type.
    
#    volume = np.rint(volume).astype(int) # redundant if using labels=True when loading volume
   
    #print "Separating ROIs with values", val_list, "into separate volumes"

    for val in val_list:
        new_volume = np.copy(volume)
        for ii in range(new_volume.shape[0]):
            for jj in range(new_volume.shape[1]):
                for kk in range(new_volume.shape[2]):
                    if (new_volume[ii, jj, kk] != val):
                        new_volume[ii, jj, kk] = 0
        volume_list.append(new_volume)
    return volume_list

def demarcateLesion(volume, coord): # this function should only be entered once per lesion
    dim_x = np.shape(volume)[0]
    dim_y = np.shape(volume)[1]
    dim_z = np.shape(volume)[2]

    coord_list = [coord] # list containing triplets of lesion coordinates
    seed_list = [coord]
    checked_list = [] # list of voxels which have been checked (>1) and used as a seed.
    while (len(seed_list) > 0):
        seed_coord = seed_list[0]
        neighbors = [] # list of 26 neighboring voxel coordinates.
        seed_ii = seed_coord[0]
        seed_jj = seed_coord[1]
        seed_kk = seed_coord[2]
        for d_ii in range(-1,2):
            for d_jj in range(-1,2):
                for d_kk in range(-1,2):
                    if (0 <= seed_ii+d_ii < dim_x) and (0 <= seed_jj+d_jj < dim_y) and (0 <= seed_kk+d_kk < dim_z): # if neighbor is inside volume
                        neighbors.append((seed_ii+d_ii, seed_jj+d_jj, seed_kk+d_kk))
        neighbors.remove(seed_coord) # the seed is not its own neighbor
        
        for neighbor in neighbors:
            voxel = volume[neighbor]
            if (voxel > 0): # if voxel is in a lesion
                if (not neighbor in checked_list): # if voxel has not been checked
                    if (not neighbor in coord_list):
                        coord_list.append(neighbor)
                    if (not neighbor in seed_list):
                        seed_list.append(neighbor)
        
        # Remove current voxel from seed list, and add it to the checked list.
        seed_list.remove(seed_coord)
        checked_list.append(seed_coord)        
    return coord_list

def findLesions(volume):
    lesions = [] # list containing lesions sets. Each lesion is a list of coordinate triplets.
    lesions_found = False
    for ii in range(np.shape(volume)[0]):
        for jj in range(np.shape(volume)[1]):
            for kk in range(np.shape(volume)[2]):
                coord = (ii,jj,kk)
                voxel = volume[ii,jj,kk]
                if (voxel > 0.0): # if voxel is in a lesion
                    if (lesions_found): # if a lesion has already been identified
                        # Look for the current voxel in known lesions.
                        lesion_known = False
                        for lesion in lesions:
                            if coord in lesion: # if voxel is in a known lesion (i.e. if the lesion has already been identified)
                                lesion_known = True
                                break
                    else: # if no lesions have been identified
                        lesion_known = False 
                        lesions_found = True

                    # Demarcate lesion if it is currently unknown.
                    if (not lesion_known): # if voxel is in an unknown lesion
                        lesion_coord_list = demarcateLesion(volume, coord) # coordinates of all voxels in lesion.
                        lesions.append(lesion_coord_list)
    return lesions

x_sep = np.abs(separations[0])
y_sep = np.abs(separations[1])
z_sep = np.abs(separations[2])
voxel_vol = x_sep*y_sep*z_sep

val_list = [1, 2] # list of values used to label ROI types.

volume_list = getVolumeSet(volume, val_list) # break input volume into separate volumes according to the number assigned to each ROI type.

max_list = [] # List storing tuple for each volume: (label used for volume, length of longest lesion, volume of largest lesion)

for volume_index in range(len(volume_list)):
    volume = volume_list[volume_index]
    label_val = val_list[volume_index]

#    label_val = np.amax(volume)

    #print
    #print "Reporting on lesions labelled with:", label_val

    # Store longest and large lesion sizes in volume:
    len_max = 0.0
    vol_max = 0.0
    vol_tot = 0.0

    lesions = findLesions(volume)
    
    # Calculate total lesion volume
    num_lesions = len(lesions)
    num_lesion_voxels = 0
    for ii in range(np.shape(volume)[0]):
            for jj in range(np.shape(volume)[1]):
                for kk in range(np.shape(volume)[2]):
                    voxel = volume[ii, jj, kk]
                    if (voxel > 0):
                        num_lesion_voxels += 1
    
    # Calculate some figures for each lesion (location, maximum length, volume).
    lesion_stats = []
    for lesion_index in range(num_lesions):
        lesion = lesions[lesion_index]
        lesion_stats.append({})
    
        # Find the volume of the lesion.
        volume = len(lesion)*voxel_vol
        lesion_stats[lesion_index]["volume"] = volume
    
        # Find the length of the lesion along its longest dimension.
        length = 0.0
        length_x = 0.0
        length_y = 0.0
        length_z = 0.0
        for coord_1 in lesion:
            for coord_2 in lesion:
                dist_x = np.abs(float(coord_1[0]) - float(coord_2[0]) + 1)*x_sep
                dist_y = np.abs(float(coord_1[1]) - float(coord_2[1]) + 1)*y_sep
                dist_z = np.abs(float(coord_1[2]) - float(coord_2[2]) + 1)*z_sep
                dist = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)
                if (dist > length):
                    length = dist
                if (dist_x > length_x):
                    length_x = dist_x
                if (dist_y > length_y):
                    length_y = dist_y
                if (dist_z > length_z):
                    length_z = dist_z
        lesion_stats[lesion_index]["length"] = length
        lesion_stats[lesion_index]["length_x"] = length_x
        lesion_stats[lesion_index]["length_y"] = length_y
        lesion_stats[lesion_index]["length_z"] = length_z
    
        # SHOULD I ADD 1/2 TO THIS CENTER OF MASS? MAYBE.
        # Find the center of mass of the lesion.
        x_tot = 0.0
        y_tot = 0.0
        z_tot = 0.0
        for coord in lesion:
            x_tot += coord[0]
            y_tot += coord[1]
            z_tot += coord[2]
        x_av = float(x_tot)/float(len(lesion))
        y_av = float(y_tot)/float(len(lesion))
        z_av = float(z_tot)/float(len(lesion))
        com = (x_av, y_av, z_av)
        lesion_stats[lesion_index]["com"] = com
    
        # Check if longest or largest lesion in volume. 
        if (len_max < length):
            len_max = length
        if (vol_max < volume):
            vol_max = volume

        # Increase total lesion volume.
        vol_tot += volume

    # Write the longest and largest lesion sizes to variables:
    max_list.append((label_val, len_max, vol_max, num_lesions, vol_tot))

    #print "Number of lesions: "+str(num_lesions)
    #print "Total number of lesion voxels: "+str(num_lesion_voxels)
    for lesion_index in range(num_lesions):
        lesion = lesions[lesion_index]
        #print 
        #print Color.bold+"Lesion "+str(lesion_index+1)+" of "+str(num_lesions)+Color.end
        #print "Lesion volume        : "+str(lesion_stats[lesion_index]["volume"])
        #print "Lesion length        : "+str(lesion_stats[lesion_index]["length"])
#        print "Lesion length (x)    : "+str(lesion_stats[lesion_index]["length_x"])
#        print "Lesion length (y)    : "+str(lesion_stats[lesion_index]["length_y"])
#        print "Lesion length (z)    : "+str(lesion_stats[lesion_index]["length_z"])
#        print "Lesion center of mass: "+str(lesion_stats[lesion_index]["com"])

# E.g. minc_file_path = BC0002_V01-CbH3.mnc
id = minc_file_path.split("_")[0].split("BC")[1].lstrip("0")
visit = minc_file_path.split("_")[1].split("-")[0]
len_max_left = max_list[1][1]
len_max_right = max_list[0][1]
max_len = max(len_max_left, len_max_right)
vol_max_left = max_list[1][2]
vol_max_right = max_list[0][2]
max_vol = max(vol_max_left, vol_max_right)
num_left = max_list[1][3]
num_right = max_list[0][3]
vol_tot_left = max_list[1][4]
vol_tot_right = max_list[0][4]
vol_tot = vol_tot_left + vol_tot_right
print id, visit, len_max_left, len_max_right, max_len, vol_max_left, vol_max_right, max_vol, num_left, num_right, vol_tot_left, vol_tot_right, vol_tot



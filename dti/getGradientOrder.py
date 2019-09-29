#!/usr/bin/env python
import os, sys
import pydicom
import numpy as np
import argparse
from collections import OrderedDict

# Create argument parser.
description = '''
Create map of volume number -> (gradient direction, b value) for UBC preterm cohort scans.
'''
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

# Define positional arguments. 
parser.add_argument("in_dir", action='store', type=str, help="path to DTI series to be mapped")

# Define optional arguments.
parser.add_argument('-o', '--out_dir', action='store', type=str, default=os.getcwd(), help="output directory (Default: current working directory)")

# Print help if no args input.
if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit()

# Parse input arguments and store them.
args = parser.parse_args()

# Generate list of all files in directory
file_list = []
for dirname, dirnames, filenames in os.walk(args.in_dir):
    for filename in filenames:
        file_list.append(os.path.abspath(os.path.join(dirname, filename)))

num_files = len(file_list)
print "Creating map for "+str(num_files)+" files."

# Create dict mapping instance numbers to gradient directions, b values, and acqusition times.
file_dict = {}
InstanceNumber_max = 0
zero_based_gradient_index = False # whether first non-zero gradient is numbered '#0' or '#1'
for file_path in file_list:
    ds = pydicom.read_file(file_path, force='force')
    InstanceNumber = int(ds[0x0020, 0x0013].value)
    if (InstanceNumber > InstanceNumber_max):
        InstanceNumber_max = InstanceNumber
    AcquisitionTime = ds[0x0008, 0x0032].value
    SequenceName = ds[0x0018, 0x0024].value
#    print InstanceNumber, type(InstanceNumber) # int (after forcing it to be)
#    print AcquisitionTime, type(AcquisitionTime) # str
#    print SequenceName, type(SequenceName) # str

    # Try to extract gradient direction and b value from SequenceName.
    # Examples of SequenceName values in the UBC preterm cohort are (including the asterisks):
    # *ep_b0
    # *ep_b600#6
    # *ep_b700#12  
    BValue = SequenceName.split('b')[1].split('#')[0]
    if ('#' in SequenceName):
        GradientNumber = int(SequenceName.split('#')[1])
        if (GradientNumber == 0):
#            raise ValueError("Gradient direction number specified as '#0' in SequenceName: '"+SequenceName+"'. Sketchily incrementing gradient number of b>0 volumes by one.")
#            print "Gradient direction number specified as '#0' in SequenceName: '"+SequenceName+"'."
            zero_based_gradient_index = True
    else:
        GradientNumber = 0

    file_dict[InstanceNumber] = {}
    file_dict[InstanceNumber]['AcquisitionTime'] = AcquisitionTime
    file_dict[InstanceNumber]['SequenceName'] = SequenceName
    file_dict[InstanceNumber]['GradientNumber'] = GradientNumber
    file_dict[InstanceNumber]['BValue'] = BValue

# Increment b>0 gradient indices if zero-based gradient direction indexing is used.
if zero_based_gradient_index:
    print "Warning: First b>0 volume gradient direction is labelled as '#0'. Incrementing GradientNumber by 1 for all b>0 volumes."
    for InstanceNumber, file in file_dict.iteritems():
        if ("#" in file['SequenceName']):
            file['GradientNumber'] += 1

# Convert dict to OrderedDict sorted based on instance number.
file_ord_dict = OrderedDict()
for InstanceNumber in range(1, InstanceNumber_max+1):
    try:
        file_ord_dict[InstanceNumber] = file_dict[InstanceNumber]
    except KeyError:
        print "Warning: InstanceNumber "+str(InstanceNumber)+" not found."

#for InstanceNumber, value in file_ord_dict.iteritems():
#    print "InstanceNumber:", InstanceNumber
#    print "  SequenceName   :", file_ord_dict[InstanceNumber]['SequenceName']
#    print "  AcquisitionTime:", file_ord_dict[InstanceNumber]['AcquisitionTime']
#    print "  GradientNumber :", file_ord_dict[InstanceNumber]['GradientNumber']
#    print "  BValue         :", file_ord_dict[InstanceNumber]['BValue']

volume_ord_dict = OrderedDict()
volume_number = 0 # Assign numbers to volumes based.
for InstanceNumber, file in file_ord_dict.iteritems():
    AcquisitionTime = file['AcquisitionTime']
    GradientNumber = file['GradientNumber']
    BValue = file['BValue']

    if (not AcquisitionTime in volume_ord_dict):
        volume_ord_dict[AcquisitionTime] = {}
        volume_ord_dict[AcquisitionTime]['files_in_volume'] = 1 # count number of files with same AcquisitionTime
        volume_ord_dict[AcquisitionTime]['GradientNumber'] = GradientNumber
        volume_ord_dict[AcquisitionTime]['BValue'] = BValue
        volume_ord_dict[AcquisitionTime]['volume_number'] = volume_number
        volume_number += 1
    else:
        volume_ord_dict[AcquisitionTime]['files_in_volume'] += 1
        if (GradientNumber != volume_ord_dict[AcquisitionTime]['GradientNumber']):
            raise ValueError("Error: Multiple GradientNumbers for AquisitionTime: "+str(AcquisitionTime))
        if (BValue != volume_ord_dict[AcquisitionTime]['BValue']):
            raise ValueError("Error: Multiple BValues for AquisitionTime: "+str(AcquisitionTime))

files_in_volume_1 = volume_ord_dict[volume_ord_dict.keys()[0]]['files_in_volume']
for AcquisitionTime, volume in volume_ord_dict.iteritems():
    if (volume['files_in_volume'] != files_in_volume_1):
        raise ValueError("Error: number of files per volume varies between volumes.")
    

for AcquisitionTime, volume in volume_ord_dict.iteritems():
    print "volume_number:", volume['volume_number']
#    print "  AcquisitionTime:", AcquisitionTime
    print "  GradientNumber :", volume['GradientNumber']
    print "  BValue         :", volume['BValue']

print "Number of volumes found:",len(volume_ord_dict)

#!/home/sufkes/Envs/venv/bin/python
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
parser.add_argument("in_dir", action='store', type=str, help="path to DICOM files for DTI series to calclate bval and bvec files for")

# Define optional arguments.
parser.add_argument('-o', '--out_dir', action='store', type=str, default=os.getcwd(), help="output directory (Default: current working directory)")
parser.add_argument('-a', '--alt', action='store_true', help='use alternative gradient sequence which was reportedly used for the first "15 or so" scans.')
parser.add_argument('-s', '--slices_per_volume', action='store', type=int, default=None, help="if the AcqusitionTime varies within volumes, you must specify the number of slices per volume.") # not sure how to handle these ones yet; they are really messed in terms of AcquisitionTime; not clear what a volume even is here.
#parser.add_arguemtn('-f', '--flip_yz', action='store_true', help='use original signs on y and z components of gradient vectors. By default, the y and z components of the *original* gradient vectors are flipped')
parser.add_argument('-g', '--guess_slices_per_volume', action='store_true', help='for early scans (those taken on or before 2006-10-12), use this option to automatically guess the number of slices per volume (based on the number of dicom files, and the assumption of exactly 25 volumes); the DICOM files will be assigned to volumes in order based on their InstanceNumber')

# Print help if no args input.
if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit()

# Parse input arguments and store them.
args = parser.parse_args()

# Check for conflicting arguments.
if (args.slices_per_volume != None) and (args.guess_slices_per_volume):
    raise ValueError('Cannot use options -s --slices_per_volume and -g --guess_slices_per_volume together.')
    sys.exit()

# Generate list of all files in directory
file_list = []
for dirname, dirnames, filenames in os.walk(args.in_dir):
    for filename in filenames:
        file_list.append(os.path.abspath(os.path.join(dirname, filename)))

num_files = len(file_list)
print "Creating map for "+str(num_files)+" files."

# Get number of slices per volume if -s or -g options used. 
slices_per_volume = None
if (args.guess_slices_per_volume):
    if (num_files % 25 != 0):
        raise ValueError('Number of files is not divisible by 25. Perhaps the number of volumes is not 25, or DICOM files are missing.')
    else:
        slices_per_volume = num_files/25
elif (args.slices_per_volume != None):
    slices_per_volume = args.slices_per_volume
if (slices_per_volume != None):
    print "Assuming 25 volumes and", slices_per_volume, "slices per volume."

# Create dict mapping instance numbers to gradient directions, b values, and acqusition times.
file_dict = {}
InstanceNumber_max = 0
zero_based_gradient_index = False # whether first non-zero gradient is numbered '#0' or '#1'
ImageOrientationPatient = None
for file_path in file_list:
    ds = pydicom.read_file(file_path, force='force')
    InstanceNumber = int(ds[0x0020, 0x0013].value)
    if (InstanceNumber > InstanceNumber_max):
        InstanceNumber_max = InstanceNumber
    AcquisitionTime = ds[0x0008, 0x0032].value
    SequenceName = ds[0x0018, 0x0024].value
    
    # Get the ImageOrientationPatient tag from one of the DICOM files (they should all be the same). Use this tag to rotate the bvecs.
    if (ImageOrientationPatient == None):
        ImageOrientationPatient = ds[0x0020, 0x0037].value
        V1=np.array([float(ImageOrientationPatient[0]),float(ImageOrientationPatient[1]),float(ImageOrientationPatient[2])]) # FSL unit vector 1
        V2=-np.array([float(ImageOrientationPatient[3]),float(ImageOrientationPatient[4]),float(ImageOrientationPatient[5])]) # FSL unit vector 2
        V3=-np.cross(V1,V2) # FSL unit vector 3 -- Yes, FSL uses a left-handed coordinate system :(
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
        

# Convert dict to OrderedDict sorted based on instance number.
file_ord_dict = OrderedDict()
for InstanceNumber in range(1, InstanceNumber_max+1):
    try:
        file_ord_dict[InstanceNumber] = file_dict[InstanceNumber]
    except KeyError:
        print "Warning: InstanceNumber "+str(InstanceNumber)+" not found."

# Increment b>0 gradient indices if zero-based gradient direction indexing is used.
if zero_based_gradient_index:
    print "Warning: First b>0 volume gradient direction is labelled as '#0'. Incrementing GradientNumber by 1 for all b>0 volumes. Also reassigning AcquisitionTime for each file to the Acquition time of the first volume (really sketchy)."
    for InstanceNumber, file in file_ord_dict.iteritems():
        if ("#" in file['SequenceName']):
            file['GradientNumber'] += 1
if (slices_per_volume != None): # Change all AcquistionTime values in a volume to the AcquisitionTime of the first slice in the volume. This is sketchy; do it for sorting purposes; we want to process scans that have intra-volume AcquistionTime variation the same as those without.
    fake_AcquisitionTime = 1.0
    print "Warning: Modifying AcquisitionTime such that the same AcquisitionTime is used within each volume."
    for InstanceNumber, file in file_ord_dict.iteritems():
        if (InstanceNumber % slices_per_volume == 1): # if "first" slice in volume (slice with lowest IntanceNumber)
            fake_AcquisitionTime += 1.0
            file_ord_dict[InstanceNumber]['AcquisitionTime'] = fake_AcquisitionTime
            print "Treating slice with InstanceNumber", InstanceNumber, "as first slice in a volume. Assuming", slices_per_volume,"slices per volume."
        else:
            file_ord_dict[InstanceNumber]['AcquisitionTime'] = fake_AcquisitionTime

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
    print "  AcquisitionTime:", AcquisitionTime
    print "  GradientNumber :", volume['GradientNumber']
    print "  BValue         :", volume['BValue']

num_volumes = len(volume_ord_dict)
print "Number of volumes found:", num_volumes

## Generate bval and bvec file using the standard list of directions, and the ordering that has just been determined.
if (not args.alt): # The gradient numbers used in the 1st and 2nd list of 'siemens_12_directions.doc' correspond to the row numbers in the following array, except that the y and z components have been sign flipped.
    grad_dirs = np.array([[0.0, 0.0, 0.0],
                          [0.862836817, -0.357430152, 0.357430152],
                          [0.862836817, 0.357430152, 0.357430152],
                          [0.862836817, 0.357430152, -0.357430152],
                          [0.862836817, -0.357430152, -0.357430152],
                          [0.357430152, -0.357430152, -0.862836817],
                          [0.357430152, -0.862836817, -0.357430152],
                          [0.357430152, -0.862836817, 0.357430152],
                          [0.357430152, -0.357430152, 0.862836817],
                          [0.357430152, 0.357430152, 0.862836817],
                          [0.357430152, 0.862836817, 0.357430152],
                          [0.357430152, 0.862836817, -0.357430152],
                          [0.357430152, 0.357430152, -0.862836817]
                          ])

# The gradient numbers used in 'Directions for P1- 15ish.txt' correspond to the row number in the following array, except that the y and z components have been sign flipped.
else: # if using alternative gradient direction series.
    grad_dirs = np.array([[0, 0, 0],
                          [0.894427191, 0, -0.447213595],
                          [0, -0.447213595, -0.894427191],
                          [0.447213595, -0.894427191, 0],
                          [0.894427191, -0.447213595, 0],
                          [0, -0.894427191, -0.447213595],
                          [0.447213595, 0, -0.894427191],
                          [0.894427191, 0, 0.447213595],
                          [0, 0.447213595, -0.894427191],
                          [-0.447213595, -0.894427191, 0],
                          [0.894427191, 0.447213595, 0],
                          [0, -0.894427191, 0.447213595],
                          [-0.447213595, 0, -0.894427191]
                          ])

bval = np.zeros(shape=(num_volumes,1))
bvec = np.zeros(shape=(num_volumes,3))

for AcquisitionTime, volume in volume_ord_dict.iteritems():
    volume_number = volume['volume_number']
    GradientNumber = volume['GradientNumber']
    BValue = volume['BValue']

    bval[volume_number, 0] = BValue

    # Choose appropriate gradient vector based on number assigned in the SequenceName, and put it in the FSL coordinate system.
    grad_dir = grad_dirs[GradientNumber, :]
    bvec[volume_number, 0] = np.dot(V1, grad_dir)
    bvec[volume_number, 1] = np.dot(V2, grad_dir)
    bvec[volume_number, 2] = np.dot(V3, grad_dir)

# Transpose the gradient directions to match FSL bvec file format
bvec = bvec.T
bval = bval.T

# Save bval and bvec files.
bvec_out_path = os.path.join(args.out_dir, 'guess.bvec')
if (os.path.exists(bvec_out_path)) and (not args.force):
    raise ValueError('Output path exists.')
else:
    with open(bvec_out_path, 'w') as fh:
        np.savetxt(fh, bvec, fmt='%.17g')
bval_out_path = os.path.join(args.out_dir, 'guess.bval')
if (os.path.exists(bval_out_path)) and (not args.force):
    raise ValueError('Output path exists.')
else:
    with open(bval_out_path, 'w') as fh:
        np.savetxt(fh, bval, fmt='%.17g')

# Inspect AcquisitionTime values for early scans in which AcquisitionTime varies within volumes. Can't seem to find a logical dependency between AcquisitionTime and InstanceNumber.
grad_num_to_time = {}
for AcquisitionTime, volume in volume_ord_dict.iteritems():
    volume_number = volume['volume_number']
    GradientNumber = volume['GradientNumber']
    BValue = volume['BValue']

    if (not GradientNumber in grad_num_to_time):
        grad_num_to_time[GradientNumber] = {}
        grad_num_to_time[GradientNumber]['AcquisitionTime_list'] = [float(AcquisitionTime)]
    else:
        grad_num_to_time[GradientNumber]['AcquisitionTime_list'].append(float(AcquisitionTime))

for GradientNumber, value in grad_num_to_time.iteritems():
    time_min = min(value['AcquisitionTime_list'])
    time_max = max(value['AcquisitionTime_list'])
    grad_num_to_time[GradientNumber]['time_min'] = time_min
    grad_num_to_time[GradientNumber]['time_max'] = time_max
    print 'GradientNumber:', GradientNumber
    print '  AcquisitionTime_list:', value['AcquisitionTime_list']
#    print '  Min AcquisitionTime:', time_min
#    print '  Max AcquisitionTime:', time_max


        
    

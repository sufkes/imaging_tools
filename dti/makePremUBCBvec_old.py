#!/usr/bin/env python
import os, sys
import pydicom
import numpy as np
import argparse
from copy import deepcopy

# Create argument parser.
description = '''
Generate FSL bvec files for UBC Preterm DTI series.
'''
parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

# Define positional arguments. 
parser.add_argument("in_file", help="path to any DICOM file in the series for which the bvec is being computed")

# Define optional arguments.
parser.add_argument('-o', '--out_dir', action='store', type=str, default=os.getcwd(), help="output directory (Default: current working directory)")
parser.add_argument('-n', '--num_reps', action='store', type=int, default=3, help="number of times the sequence of directions was repeated (Default 3).")
parser.add_argument('-r', '--reverse', action='store_true', help='assume sequence of directions was reversed.')
parser.add_argument('-f', '--force', action='store_true', help='overwrite bval and bvec files if the path exists.')
parser.add_argument('-g', '--group_gradients', action='store_true', help='assume volumes are ordered such that all volumes corresponding to a particular b value and gradient direction are consecutive (i.e. grad1, grad1, grad2, grad2, grad3, grad3; instead of the grad1, grad2, grad3, grad1, grad2, grad3). This arrangement appears to be uncommon, but happens a few times.')
parser.add_argument('-a', '--alt', action='store_true', help='use alternative gradient sequence which was reportedly used for the first "15 or so" scans.')
parser.add_argument('-s', '--shift_zero', action='store_true', help='move the b=0 volume from the start of the standard series to the end (this option does not reverse the entire series)')

# Print help if no args input.
if (len(sys.argv) == 1):
    parser.print_help()
    sys.exit()

# Parse input arguments and store them.
args = parser.parse_args()

# Open sample dicom in series (to get the ImageOrientationPatient tag for the series).
if (not os.path.exists(args.in_file)):
    print "Input file does not exist."
    sys.exit()
else:
    ds=pydicom.read_file(args.in_file, force='force')
img_plane_position=ds[0x0020, 0x0037].value

# Get FSL-frame unit vectors?
V1=np.array([float(img_plane_position[0]),float(img_plane_position[1]),float(img_plane_position[2])])
V2=-np.array([float(img_plane_position[3]),float(img_plane_position[4]),float(img_plane_position[5])])
V3=-np.cross(V1,V2) # Yes, it is a left-handed coordinate system.

# The following sequence of gradient directions was guessed using the following method:
# (1) Look at the directions in the siemens_12_directions.doc file that I found buried in /hpf/largeprojects/smiller/images/prem_ubc/admin
# (2) Compare the sequence of directions with those listed in the few UBC preterm DTI series that have the gradient directions included.
# (3) Normalize the directions in the siemens_12_directions.doc, and flip the Y and Z coordinate signs simply to match the directions in the header files.
# At this point, the guessed gradient directions (in array below) match the gradient directions in the header file approximately (maybe to 1 or 2 decimal places).
# Note: In Justin Foong's script, he said something about flipping the X direction due to an inconsistency between FSL and GE's definition of the X direction (although, it appears that all the UBC Preterm scans were done on Siemens scanners). I flip both the Y and Z direction in order to match the bvec files output by dcm2niix. It could be that Justin's bvecs were the negatives of my bvecs, but it might not matter if the diffusion tensor is symmetric (I'm not sure about this).

# Load gradient directions from text file (old method).
#with open('/Users/steven ufkes/scripts/image_processing/guess_gradient_dirs_from_doc_after_norm_and_sign_flip_matches_20121212_BC0268_V02-009.txt','r') as fh:
#    grad_dirs_list = fh.readlines()
#grad_dirs = np.array(grad_dirs_list[0].split()).reshape(1,3)
#for line in grad_dirs_list[1:]:
#    grad_dir = np.array(line.split()).reshape(1,3)
#    grad_dirs = np.append(grad_dirs, grad_dir, axis=0)
#grad_dirs = grad_dirs.astype(np.float)

# Write the directions directly into a NumPy array for portability.
if (not args.alt):
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
                          [0.357430152, 0.357430152, -0.862836817],
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

# following alternative sequence also had y and z components flipped.
else: # if using alternative gradient direction series 
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
                          [-0.447213595, 0, -0.894427191],
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

bval = np.array([[0],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [600],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700],
                 [700]
                 ])

bvec = np.zeros(shape=grad_dirs.shape)

# Get the gradient directions in the FSL frame (hopefully).
for row in range(len(grad_dirs[:,0])):
    grad_dir = grad_dirs[row,:] # gradient direction in scanner bore frame?
    bvec[row,0] = np.dot(V1, grad_dir)
    bvec[row,1] = np.dot(V2, grad_dir)
    bvec[row,2] = np.dot(V3, grad_dir)

# Shift the b=0 volume from the start of the series to the end (without reversing the whole sequence).
if (args.shift_zero):
    bval = np.roll(bval, -1, axis=0)
    bvec = np.roll(bvec, -1, axis=0)

# Duplicate gradient direction sequence if the sequence was performed many times (usually 3 times)
if (args.num_reps > 1):
    bvec_org = deepcopy(bvec)
    bval_org = deepcopy(bval)
    if (not args.group_gradients):
        for ii in range(args.num_reps-1):
            bvec = np.append(bvec, bvec_org, axis=0)
            bval = np.append(bval, bval_org, axis=0)
    else: # if volumes with same gradient direction are consecutive
        bvec = np.zeros(shape=(3*bvec_org.shape[0], bvec_org.shape[1]))
        bval = np.zeros(shape=(3*bval_org.shape[0], bval_org.shape[1]))
        for start in range(args.num_reps):
            bvec[start::args.num_reps, :] = bvec_org
            bval[start::args.num_reps, :] = bval_org

# Reverse the sequence if the NIFTI file is suspected to have the volumes in reverse order (e.g. B=0 is the 25th image).
if (args.reverse):
    bvec = bvec[::-1,:]
    bval = bval[::-1,:]
    
# Transpose the gradient directions to match FSL bvec file format
grad_dirs = grad_dirs.T
bvec = bvec.T
bval = bval.T

#print grad_dirs
#print bvec

bvec_out_path = os.path.join(args.out_dir, 'makePremUBCBvec.bvec')
if (os.path.exists(bvec_out_path)) and (not args.force):
    raise ValueError('Output path exists.')
else:
    with open(bvec_out_path, 'w') as fh:
        np.savetxt(fh, bvec, fmt='%.17g')
bval_out_path = os.path.join(args.out_dir, 'makePremUBCBvec.bval')
if (os.path.exists(bval_out_path)) and (not args.force):
    raise ValueError('Output path exists.')
else:
    with open(bval_out_path, 'w') as fh:
        np.savetxt(fh, bval, fmt='%.17g')

#!/usr/bin/env python

import pydicom
from struct import unpack
import numpy as np

import sys, os
import glob
from natsort import natsorted

"""Get bvec and bval file from DTI series in DICOM format.

Usage:
./getBvecBval.py <input dir> <scanner type> [step size] [num files]

Input the character corresponding to the scanner manufacturer:
  Siemens: s
  Philips: p
"""

input_dir = str(sys.argv[1])
man = str(sys.argv[2]) # Philips = 'p', Siemens = 's'.

# Specify step size
if (len(sys.argv) >= 4):
    step_size = int(sys.argv[3]) # specify step size for files (i.e. read every [step_size]th file)
else:
    step_size = 1

# Specify the number of different B fields
if (len(sys.argv) == 5):
    num_files = int(sys.argv[4])
else:
    num_files = None

# Find files in input directory:
file_list = glob.glob(input_dir+"/*dcm")
file_list = natsorted(file_list)

if (num_files != None):
    file_list = file_list[:num_files]

print "Looking at",len(file_list),"files with step size",step_size

first = True
for input_file in file_list[::step_size]:
    
    ds=pydicom.read_file(input_file, force='force')
#    print ds
    
    if man == 's':
        bval = ds[0x0019,0x100c].value 
    elif man == 'p':
        bval = ds[0x0018,0x9087].value
#    print "bval (orig):", bval
        
    try:
        if man == 's':
            bvec = ds[0x0019,0x100e].value
        elif man == 'p':
            bvec = ds[0x0018,0x9089].value
    except KeyError:
        bvec = np.zeros(3)
#    print "bvec (orig):", bvec
            
    if (man == 's') or (man == 'p'):
        img_plane_position=ds[0x0020, 0x0037].value
#    print "img_plane_position:", img_plane_position
                
    V1=np.array([float(img_plane_position[0]),float(img_plane_position[1]),float(img_plane_position[2])])
    V2=np.array([float(img_plane_position[3]),float(img_plane_position[4]),float(img_plane_position[5])])
    V3=np.cross(V1,V2)
                
    pbvec=np.zeros((3,1))
    pbvec[0][0] = np.dot(V1,bvec)
    pbvec[1][0] = -np.dot(V2,bvec) # 2019-05-15 Steven Ufkes: This sign flip is probably here to comply with FSL's absurd left-handed cooredinate system.
    pbvec[2][0] = np.dot(V3,bvec) 
                
#    print "bvec:"
#    print pbvec
    if first:
        bvec_file = np.zeros((3,1))
        bvec_file[:,0] = pbvec[:,0]
        first = False
        #print bvec_file

        bval_file = np.zeros((1,1))
        bval_file[0][0] = bval
    else:
        bvec_file = np.append(bvec_file, pbvec, axis=1)
        bval_file = np.append(bval_file, [[bval]], axis=1)
print bvec_file
print bval_file

with open(os.path.join(input_dir,"output.bvec"), 'wb') as fh:
    np.savetxt(fh, bvec_file, fmt='%.17g')

with open(os.path.join(input_dir,"output.bval"), 'wb') as fh:
    np.savetxt(fh, bval_file, fmt='%.17g')


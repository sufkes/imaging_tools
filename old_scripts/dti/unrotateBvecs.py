#!/usr/bin/env python
import sys
import pydicom
import numpy as np

# Unrotate gradient directions for 20130116_BC0269_V02/008-DIFF-12-6007003-AVE-UID_1-3-12-2-1107-5-99-2-1561-30000012120319052546800006676

x_hat = np.array([1,0,0])
y_hat = np.array([0,1,0])
z_hat = np.array([0,0,1])

x_new = np.array([0.98340339391276, 0.09564653800052, -0.1541736183962])
y_new = np.array([-0.1207467531479, 0.97926424502878, -0.1626707103705])
z_new = np.cross(x_new, y_new)

#print np.sqrt(x_new[0]**2+x_new[1]**2+x_new[2]**2)
#print np.sqrt(y_new[0]**2+y_new[1]**2+y_new[2]**2)
#print np.sqrt(z_new[0]**2+z_new[1]**2+z_new[2]**2)

# Define rotation matrix
#r1 = np.zeros(shape=(3,3))
#r1 = np.zeros(shape=(3,3))
r1 = np.array([[np.dot(x_hat, x_new), np.dot(x_hat, y_new), np.dot(x_hat, z_new)],
          [np.dot(y_hat, x_new), np.dot(y_hat, y_new), np.dot(y_hat, z_new)],
          [np.dot(z_hat, x_new), np.dot(z_hat, y_new), np.dot(z_hat, z_new)]
          ])
#print r1

r2 = r1.T
#print r2

#r1_inv = np.linalg.inv(r1)# just the transpose.
#print r1_inv

#r2_inv = np.linalg.inv(r2)# just the transpose.
#print r2_inv

#print np.matmul(r1,r2)
#print np.matmul(r1,r1_inv)
#print np.matmul(r2,r1)
#print np.matmul(r2,r2_inv)

grad_dir = pydicom.read_file(str(sys.argv[1]), force='force')[0x0019, 0x100e].value
print grad_dir
print np.matmul(r1,grad_dir)
print np.matmul(r2,grad_dir)


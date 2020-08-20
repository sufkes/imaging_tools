#!/usr/bin/env python
import numpy as np
from scipy.spatial.distance import pdist, squareform
# Distance covariance is defined as in 2014 - Szekely et al.

def euclideanDistance(X):
    a = np.zeros((X.shape[0], X.shape[0]))
    for ii in range(X.shape[0]):
        for jj in range(X.shape[0]):
            for kk in range(X.shape[1]):
                a[ii, jj] += np.power(X[ii,kk] - X[jj,kk], 2)
            a[ii, jj] = np.sqrt(a[ii,jj])
    return a

def uCenteredDistance(X): # same result as R energy
    a = euclideanDistance(X)
    
    A = np.zeros(a.shape)
    n = a.shape[0]
    for ii in range(n):
        for jj in range(n):
            if ii == jj:
                continue # A[ii,jj] already set to zero
            else:
                A[ii,jj] = a[ii,jj] - 1/float(n-2)*np.sum(a[ii,:]) - 1/float(n-2)*np.sum(a[:,jj]) + 1/float(n-1)/float(n-2)*np.sum(a)
    return A

def doubleCenteredDistance(X): # same result as R energy
    a = euclideanDistance(X)
    
    n = a.shape[0]
    A = np.zeros(a.shape)
    for ii in range(n):
        for jj in range(n):
            A[ii, jj] = a[ii, jj] - np.mean(a[ii, :]) - np.mean(a[:, jj]) + np.mean(a)
    return A

def center(X, double=False):
    if double:
        return doubleCenteredDistance(X)
    else:
        return uCenteredDistance(X)

def distanceCovarianceSquared_slow(X, Y, double=False):
    A = center(X, double=double)
    B = center(Y, double=double)

    n = A.shape[0]
    if double:
        K = n**2
    else:
        K = n*(n-3)
    K = float(K)
    val = 0
    for ii in range(n):
        for jj in range(n):
            val += A[ii,jj]*B[ii,jj]
    val = val/K
    return val

def distanceCorrelation_slow(X, Y, double=False, verbose=False):

    d_cov_squared = distanceCovarianceSquared_slow(X, Y, double=double)
    if verbose:
        print "distanceCovariance(X,Y)^2:", d_cov_squared
        print "sqrt(|distanceCovariance(X,Y)|^2):", np.sqrt(np.abs(d_cov_squared))

    d_var_x = np.sqrt(distanceCovarianceSquared_slow(X, X, double=double))
    d_var_y = np.sqrt(distanceCovarianceSquared_slow(Y, Y, double=double))
    if verbose:
        print "distanceCovariance(X):", d_var_x
        print "distanceCovariance(Y):", d_var_y

    # Handle zero variance cases as in 2014-Szekely paper.
    try:
        d_cor_squared = d_cov_squared/(d_var_x*d_var_y)
    except ZeroDivisionError:
        d_cor_squared = 0.0

    geerligs_measure = np.sqrt(max(0.0, d_cor_squared)) # This is the measure used in 2016-L Geerligs, where they state that negative values correspond to independence.
    if verbose:
        print "distanceCorrelation(X, Y)^2:", d_cor_squared
        print "sqrt(|distanceCorrelation(X, Y)^2|):", np.sqrt(np.abs(d_cor_squared))
        print "sqrt(max(0, distanceCorrelation(X, Y)^2)):", geerligs_measure
    return geerligs_measure

# This function for the double-centered distance correlation was taken from https://gist.github.com/satra/aa3d19a12b74e9ab7941
#def distanceCorrelation2_double(X, Y):
#    X = np.atleast_1d(X)
#    Y = np.atleast_1d(Y)
#    if np.prod(X.shape) == len(X):
#        X = X[:, None]
#    if np.prod(Y.shape) == len(Y):
#        Y = Y[:, None]
#    X = np.atleast_2d(X)
#    Y = np.atleast_2d(Y)
#    n = X.shape[0]
#    if Y.shape[0] != X.shape[0]:
#        raise ValueError('Number of samples must match')
#    a = squareform(pdist(X))
#    b = squareform(pdist(Y))
#    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
#    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()
#    
#    dcov2_xy = (A * B).sum()/float(n * n)
#    dcov2_xx = (A * A).sum()/float(n * n)
#    dcov2_yy = (B * B).sum()/float(n * n)
#    dcor = np.sqrt(dcov2_xy)/np.sqrt(np.sqrt(dcov2_xx) * np.sqrt(dcov2_yy))
#    print "dcor:", dcor
#    return

def distanceCorrelation(X, Y, double=False):
    """Compute the squared distance correlation as defined in Szekely et al. (2014). Returns a tuple of four values:
    dcor^2(X,Y), dcov^2(X), dcov^2(Y), dcov^2(X, Y)

    Input matrices must have same lenght of first dimension.
    """

    a = squareform(pdist(X))
    b = squareform(pdist(Y))

    if double:
        def center(a):
            n = float(a.shape[0])
            A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
            #np.fill_diagonal(A, 0.0)
            return A
    else:
        def center(a):
            n = float(a.shape[0])
            A = a - n/(n-2)*a.mean(axis=0)[None, :] - n/(n-2)*a.mean(axis=1)[:, None] + n*n/(n-1)/(n-2)*a.mean()
            np.fill_diagonal(A, 0.0)
            return A
    
    A = center(a)
    B = center(b)

    # Calculate covariance
    if double:
        def dcov2(A,B):
            n = A.shape[0]
            K = float(n**2)
            return (A*B).sum()/K
    else:
        def dcov2(A,B):
            n = A.shape[0]
            K = float(n*(n-3))
            return (A*B).sum()/K

    dcov2_XX = dcov2(A,A)
    dcov2_YY = dcov2(B,B)
    dcov2_XY = dcov2(A,B)
    #print "cov2_XY:", cov2_XY
    
    dcor2 = dcov2_XY/np.sqrt(dcov2_XX*dcov2_YY)
    
    #print "dcor^2:", dcor2
    return dcor2, dcov2_XX, dcov2_YY, dcov2_XY

def normalizeTimeSeries(X):
    return (X-X.mean(axis=0))/X.std(axis=0, ddof=1)
    
#    n = float(X.shape[0])
#    for jj in range(X.shape[1]):
#        col_mean = np.mean(X[:, jj])
#        col_std = np.std(X[:, jj])*np.sqrt(n/(n-1))
#        X[:, jj] = (X[:, jj] - col_mean)/col_std
#    X = X - X.mean(axis=0)[]
#    return X

if (__name__ == "__main__"):
    # X and Y are matrices of n time points by v voxels.
    X = np.array([5.1, 3.5, 4.9, 3.0, 4.7, 3.2, 4.6, 3.1, 5.0, 3.6]).reshape(5,2)
    Y = np.array([7.0, 3.2, 6.4, 3.2, 6.9, 3.1, 5.5, 2.3, 6.5, 2.8]).reshape(5,2)
    #Y=X
    #X = normalizeTimeSeries(X)
    #Y = normalizeTimeSeries(Y)

    if 1:
        print "Double-centered distanceCorrelation(X, Y):"
        print "="*60
        print "SLOW:"
        distanceCorrelation_slow(X, Y, double=True, verbose=True)
        print "-"*40
        print "FAST:"
        print distanceCorrelation(X, Y, double=True)
        print
    if 1:
        print "U-centered distanceCorrelation(X, Y):"
        print "="*60
        print "SLOW:"
        distanceCorrelation_slow(X, Y, verbose=True)
        print "-"*40
        print "FAST:"
        print distanceCorrelation(X, Y)


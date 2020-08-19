#!/usr/bin/env python
import numpy as np

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

def distanceCovarianceSquared(X, Y, double=False):
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

#def dVar_squared(X, double=False): # can just use the distanceCovarianceSquared function with both args set to the same array.
    #A = center(X, double=double)
    
    #n = A.shape[0]
    #if double:
        #K = float(np.power(n,2))
    #else:
        #K = float(n*(n-3))
    #val = 1.0/K*np.sum(np.power(A, 2))
    #return val
    
def distanceCorrelation(X, Y, double=False, verbose=False):

    d_cov_squared = distanceCovarianceSquared(X, Y, double=double)
    if verbose:
        print "distanceCovariance(X,Y)^2:", d_cov_squared
        print "sqrt(|distanceCovariance(X,Y)|^2):", np.sqrt(np.abs(d_cov_squared))

    d_var_x = np.sqrt(distanceCovarianceSquared(X, X, double=double))
    d_var_y = np.sqrt(distanceCovarianceSquared(Y, Y, double=double))
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

def normalizeTimeSeries(X):
    n = float(X.shape[0])
    for jj in range(X.shape[1]):
        col_mean = np.mean(X[:, jj])
        col_std = np.std(X[:, jj])*np.sqrt(n/(n-1))
        X[:, jj] = (X[:, jj] - col_mean)/col_std
    return X

if (__name__ == "__main__"):
    # X and Y are matrices of n time points by v voxels.
    X = np.array([5.1, 3.5, 4.9, 3.0, 4.7, 3.2, 4.6, 3.1, 5.0, 3.6]).reshape(5,2)
    Y = np.array([7.0, 3.2, 6.4, 3.2, 6.9, 3.1, 5.5, 2.3, 6.5, 2.8]).reshape(5,2)
    #X = normalizeTimeSeries(X)
    #Y = normalizeTimeSeries(Y)
    
    print "U-centered distanceCorrelation(X, Y):"
    distanceCorrelation(X, Y, verbose=True)
    print
    print "Double-centered distanceCorrelation(X, Y):"
    distanceCorrelation(X, Y, double=True, verbose=True)

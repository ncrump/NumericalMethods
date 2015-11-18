"""
Created on Sun Nov 15 14:49:36 2015
CSI 690, Assignment 9
Nick Crump
"""

# Problem 27.11/12 Power method to solve eigenvalue equation
# Finds the highest and lowest eigenvalues and their eigenvectors
"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import numpy as np


# input parameters
# -------------------------------------------------
# input matrix
A = np.array([[2, 8, 10],\
              [8, 4, 5], \
              [10,5, 7]])

# input error tolerance
tol = 1e-6
# -------------------------------------------------

# solves eigenvalue equation Ax = Lx

# get matrix dimension
dim   = len(A)
# set initial guess for x
xMax  = np.ones((dim,1),float)
xMin  = np.ones((dim,1),float)
# set initial error
errMx = 1
errMn = 1
# set iteration variables
itrMx = 0
itrMn = 0
tmpMx = 0
tmpMn = 0

# find max eigenvalue (L) & eigenvector (x)
# iterate over matrix A
while errMx > tol:
    # multiply Ax
    xMax = np.dot(A,xMax)
    # get eigenvalue
    mx = np.max(xMax)
    mn = np.min(xMax)
    if abs(mx) >= abs(mn):
        LMax = mx
    else:
        LMax = mn
    # normalize eigenvector
    xMax  = xMax/LMax
    # get relative error
    errMx = abs((LMax-tmpMx)/LMax)
    tmpMx = LMax
    # count iterations
    itrMx += 1

# find min eigenvalue (L) & eigenvector (x)
# this time iterate over inverse of A
Ainv = np.linalg.inv(A)
while errMn > tol:
    # multiply inverse(A)x
    xMin = np.dot(Ainv,xMin)
    # get 1/eigenvalue
    mx = np.max(xMin)
    mn = np.min(xMin)
    if abs(mx) >= abs(mn):
        LMin = mx
    else:
        LMin = mn
    # normalize eigenvector
    xMin  = xMin/LMin
    # get relative error
    errMn = abs((LMin-tmpMn)/LMin)
    tmpMn = LMin
    # count iterations
    itrMn += 1
# get eigenvalue from reciprocal
LMin = 1.0/LMin

# check results Ax-Lx = 0
diffMx = np.round(np.dot(A,xMax) - np.dot(LMax,xMax),3)
diffMn = np.round(np.dot(A,xMin) - np.dot(LMin,xMin),3)

# print max results
print '\n'
print 'max eigenvalue:'
print '%11.6f' % LMax
print 'eigenvector:'
print xMax
print 'iterations:'
print '%4i' % itrMx
# print min results
print '\n'
print 'min eigenvalue:'
print '%11.6f' % LMin
print 'eigenvector:'
print xMin
print 'iterations:'
print '%4i' % itrMn
# print check results
print '\n'
print 'check Ax-Lx=0:'
print 'max:'
print diffMx
print 'min:'
print diffMn
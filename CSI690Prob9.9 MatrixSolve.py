"""
Created on Thu Oct 08 18:46:50 2015
CSI 690, Assignment 4
Nick Crump
"""

# Problem 9.9: Gaussian elimination to solve matrix equation

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import numpy as np


# input parameters
# ---------------------------------------------------------------------
# input matrix as row elements
# note: precision must be float or better to get good results!
a = np.array([[4,1,-1],[5,1,2],[6,1,1]],float)
b = np.array([[-2],[4],[6]],float)

# input tolerance for detecting singularity (near-zero determinant)
tol = 10e-6
# ---------------------------------------------------------------------

# get matrix dimensions
arow = np.size(a,axis=0)
acol = np.size(a,axis=1)
brow = np.size(b,axis=0)
bcol = np.size(b,axis=1)

# check input is valid
if arow != acol or arow != brow or bcol != 1:
    print '\nwarning: matrix is not square. check input.'
    exit()

# initialize arrays
acpy = np.copy(a)
xsol = np.zeros([arow,1])
atmp = np.zeros(arow)
btmp = np.zeros(1)
sign = 0

# forward elimination to reduce matrix
for i in range(0,arow-1):
    # partial pivoting to avoid divide by zero
    indx = i
    pmax = abs(a[i,i])
    for p in range(i+1,acol):
        ai = abs(a[p,i])
        if ai > pmax:
            pmax = ai
            indx = p
    # exchange rows if pivot value not max
    if indx != i:
        atmp[:]   = a[indx,:]
        btmp[:]   = b[indx,:]
        a[indx,:] = a[i,:]
        b[indx,:] = b[i,:]
        a[i,:]    = atmp[:]
        b[i,:]    = btmp[:]
        sign += 1
    # main Gaussian elimination
    # final matrix is upper right triangular
    for j in range(i+1,arow):
        f = a[j,i]/float(a[i,i])
        for k in range(i+1,arow):
            a[j,k] = a[j,k] - f*a[i,k]
        b[j]   = b[j] - f*b[i]
        
# check that matrix is not singular (zero determinant)
# determinant for triangular matrix is product of diagonals
det = np.prod(np.diag(a))*(-1)**sign
if abs(det) < tol:
    print '\nwarning: singular matrix. check input.'
    exit()

# back substitution to solve matrix
mx = arow-1
xsol[mx] = b[mx]/float(a[mx,mx])
for i in range(arow-1,-1,-1):
    s = b[i]
    for j in range(i+1,arow):
        s = s - a[i,j]*xsol[j]
    xsol[i] = s/float(a[i,i])
    
# check solution by substituting back into original equations
xchk = np.sum(acpy*np.transpose(xsol),axis=1).reshape(arow,1)

# print solution
print '\nsolution:'
print xsol
print 'determinant:'
print det
print 'solution check:'
print xchk
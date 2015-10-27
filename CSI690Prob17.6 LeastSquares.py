"""
Created on Thu Oct 22 21:35:13 2015
CSI 690, Assignment 6
Nick Crump
"""

# Problem 17.6 Use least squares polynomial regression to fit a data set

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


# input parameters
# -------------------------------------------------
# input data set
x = np.array([1,2,3,4,5,6,7,8,9],float)
y = np.array([1,1.5,2,3,4,5,8,10,13],float)

# input degree of polynomial fit
n = 2
# -------------------------------------------------

# fit model is polynomial: y = c0 + c1(x) + c2(x**2) +...+ cn(x**n)

pts     = len(x)  # number of poins in data set
polyval = []      # array to store poly fit values

# if degree of polynomial is not less than number of points print error
if n+1 > pts:
    print '\nWarning: Polynomial degree must be less than number of data points'
    print 'Check input'
    sys.exit()

# compute fit function using polynomial least squares
else:
    matrixA = np.zeros((pts,n+1))             # initialize matrix A
    matrixB = np.zeros((pts,1))               # initialize matrix b

    # loop to populate arrays
    for i in range(pts):
        matrixB[i][0] = y[i]                  # b gets y values
        for j in range(n+1):
            matrixA[i][j] = (x[i])**j         # A gets x values^j

    # create normal equations (A^T)Ax=(A^T)b
    At = np.transpose(matrixA)                # A transpose (A^T)
    AtA = np.dot(At,matrixA)                  # matrix multiply (A^T)A
    AtB = np.dot(At,matrixB)                  # matrix multiply (A^T)b
    coeff = np.linalg.solve(AtA,AtB)          # solve normal equation (A^T)Ax=(A^T)b

    # get y-fit values
    for i in range(pts):
        yval = 0
        for j in range(n+1):
            yval = yval + coeff[j,0]*(x[i])**j  # gets new y value using fit coeffs
        polyval.append(yval)                    # collect new y values

    # get standard deviation of fit as "goodness of fit"
    error  = y-polyval
    sumval = np.sum(error**2)
    stddev = (sumval/(pts-(n+1)))**0.5

    # get correlation coefficient of fit as alternate "goodness of fit"
    ymean  = np.sum(y)/pts
    sumfit = np.sum((y-ymean)**2)
    rvalue = ((sumfit-sumval)/sumfit)**0.5

    # plot data and fit
    plt.figure()
    plt.plot(x,y,'bo-',label='data set')
    plt.plot(x,polyval,'r-',label='polynomial fit')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc=2)

    # print results
    print '\npolynomial coefficients'
    for i in coeff.flatten():
        print '%10.6f' % i
    print '\nstandard error'
    print '%10.6f' % stddev
    print '\ncorrelation coefficient'
    print '%10.6f' % rvalue
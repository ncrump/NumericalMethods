"""
Created on Thu Jun 20 06:14:19 2013
PHYS 510, Assignment 5 (This part just for fun)
"""

# Exponential Fitting by Linearization
"""
The program below uses the Least Squares method to fit an exponential to a data set using the method
of linearization. The exponential model is first transormed into a log model in order to linearize
the data and compute the best coefficients which can then be placed back into the exponential model.
This is a common technique of linearizing a nonlinear problem so it can be solved easier. 
NOTE: THIS ONLY WORKS FOR DATA SETS WITH POSITIVE Y VALUES SINCE IT INVOLVES TAKING A LOG(Y)!
"""


import numpy as np
import matplotlib.pyplot as plt
import math

# this part just for fun
# Exponential fitting method for nonlinear data
#*******************************************************************
def ExponentialFit(x,y):
    # ** NOTE: this model only works on data sets for POSITIVE y values since uses log(y) **
    # fit model is exponential: y = c0exp(c1*x), but we have to linearize it to use least squares
    # instead we change the fit model to: ln(y) = ln(c0) + c1(x) = k0 + c1(x)
    # where k0 = ln(c0) to make linear
    # solve for these coeffients k0 and c1 then plug back into exponential model
        
    pts = len(x)                              # number of points in data set
    expval = []                               # array to store exponential fit values
    
    matrixA = np.zeros((pts,2))               # initialize matrix A
    matrixB = np.zeros((pts,1))               # initialize matrix b
    
    # loop to populate arrays
    for i in range(pts):
        matrixB[i][0] = math.log(y[i])        # b gets ln(y) values to linearize fit model
        matrixA[i][0] = 1                     # 1st column of A for k0 coeffs are 1's
        matrixA[i][1] = x[i]                  # 2nd column of A for c1 coeffs are x values
    
    # create normal equation (A^T)Ax=(A^T)b 
    At = np.transpose(matrixA)                # A transpose (A^T)
    AtA = np.dot(At,matrixA)                  # matrix multiply (A^T)A
    AtB = np.dot(At,matrixB)                  # matrix multiply (A^T)b
    
    coeff = np.linalg.solve(AtA,AtB)          # solve normal equation (A^T)Ax=(A^T)b for x
    k0 = coeff[0]                             # name nonlinear coefficient k0 for use below
    c0 = math.exp(k0)                         # name linear coefficient c0 for use below
    c1 = coeff[1]                             # name nonlinear coefficient c1 for use below
    
  
    # we now have the coefficients for the exponential fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        yval = c0*math.exp(c1*x[i])
        expval.append(yval)
                                
        # compute variance from erorr (error=actualY-fitY)
        # NOTE: variance here is calculated on the LINEARIZED equation
        error = math.log(y[i]) - k0 - c1*x[i]
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = round(varsum/(pts-2-1), 4)
    
        
    # plot data set with exponential fit and annotate variance
    plt.title('Exponential Fit by Linearization')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(x,expval, 'r', linewidth=2, label='Exponential')
    plt.annotate('Variance='+str(variance),fontsize=14,xy=(0.15,0.83),xycoords='figure fraction')
    plt.legend()  

    
#ExponentialFit(x,y)
#-------------------------------------------------------------------
# x = x data set
# y = y data set  (y DATA MUST BE POSITIVE)
#------------------------------------------------------------------- 
#*******************************************************************
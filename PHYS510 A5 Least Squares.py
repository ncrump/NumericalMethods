"""
Created on Thu Jun 20 06:14:19 2013
PHYS 510, Assignment 5
"""

# Least Squares Fitting
"""
1. Using the file 'data.txt', find the least squares approximation for polynomials 
   of degree 1-7 and compare these techniques by plotting the variance vs. the degree 
   of the polynomial. Determine which technique gives the best result.
2. Repeat the first problem for an exponential fit. How does the approximation 
   in problem 1 compare to the result just found?
"""


import numpy as np
import matplotlib.pyplot as plt
import math

# import sample data set "data.txt"
#x, y = np.loadtxt('data.txt', skiprows=1, unpack=True)


# Part 1
# Least Squares fitting method for linear data
#*******************************************************************
def LeastSquaresFit(x,y,n):
    # fit model is polynomial: y = c0 + c1(x) + c2(x**2) +...+ cn(x**n)
    
    pts = len(x)  # number of poins in data set
    polyval = []  # array to store poly fit values
    
    # if degree of polynomial is greater than number of points print error
    if n+1 >= pts:
        print 'INVALID ENTRY'
        print 'POLYNOMIAL DEGREE CANNOT BE GREATER THAN NUMBER OF DATA POINTS'
    
    # otherwise compute coefficients for polynomial fit
    # solves for coefficients using matrix solve of 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = y values, A = x values^power
    # this needed because system of equations is overdetermined (#equations > #unknowns)
    else:
        matrixA = np.zeros((pts,n+1))             # initialize matrix A
        matrixB = np.zeros((pts,1))               # initialize matrix b
        
        # loop to populate arrays
        for i in range(pts):
            matrixB[i][0] = y[i]                  # b gets y values
            
            for j in range(n+1):
                matrixA[i][j] = (x[i])**j         # A gets x values^j
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)                # A transpose (A^T)
        AtA = np.dot(At,matrixA)                  # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)                  # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)          # solve normal equation (A^T)Ax=(A^T)b for x
        
        
        # we now have the coefficients for the polynomial fit
        # next plug x values back into fit model and compute new y fit values
        varsum = 0     
        # loop to get y fit values and compute variance on fit                           
        for i in range(pts):
            yval = 0
            for j in range(n+1):
                yval = yval + coeff[j]*(x[i])**j  # gets new y value for each x using fit coeffs
                            
            polyval.append(yval)                  # collect new y values
            
            # compute variance from erorr (error=actualY-fitY)
            error = y[i]-yval
            varsum = varsum + error**2
        
        # using variance as measure for 'goodness of fit'
        # computed as variance = (sum of squares of errors) / (nDataPts-PolyDegree-1)
        variance = round(varsum/(pts-n-1), 4)
        
            
        # plot data set with poly fit and annotate variance
        plt.title('Least Squares Polynomial Fit')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.plot(x,y,'bo', label='Data')
        plt.plot(x,polyval, 'r', linewidth=2, label='Order=' + str(n))
        plt.annotate('Variance='+str(variance),fontsize=14,xy=(0.15,0.83),xycoords='figure fraction')
        plt.legend()   
        
#LeastSquaresFit(x,y,4)
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# n = degree of polynomial fit
#------------------------------------------------------------------- 
#*******************************************************************


# Part 2
# Gauss-Newton Exponential fitting method for nonlinear data
#*******************************************************************
def GaussNewtonFit(x,y,c0,c1,n):
    # fit model is exponential: y = c0exp(c1*x)
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = errors to minimize (actualY-fitY), A = numeric derivative of b
    # the initial guess of c0 and c1 is used to construct both b and A
        
    pts = len(x)                                       # number of points in data set
    expval = []                                        # array to store exponential fit values
    
    # loop over specified number of iterations to update coefficients c0 and c1
    for k in range(n):
        matrixA = np.zeros((pts,2))                    # initialize matrix A
        matrixB = np.zeros((pts,1))                    # initialize matrix b
   
       # loop to populate arrays
        for i in range(pts):
            matrixB[i][0] = y[i]-c0*math.exp(c1*x[i])  # b gets residual errors y[i]-c0exp(c1*x[i])
            matrixA[i][0] = math.exp(c1*x[i])          # 1st column of A gets derivative of b wrt c0
            matrixA[i][1] = c0*x[i]*math.exp(c1*x[i])  # 2nd column of A gets derivative of b wrt c1
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)                     # A transpose (A^T)
        AtA = np.dot(At,matrixA)                       # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)                       # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)               # solve normal equation (A^T)Ax=(A^T)b for x
        c0 = c0 + coeff[0]                             # update c0 for next iteration
        c1 = c1 + coeff[1]                             # update c1 for next iteration
        
    
    # we now have the coefficients for the exponential fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        yval = c0*math.exp(c1*x[i])
        expval.append(yval)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - c0*math.exp(c1*x[i])
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = round(varsum/(pts-2-1), 4)
    
    print c0,c1    
        
    # plot data set with exponential fit and annotate variance
    plt.title('Exponential Fit by Gauss-Newton Method')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(x,expval, 'r', linewidth=2, label='Exponential')
    plt.annotate('Variance='+str(variance),fontsize=14,xy=(0.15,0.83),xycoords='figure fraction')
    plt.legend()
    
    
#GaussNewtonFit(x,y,3.5,-1.5,5)
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# c0 = initial first guess at coefficient c0 in y = c0exp(c1*x)
# c1 = initial first guess at coefficient c1 in y = c0exp(c1*x)
# n = number of iterations to converge on c0 and c1
#------------------------------------------------------------------- 
#******************************************************************* 
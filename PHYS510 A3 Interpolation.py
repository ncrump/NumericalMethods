"""
Created on Tue Jun 04 17:40:50 2013
PHYS 510, Assignment 3 
"""

# Interpolation
"""
Below is the cross section for a resonant scattering experiment of a neutron from a nucleus
-----------------------------------------------------------------------
 E(MeV) |  0   |  25  |  50  |  75  | 100  | 125  | 150  | 175  | 200 | 
Sig(mb) | 10.6 | 16.0 | 45.0 | 83.5 | 52.8 | 19.9 | 10.8 | 8.25 | 4.7 |
-----------------------------------------------------------------------
You are asked to determine values for the cross section at energies between 
the measured values.  we assume that we don't have a theoretical model for 
this process and we simply want to interpolate the measured data. 

Complete the following tasks:
1. Use the improved Neville's algorithm to fit the entire experimental spectrum with 
   one polynomial. This means using an eight degree polynomial to interpolated the nine 
   data points. Use the interpolated values to plot the cross section in steps of 5 Mev.
2. A more reasonable use of interpolation algorithms is to model LOCAL behavior with a 
   small number of data points. Redo the interpolation of the cross section data with 
   only 3 data points at a time. Plot the resulting piece-wise smooth interpolating 
   polynomials in 5 Mev steps.
3. Use the Cubic Spline Interpolation and repeat step 3 (use natural spline condition 
   for the end points).
4. The actual theoretical model for this scattering process is given by the Breit-Wigner 
   function. Compare these values with the following theoretical values:
   Sig0 = 6.389x10**4, Er = 78, and Gamma = 55 
   Plot the theoretical prediction on top of the two interpolated curves
"""

import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

# data set to interpolate
E = [0, 25, 50, 75, 100, 125, 150, 175, 200]
Sig = [10.6, 16.0, 45.0, 83.5, 52.8, 19.9, 10.8, 8.25, 4.7] 


# Warmup
# Lagrange Interpolation method
#*******************************************************************
def LagrangeInterp(x, y, nPoints, nxIncrement):
    t = sym.Symbol('t')  # define symbolic variable
    n = range(nPoints)   # index for interpolation 
    p = 0
    
    for i in n:
        term = 1
        
        # Lagrange Method
        for j in [k for k in n if k != i]:
            term = term * ((t - x[j]) / (x[i] - x[j]))
            
        p = p + y[i]*term  # build symbolic interpolating polynomial
    
    fInterp = sym.lambdify(t,p)                                  # turn symbolic poly into real poly    
    xx = np.arange(x[0], x[nPoints-1]+nxIncrement, nxIncrement)  # create x interped values
    yy = [fInterp(z) for z in xx]                                # get y interped values
    
    # plot
    plt.figure(1)
    plt.title('Data Interpolation' + '\nInterpolated in 5 Mev Steps')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(xx,yy, label=str(nPoints)+' Pts')
    plt.legend()
    

#LagrangeInterp(E,Sig,9,5)
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------    
#*******************************************************************
    
    
# Part 1
# Neville Interpolation method
#*******************************************************************
def NevilleInterp(x, y, nPoints, nxIncrement):
    xx = np.arange(x[0], x[nPoints-1]+nxIncrement, nxIncrement)  # create x interped values
    yy = []  # for storing y interped values
    
    for xval in xx:
    
        arry = y[0:nPoints]                # get intial y values
        arr = np.zeros((nPoints,nPoints))  # create Neville matrix
        arr[0] = arry                      # write y values to matrix
        
        xpos = 0
        
        for i in range(nPoints - 1): 
            xpos = xpos + 1
             
            # Neville Method & build matrix
            for j in range(nPoints - 1 - i):
                term = ((xval-x[j+xpos])*arr[i][j] - (xval-x[j])*arr[i][j+1]) / (x[j]-x[j+xpos])
                arr[i+1][j] = term         # build Neville matrix
                   
        yy.append(term) 
    
    print arr
        
    # plot
    plt.figure(2)
    plt.title('Neville Polynomial Interpolation'+'\n'+str(nPoints)+' Points Interpolated in 5 Mev Steps')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(xx,yy, label=str(nPoints)+' Pts')
    plt.legend()
    

NevilleInterp(E,Sig,4,5)
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------
#*******************************************************************
    

# Part 2
# Lagrange Piecewise Interpolation method
#*******************************************************************
def LagrangeInterpPWise(x, y, nPoints, nLocal, nxIncrement):
    t = sym.Symbol('t')                       # define symbolic variable
    n = range(0, nPoints-nLocal+1, nLocal-1)  # create x interped values
    
    xx = []  # for storing x interped values
    yy = []  # for storing y interped values
    
    # iterate over data points
    for loc in n:
        p = 0
        pts = range(loc, loc+nLocal)  # define local points to iterate over
    
        # iterate over local points
        for i in pts:
            term = 1
            
            # Lagrange Method
            for j in [k for k in pts if k != i]:
                term = term * ((t - x[j]) / (x[i] - x[j]))
                
            p = p + y[i]*term  # build symbolic interpolating polynomial
        
        fInterp = sym.lambdify(t,p)                                        # turn symbolic poly into real poly
        xx0 = np.arange(x[loc], x[loc+nLocal-1]+nxIncrement, nxIncrement)  # get local interped x values
        yy0 = [fInterp(z) for z in xx0]                                    # get local interped y values
        # write local interped x and y values to array
        # NOTE: duplicate data points exist at locations where piecewise function comes together
        # NOTE: but since they're equal it doesn't matter
        xx = np.concatenate((xx,xx0))
        yy = np.concatenate((yy,yy0))

            
    # plot
    plt.figure(3)
    plt.title('Lagrange PieceWise Interpolation'+'\n'+str(nLocal)+' Point Local Interpolation in 5 Mev Steps')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(xx,yy, label=str(nLocal)+' Pts')
    plt.legend()


#LagrangeInterpPWise(E,Sig,9,3,5)
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nLocal = number of local points to use for piecewise
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------
#*******************************************************************


# Part 3
# Cubic Spline Interpolation method
# NOTE: uses 'natrual spline' condition (set 2nd derivative to zero at endpoints)
#*******************************************************************
def CubicSpline(x, y, nPoints, nxIncrement):
    # start by solving for double primes to use in cubic spline equation
    # create matrices for Ax=b matrix solver to get values of double primes
    dprimesA = np.zeros((nPoints,nPoints))  # create NxN matrix A for double primes
    dprimesB = np.zeros((nPoints))          # create Nx1 matrix b for double primes
    
    # build matrix A and b to solve for double primes
    for i in range(1, nPoints-1):
        # build matrix A of coefficients
        dprimesA[i][i-1] = x[i]-x[i-1]
        dprimesA[i][ i ] = 2*(x[i]-x[i-1]) + 2*(x[i+1]-x[i])
        dprimesA[i][i+1] = x[i+1]-x[i]
        # build matrix b
        dprimesB[i] = 6*((y[i+1]-y[i])/(x[i+1]-x[i])) - 6*((y[i]-y[i-1])/(x[i]-x[i-1]))
        
    # set endpoint values for spline 'natural condition' 
    dprimesA[0][0] = 1
    dprimesA[nPoints-1][nPoints-1] = 1
    dprimesA[1][0] = 0
    
    # solve for double primes using matrix Ax=b solver to get x
    dprimesX = np.linalg.solve(dprimesA,dprimesB)    
    
    #--------------------------------------------------------------
    
    xx = []  # for storing x interped values
    yy = []  # for storing y interped values
    
    t = sym.Symbol('t')    # define symbolic variable
    n = range(nPoints-1)   # index for interpolation 
    p = 0
    
    for i in n:   
        # Cubic Spline Method
        p0 = y[i]
        p1 = (y[i+1]-y[i])/(x[i+1]-x[i]) - (x[i+1]-x[i])*dprimesX[i+1]/6 - (x[i+1]-x[i])*dprimesX[i]/3
        p2 = dprimesX[i]/2
        p3 = (dprimesX[i+1]-dprimesX[i])/(6*(x[i+1]-x[i]))
        
        # build symbolic cubic function at each point
        p = p0 + p1*(t-x[i]) + p2*(t-x[i])**2 + p3*(t-x[i])**3
        
        fInterp = sym.lambdify(t,p)                             # turn symbolic p into real p
        xx0 = np.arange(x[i], x[i+1]+nxIncrement, nxIncrement)  # get interped values from x(i) to x(i+1)
        yy0 = [fInterp(z) for z in xx0]                         # get interped y values
        # write local interped x and y values to array
        # NOTE: duplicate data points exist at locations where piecewise function comes together
        # NOTE: but since they're equal it doesn't matter
        xx = np.concatenate((xx,xx0))
        yy = np.concatenate((yy,yy0))
                
    # plot
    plt.figure(4)
    plt.title('Cubic Spline Interpolation'+'\n'+ 'Interpolated in 5 Mev Steps')
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.plot(x,y,'bo', label='Data')
    plt.plot(xx,yy, label='Cubic Spline')
    plt.legend()   
    
    
#CubicSpline(E,Sig,9,5)
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------    
#*******************************************************************  
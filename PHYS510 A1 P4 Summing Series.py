"""
Created on Mon May 27 13:54:27 2013
PHYS 510, Assignment 1, Problem 4
"""

# Summing Series
"""
1. Write a program that plots the partial sum S(x,N) of the Taylor expansion 
   for exp(x) and its relative fractional error vs N up to N=30 for a given x.
   Test your program for x=10,2,-2, and -10. From the plots, explain why 
   this is not a good way to evaluate exp(x) when x<0. [Note: when x<0, 
   the Taylor expansion of the exponential is an alternating series.]
2. Modify your program so that it uses the identity exp(x)=1/exp(-x)~1/S(-x,N)
   to evaluate the exponential when x is negative. Explain why this 
   approach works better.
"""

import math
from matplotlib import pyplot

# Part 1
#*****************************************************************
# Define function to approximate exp(x) using a Taylor expansion
# Calculates partial sums and absolute fractional error 
# Then plots partial sums and errors vs n
def TaylorExp(x,n): 
    actual = math.exp(x)
    psum = 0
    psumArray = []
    relError = []
    iteration = range(n)
    
    for i in iteration:
        term = (float(x)**i) / math.factorial(i)  # compute each term
        psum = psum + term                 # keep a rolling partial sum
        
        err = abs((actual - psum)/actual)  # compute error abs(x-x' / x) 
        
        psumArray.append(psum)             # write values to arrays
        relError.append(err)
        
    # generate plot for partial sum    
    pyplot.subplot(2,1,1)
    pyplot.plot(iteration, psumArray, 'bo')
    pyplot.title('Taylor Approximation of Exp(' + str(x) + ')')
    pyplot.ylabel('Partial Sum')
    
    # plot actual exp(x) value as a line on partial sum plot
    pyplot.subplot(2,1,1)
    pyplot.plot([0, n-1], [actual, actual], 'r')
    
    # generate plot for errors
    pyplot.subplot(2,1,2)
    pyplot.plot(iteration, relError, 'ro')
    pyplot.xlabel('Iteration n')
    pyplot.ylabel('Relative Error')
#*****************************************************************


# Part 2
#*****************************************************************
# Define a function that modifies the Taylor approximation function above
# This function uses the identity described in the intro to evaluate the 
# exponential when x is negative
def ModTaylorExp(x,n): 
    actual = math.exp(x)
    psum = 0
    psumArray = []
    relError = []
    iteration = range(n)
    
    for i in iteration:
        term = (float(-x)**i) / math.factorial(i)  # modified so x = -x
        psum = psum + term                         # keep a rolling partial sum
        
        err = abs((actual - (1/psum))/actual) # modified so psum = 1/psum 
        
        psumArray.append(1/psum)              # modified so psum = 1/psum
        relError.append(err)
        
    # generate plot for partial sum    
    pyplot.subplot(2,1,1)
    pyplot.plot(iteration, psumArray, 'bo')
    pyplot.title('Modified Taylor Approximation of Exp(' + str(x) + ')')
    pyplot.ylabel('Modified Partial Sum')
    #pyplot.axis([0,30,-0.1,1])
    
    # plot actual exp(x) value as a line on partial sum plot
    pyplot.subplot(2,1,1)
    pyplot.plot([0, n-1], [actual, actual], 'r')
    
    # generate plot for errors
    pyplot.subplot(2,1,2)
    pyplot.plot(iteration, relError, 'ro')
    pyplot.xlabel('Iteration n')
    pyplot.ylabel('Relative Error')
#*****************************************************************

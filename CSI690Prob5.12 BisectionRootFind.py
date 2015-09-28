"""
Created on Sun Sep 27 17:46:41 2015
CSI 690, Assignment 3
Nick Crump
"""

# Problem 5.12: Bisection method to find max of equation

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

from sys import exit
from math import isnan,isinf
import numpy as np
import matplotlib.pyplot as plt


# input parameters
# -------------------------------------------------
# input function
def f(x):
    f  = -2*x**6 - 1.6*x**4 + 12*x + 1
    return f

# input function derivative
def df(x):
    df = -12*x**5 - 6.4*x**3 + 12
    return df

# input bracket parameters and stop criteria
xL   = 0       # lower bracket value
xH   = 1       # upper bracker value
pTol = 5       # percent error tolerance to stop
nMax = 20      # max iterations
# -------------------------------------------------


# error checking
notnum1 = isnan(df(xL)*df(xH))
notnum2 = isinf(df(xL)*df(xH))
if df(xL)*df(xH) >= 0:
    print '\nEITHER NO ROOTS OR MULTIPLE ROOTS EXIST BETWEEN BRACKET POINTS.'
    print 'ADJUST BRACKET VALUES.\n'
    exit()
if notnum1 == True or notnum2 == True:
    print '\nFUNCTION IS UNDEFINED AT BRACKET POINTS.'
    print 'ADJUST BRACKET VALUES.\n'
    exit()

# initialize variables
err = 100
ni  = 0
xi  = 0
xM  = 0

# print header
print '\n'
print 'iteration    x-value    f-max    relative error %'

# loop until error is less than input tolerance
while err > pTol:

    # stop if max iterations reached
    if ni >= nMax:
        print '\nMAX ITERATIONS REACHED'
        print 'Root =          ',xM
        print 'Iterations =    ',ni
        print 'Relative Error =',round(err,3),'%'
        exit()

    # bisection method
    xM = 0.5*(xL+xH)               # define midpoint

    if df(xL)*df(xM) > 0:
        xL = xM
        err = abs(xM-xi)/xM*100    # calculate approx error
        ni += 1                    # increment iteration
        xi = xM                    # store the n-1 midpoint

    elif df(xL)*df(xM) < 0:
        xH = xM
        err = abs(xM-xi)/xM*100    # calculate approx error
        ni += 1                    # increment iteration
        xi = xM                    # store the n-1 midpoint

    # print results
    print '    %i       %8.5f  %8.5f     %8.5f' % (ni,xM,f(xM),err)

# plot functions
x = np.arange(0,1.2,0.1)
n = len(x)
plt.figure()
plt.subplot(211)
plt.plot(x,f(x),'b',label='function')
plt.plot([xM,xM],[0,10],'k--')
plt.legend(loc=2)
plt.subplot(212)
plt.plot(x,df(x),'r',label='derivative')
plt.plot(x,np.zeros(len(x)),'k--')
plt.plot([xM,xM],[-20,15],'k--')
plt.legend(loc=3)
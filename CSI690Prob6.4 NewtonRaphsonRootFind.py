"""
Created on Sun Sep 27 17:46:41 2015
CSI 690, Assignment 3
Nick Crump
"""

# Problem 6.4: Newton-Raphson method to find roots of equation

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import numpy as np
import matplotlib.pyplot as plt


# input parameters
# -------------------------------------------------
# input function
def f(x):
    f  = -1 + 5.5*x -4*x**2 + 0.5*x**3
    return f

# input function derivative
def df(x):
    df = 5.5 - 8*x + 1.5*x**2
    return df

# input bracket parameters and stop criteria
x0   = 6       # initial guess
pTol = 0.01    # percent error tolerance to stop
nMax = 20      # max iterations
# -------------------------------------------------


# initialize variables
err = 100
ni  = 0
xi  = 0

# print header
print '\n'
print 'iteration      xi        f(xi)    relative error %'

# loop until error is less than input tolerance
while err > pTol:

    # stop if max iterations reached
    if ni >= nMax:
        print '\nMAX ITERATIONS REACHED'
        print 'Root =          ',xi
        print 'Iterations =    ',ni
        print 'Relative Error =',round(err,3),'%'
        exit()

    # Newton-Raphson method
    xi = x0 - (f(x0)/df(x0))   # predict root from derivative
    err = abs(xi-x0)/xi*100    # calculate approx error
    ni += 1                    # increment iteration
    x0 = xi                    # store the n-1 root


    # print results
    print '    %i       %8.5f   %8.5f     %8.5f' % (ni,xi,f(xi),err)

# plot functions
x = np.arange(-2,8,0.1)
n = len(x)
plt.figure()
plt.subplot(211)
plt.plot(x,f(x),'b',label='function')
plt.plot(x,np.zeros(len(x)),'k--')
plt.legend(loc=2)
plt.subplot(212)
plt.plot(x,df(x),'r',label='derivative')
plt.legend(loc=2)
"""
Created on Tue Oct 13 21:02:39 2015
CSI 690, Assignment 5
Nick Crump
"""

# Problem 14.8 Use steepest descent/ascent to minimize/maximize a 2D function

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

# input parameters
# -------------------------------------------------
# input function of 2-variables
def f(x,y):
    f = -8*x + x**2 + 12*y + 4*y**2 - 2*x*y
    return f

# input function gradient
def gradf(x,y):
    dfx = -8 + 2*x - 2*y
    dfy = 12 + 8*y - 2*x
    return dfx,dfy

# select whether finding max or min
opt   = 'min'  # select 'max' or 'min'

# input starting parameters
x0,y0 = 0,0    # initial x,y point
h     = 0.2    # step size
tol   = 1e-8   # stop tolerance
# -------------------------------------------------

# get sign for max/min selector
if opt == 'max': sign =  1
if opt == 'min': sign = -1

# get initial values
f0 = f(x0,y0)
dfx0,dfy0 = gradf(x0,y0)
err = 1
itr = 0

# print header
print '\niteration     x        y      f(x,y)  relative error'
print '    %i    %8.3f %8.3f  %8.3f    %8.6f' % (itr,x0,y0,f0,err)

# iterate to optimize function
while err > tol:
    # get points on line along gradient direction
    xi,yi = x0+sign*dfx0*h, y0+sign*dfy0*h
    # get new function value and error
    fi = f(xi,yi)
    err = abs((fi-f0)/fi)
    # update values
    dfx0,dfy0 = gradf(xi,yi)
    x0,y0 = xi,yi
    f0 = fi
    itr += 1

# print final results
print '   %i    %8.3f %8.3f  %8.3f    %8.6f' % (itr,x0,y0,f0,err)
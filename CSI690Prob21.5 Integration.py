"""
Created on Sat Oct 31 22:37:14 2015
CSI 690, Assignment 7
Nick Crump
"""

# Problem 21.5 Numerically integrate a function
# Uses Trapezoid, Simpson 1/3, Romberg, Gauss-Legendre methods

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import numpy as np
from scipy.special import p_roots


# input parameters
# -------------------------------------------------
# input function
def f(x):
    f = (4*x - 3)**3
    return f

# input limits of integration
a = -3
b =  5

# input number of intervals
n = 10
# -------------------------------------------------

# analytical solution
def intf(a,b):
    intf = ((4*b-3)**4 - (4*a-3)**4)/16
    return intf
asol = intf(a,b)


# trapezoid method
# -------------------------------------
def trapezoid(f,a,b,n):
    # get intervals
    h = (b-a)/float(n)
    x = np.arange(a,b+h,h)
    # sum intervals
    trap = f(x[0])+f(x[n])
    for i in range(1,n):
        trap += 2*f(x[i])
    trap = (h/2.0)*trap
    return trap
# -------------------------------------

# simpson 1/3 method
# -------------------------------------
def simpson(f,a,b,n):
    # intervals must be even for simpson 1/3
    if n%2 != 0:
        n  += 1
    h = (b-a)/float(n)
    x = np.arange(a,b+h,h)
    # sum intervals
    simp = f(x[0])+f(x[n])
    for i in range(1,n):
        if i%2 == 0: w = 2
        else:        w = 4
        simp += w*f(x[i])
    simp = (h/3.0)*simp
    return simp
# -------------------------------------

# romberg method
# -------------------------------------
def romberg(f,a,b,n):
    # initialize romberg matrix
    steps = []
    R      = np.zeros((n,n))
    R[0,0] = 0.5*(b-a)*(f(a)+f(b))
    # fill first column of romberg matrix
    for i in range(1,n):
        pts  = 2**i
        h    = (b-a)/float(pts)
        steps.append(h)
        psum = 0
        for j in range(1, pts):
            if j % 2 != 0:
                psum += f(a + j*h)
        # populate matrix
        R[i,0] = 0.5*R[i-1,0] + h*psum
    # get integral from romberg matrix
    for i in range(1,n):
        for j in range(i,n):
            r1 = R[j,i-1]
            r0 = R[j-1,i-1]
            k  = 1.0/(4.0**i-1)
            R[j,i] = r1 + k*(r1 - r0)
    romb = R[n-1,n-1]
    return romb
# -------------------------------------

# gauss-legendre method
# -------------------------------------
def gausslegendre(f,a,b,n):
    # uses roots of legendre polynomials
    x,w  = p_roots(n)
    psum = 0
    # valid only over interval [-1,1]
    # map input [a,b] to interval [-1,1]
    for i in range(n):
        # transform interval
        psum += w[i]*f(0.5*(b-a)*x[i] + 0.5*(b+a))
    gaus = 0.5*(b-a)*psum
    return gaus
# -------------------------------------


# call integration routines
trap = trapezoid(f,a,b,n)
simp = simpson(f,a,b,n)
romb = romberg(f,a,b,n)
gaus = gausslegendre(f,a,b,n)
# get integration error
terr = (abs(asol-trap)/asol)*100
serr = (abs(asol-simp)/asol)*100
rerr = (abs(asol-romb)/asol)*100
gerr = (abs(asol-gaus)/asol)*100


# print results
print '\nMethod            Solution    Error'
print '--------------------------------------'
print 'analytical     %12.6f   %6.3f'  % (asol,0),'%'
print 'trapezoid      %12.6f   %6.3f'  % (trap,terr),'%'
print 'simpson        %12.6f   %6.3f'  % (simp,serr),'%'
print 'romberg        %12.6f   %6.3f'  % (romb,rerr),'%'
print 'gauss-legendre %12.6f   %6.3f'  % (gaus,gerr),'%'
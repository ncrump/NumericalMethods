"""
Created on Tue Nov 03 16:34:18 2015
CSI 690, Assignment 8
Nick Crump
"""

# Problem 25.5 Numerically solve ODE
# Uses Euler, Heun, Midpoint, RK2 and RK4 methods

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

import numpy as np
import matplotlib.pyplot as plt


# input parameters
# -------------------------------------------------
# input ode as function
def func(x,y):
    f = y*np.sin(x)**3
    return f

# input initial conditions
y0 = 1    # initial condition on y
x0 = 0    # initial condition on x
xf = 3    # upper limit on x
h  = 0.1  # step size
# -------------------------------------------------

# analytical solution
def y(x):
    y = 1.947734*np.exp((np.cos(3*x)-9*np.cos(x))/12.0)
    return y

# simple euler method
# -------------------------------------
def simpleEuler(func,x0,y0,xf,h):
    # uses derivative at start of interval
    pts  = np.int(abs(xf-x0)/float(h))+1
    x    = np.linspace(x0,xf,pts)
    y    = np.zeros(pts)
    y[0] = y0
    # integrate ode
    for i in range(pts-1):
        # get startpoint derivative
        dy     = func(x[i],y[i])
        y[i+1] = y[i] + h*dy
    return x,y
# -------------------------------------

# improved euler method (Heun)
# -------------------------------------
def improvedEuler(func,x0,y0,xf,h):
    # uses average derivative at endpoints of interval
    pts  = np.int(abs(xf-x0)/float(h))+1
    x    = np.linspace(x0,xf,pts)
    y    = np.zeros(pts)
    y[0] = y0
    # integrate ode
    for i in range(pts-1):
        # get average derivative
        dy1    = func(x[i],y[i])
        dy2    = func(x[i]+h,y[i]+h*dy1)
        dyAvg  = 0.5*(dy1+dy2)
        y[i+1] = y[i] + h*dyAvg
    return x,y
# -------------------------------------

# modified euler method (Midpoint)
# -------------------------------------
def modifiedEuler(func,x0,y0,xf,h):
    # uses derivative at midpoint of interval
    pts  = np.int(abs(xf-x0)/float(h))+1
    x    = np.linspace(x0,xf,pts)
    y    = np.zeros(pts)
    y[0] = y0
    # integrate ode
    for i in range(pts-1):
        # get midpoint derivative
        xmid   = x[i]+0.5*h
        ymid   = y[i]+0.5*h*func(x[i],y[i])
        dyMid  = func(xmid,ymid)
        y[i+1] = y[i] + h*dyMid
    return x,y
# -------------------------------------

# 2nd-order runge-kutta (Ralston RK2)
# -------------------------------------
def rungekutta2(func,x0,y0,xf,h):
    # 2nd order accurate
    pts  = np.int(abs(xf-x0)/float(h))+1
    x    = np.linspace(x0,xf,pts)
    y    = np.zeros(pts)
    y[0] = y0
    # integrate ode
    for i in range(pts-1):
        y1     = func(x[i],y[i])
        y2     = func(x[i]+0.75*h,y[i]+0.75*h*y1)
        y[i+1] = y[i] + (h/3.0)*(y1 + 2*y2)
    return x,y
# -------------------------------------

# 4th-order runge-kutta (Standard RK4)
# -------------------------------------
def rungekutta4(func,x0,y0,xf,h):
    # 4th order accurate
    pts  = np.int(abs(xf-x0)/float(h))+1
    x    = np.linspace(x0,xf,pts)
    y    = np.zeros(pts)
    y[0] = y0
    # integrate ode
    for i in range(pts-1):
        y1     = func(x[i],y[i])
        y2     = func(x[i]+0.5*h,y[i]+0.5*h*y1)
        y3     = func(x[i]+0.5*h,y[i]+0.5*h*y2)
        y4     = func(x[i]+h,y[i]+h*y3)
        y[i+1] = y[i] + (h/6.0)*(y1 + 2*y2 + 2*y3 + y4)
    return x,y
# -------------------------------------


# call functions
x1,y1 = simpleEuler(func,x0,y0,xf,h)
x2,y2 = improvedEuler(func,x0,y0,xf,h)
x3,y3 = modifiedEuler(func,x0,y0,xf,h)
x4,y4 = rungekutta2(func,x0,y0,xf,h)
x5,y5 = rungekutta4(func,x0,y0,xf,h)
ya    = y(x1)

# plot results
plt.figure()
plt.plot(x1,ya,'o-',label='Analytical')
plt.plot(x1,y1,label='Euler')
plt.plot(x2,y2,label='Heun')
plt.plot(x3,y3,label='Midpoint')
plt.plot(x4,y4,label='RK2')
plt.plot(x5,y5,label='RK4')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc=2)
plt.annotate('h='+str(h),fontsize=16,xy=(0.20,0.54),xycoords='figure fraction')
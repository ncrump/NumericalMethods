"""
Example of fifth-order Runge-Kutta-Fehlberg method (RKF45) with 
adaptive step size to solve a first-order linear homogeneous ODE.

ODE:                y'(x) - y = x
Initial Condition:  y(1) = 3
Boundary:           x = [-2, 2] 
"""


import math as m
import numpy as np
import matplotlib.pyplot as plt


# ODE function f = y'(x) = x + y
#*********************************
def ODEfunc(xi, yi):
    fi = xi + yi
    
    return fi
#*********************************


# Fifth-Order Runge-Kutta-Fehlberg Method (RKF45) with adaptive step size
#*************************************************************************************
def RKF45(x0,xf,y0,h,eps):

    # create x and y arrays for storing solution values
    x = [x0]
    y = [y0]
    
    # set initial values
    xi = x0
    yi = y0
        
    # RKF45 coefficients
    # for intermediate function f1
    a1, b1, = 0.25, 0.25
    # for intermediate function f2
    a2, b2, c2 = 3.0/8.0, 3.0/32.0, 9.0/32.0
    # for intermediate function f3
    a3, b3, c3, d3 = 12.0/13.0, 1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0
    # for intermediate function f4
    a4, b4, c4, d4, e4 = 1, 439.0/216.0, 8, 3680.0/513.0, 845.0/4104.0
    # for intermediate function f5
    a5, b5, c5, d5, e5, f5 = 0.5, 8.0/27.0, 2, 3544.0/2565.0, 1859.0/4104.0, 11.0/40.0
    
    # for fifth-order solution
    Y1, Y2, Y3, Y4, Y5 = 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0
    # for error estimate 
    E1, E2, E3, E4, E5 = 1.0/360.0, 128.0/4275.0, 2197.0/75240.0, 1.0/50.0, 2.0/55.0
    
    # RKF45 method to solve ODE over each adaptive step h
    while xi <= xf:
        
        f0 = ODEfunc(xi,yi)
        
        xI = xi + a1*h
        yI = yi + b1*h*f0
        f1 = ODEfunc(xI,yI)
        
        xI = xi + a2*h
        yI = yi + b2*h*f0 + c2*h*f1
        f2 = ODEfunc(xI,yI)
        
        xI = xi + a3*h
        yI = yi + b3*h*f0 - c3*h*f1 + d3*h*f2
        f3 = ODEfunc(xI,yI)
        
        xI = xi + a4*h
        yI = yi + b4*h*f0 - c4*h*f1 + d4*h*f2 - e4*h*f3
        f4 = ODEfunc(xI,yI)
        
        xI = xi + a5*h
        yI = yi - b5*h*f0 + c5*h*f1 - d5*h*f2 + e5*h*f3 - f5*h*f4
        f5 = ODEfunc(xI,yI)
        
        
        yy = yi + h*(Y1*f0 + Y2*f2 + Y3*f3 - Y4*f4 + Y5*f5)        
        error = h*(E1*f0 - E2*f2 - E3*f3 + E4*f4 + E5*f5)
        RelError = abs(error/yy)
    
        hNew = h*(eps/RelError)**0.2
        if hNew >= h: 
            xi = xi + h  
            yi = yy
            x.append(xi)        
            y.append(yy)
        
        h = 0.9*hNew
                
    return x,y
#*************************************************************************************


# Main program that calls methods and plots solutions
#*************************************************************************************
# set up analytic solution of ODE y'(x)= y(x) + x for plotting
# analytic solution is y(x) = C*exp(x) - x - 1, for C = 5/exp(1) from y(1) = 3
# plotting over interval x = [-2,2] so y(-2) used as initial y0 in methods
x = np.linspace(-2,2,81)
y = []
c = 5.0/m.exp(1)
y0 = c*m.exp(-2) + 2 - 1

# get exact solution for y
for i in x:
    yActual = c*m.exp(i) - i - 1
    y.append(yActual)
        
# plot exact solution for comparison with method solutions
plt.title("Runge-Kutta-Fehlberg (RKF45) Solution to ODE $y$"+"'"+"$(x) - y = x$")
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x,y, linewidth=4, label='Exact Solution')

# call and plot method solution
xx,yy = RKF45(-2,2,y0,0.1,1.0e-5)

# plot method solutions
plt.plot(xx,yy, 'r', linewidth=1.5, label='RKF45')
plt.legend(loc=2)
#*************************************************************************************    
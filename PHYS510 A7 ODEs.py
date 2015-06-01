"""
Created on Wed Jul 10 20:41:20 2013
PHYS 510, Assignment 7

"""

# Ordinary Differential Equations (ODEs)
"""
Write a program to solve the differential equation y'(x) = y + x using the following methods:
    
    1. Simple Euler Method
    2. Modified Euler Method
    3. Improved Euler Method
    4. Fourth-Order Runge-Kutta Method (RK4)
    
Use x in the range [-2,2] with the initial condition y(1) = 3. 
Plot the results for step size h = 0.05, 0.10, 0.15, and 0.20, along with the exact result. 
Compare and discuss the four methods.
"""


import math as m
import numpy as np
import matplotlib.pyplot as plt

# Part 1
# Simple Euler Method
#*************************************************************************************
def simpleEuler(x0,xf,y0,h):
    # this solves the ODE y'(x)= y(x) + x with initial condition y(1) = 3
    # [x0,xf] = boundary of function domain, y0 = initial condition, h = step size
    # this method uses the derivative at the start of the interval across each step  
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Simple Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        y[i+1] = y[i] + h*(y[i] + x[i])  # Euler approx to function
        
    return x,y
#*************************************************************************************


# Part 2
# Modified Euler Method
#*************************************************************************************
def modifiedEuler(x0,xf,y0,h):
    # this solves the ODE y'(x)= y(x) + x with initial condition y(1) = 3
    # this method uses the derivative at the midpoint of the interval across each step
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1
    
    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Modified Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        xmid = x[i] + 0.5*h                # get x midpoint across interval
        ymid = y[i] + 0.5*h*(y[i] + x[i])  # get y at x midpoint
                
        y[i+1] = y[i] + h*(ymid + xmid)    # Euler approx to function at midpoint
        
    return x,y
#*************************************************************************************


# Part 3
# Improved Euler Method
#*************************************************************************************
def improvedEuler(x0,xf,y0,h):
    # this solves the ODE y'(x)= y(x) + x with initial condition y(1) = 3
    # this method uses the average derivative value over the interval across each step
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Improved Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        dy1 = y[i] + x[i]                  # get derivative at interval start
        dy2 = (y[i] + h*dy1) + (x[i] + h)  # get derivative at interval end
        dyAvg = 0.5*(dy1 + dy2)            # take average of derivatives
                
        y[i+1] = y[i] + h*(dyAvg)          # Euler approx to function using dyAvg
        
    return x,y
#*************************************************************************************


# Part 4
# Fourth-Order Runge-Kutta Method (RK4)
#*************************************************************************************
def RungeKutta(x0,xf,y0,h):
    # this solves the ODE y'(x)= y(x) + x with initial condition y(1) = 3
    # this method uses the fourth-order Runge-Kutta Method (RK4)
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Runge-Kutta Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        y1 = y[i] + x[i]
        y2 = (y[i] + 0.5*h*y1) + (x[i] + 0.5*h)
        y3 = (y[i] + 0.5*h*y2) + (x[i] + 0.5*h)
        y4 = (y[i] + h*y3) + (x[i] + h)
        
        y[i+1] = y[i] + (h/6.0)*(y1 + 2*y2 + 2*y3 + y4)
        
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
plt.title("Runge-Kutta (RK4) Solution to ODE $y$"+"'"+"$(x) - y = x$")
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x,y, linewidth=2.5, label='Exact Solution')

# call and plot method solutions for different step sizes below
h = [0.05, 0.10, 0.15, 0.20]
for i in h:
    #xi,yi = simpleEuler(-2,2,y0,i)
    #xi,yi = modifiedEuler(-2,2,y0,i)
    #xi,yi = improvedEuler(-2,2,y0,i)
    xi,yi = RungeKutta(-2,2,y0,i)
    
    # plot method solutions
    plt.plot(xi,yi, '.-', linewidth=1.5, label='h =' + str(i))
    plt.legend(loc=2)
#*************************************************************************************    
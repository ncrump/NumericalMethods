"""
Created on Thu Jun 13 20:47:04 2013
PHYS 510, Assignment 4 Parts 1,2,3
"""

# Numerical Integration
"""
This assignment is to develop algorithms for numerical integration using the following methods:
    1. Trapezoid
    2. Composite Simpson 1/3 
    3. Romberg 
Use the algorithms to evaluate the following integrals:
    1. Int[sin(x)*exp(x)], x = -pi to pi
    2. Int[x*exp(x**2)], x = -1 to 1
    3. Int[1/(1+x**2)], x= -4 to 4
Make a result table and compare your numerical result to the analytical solution by finding the relative 
percent error for each method.
"""

import numpy as np

# Part 1
# Trapezoid method
#*******************************************************************
def Trapezoid(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    h = (b-a)/n              # interval size (using n rather than n-1 for Python)
    x = np.arange(a,b,h)     # x array defined by interval size
    fx = [f(z) for z in x]   # y array of function values evaluated at x
    
    m = len(x)               # number of points to integrate over
    intSum = 0               # collects sum as integration iterates
    
    # loop over points to sum trapezoids in each interval to get integral
    for i in range(m):
        # check for using correct weight value in summation
        if i == 0 or i == n-1:
            w = 0.5
        else: 
            w = 1
            
        # Trapezoid routine
        intSum = intSum + fx[i]*h*w
        
    print 'Trapezoid Integral = ', intSum
        

#Trapezoid('sin(x)*exp(x)',-pi, pi, 50)
#Trapezoid('x*exp(x**2)',-1.0, 1.0, 50)        
#Trapezoid('1/(1+x**2)', -4.0, 4.0, 50)
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of points to evaluate
#-------------------------------------------------------------------    
#*******************************************************************    


# Part 2
# Simpson method
#*******************************************************************
def Simpson(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    # Simpson requires that n be odd
    # check if n is odd and if not then add 1
    if n % 2 == 0:
        n = n + 1
    
    h = (b-a)/n              # interval size (using n rather than n-1 for Python)
    x = np.arange(a,b,h)     # x array defined by interval size
    fx = [f(z) for z in x]   # y array of function values evaluated at x
    
    m = len(x)               # number of points to integrate over
    intSum = 0               # collects sum as integration iterates
    
    # loop over points to sum parabolas in each interval to get integral
    for i in range(m):
        # check for using correct weight value in summation
        if i == 0 or i == n-1:  # for endpoints use w = 1/3
            w = 1.0/3.0
        elif i % 2 != 0:        # for odd index use w = 4/3
            w = 4.0/3.0
        elif i % 2 == 0:        # for even index use w = 2/3
            w = 2.0/3.0
                        
        # Simpson routine
        intSum = intSum + fx[i]*h*w
        
    print 'Simpson Integral = ', intSum

        
#Simpson('sin(x)*exp(x)',-pi, pi, 50)
#Simpson('x*exp(x**2)',-1.0, 1.0, 50)        
#Simpson('1/(1+x**2)', -4.0, 4.0, 50)
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of points to evaluate
#------------------------------------------------------------------- 
#*******************************************************************


# Part 3
# Romberg method
#*******************************************************************
def Romberg(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    Rnm = np.zeros((n,n))    # create Romberg matrix
        
    # this part populates the initial column of the Romberg matrix 
    # using modified Trapezoid rule
    # ---------------------------------------------------
    # iterate over interval step size up to n
    for s in range(n):              
        
        h = (b-a)/(2**s)  # interval size
        pSum = 0          # partial sum in integral routine
        
        if s == 0:
            R = 0.5*(b-a)*(f(a)+f(b))  # R(0,0) of Romberg method
            Rnm[s][0] = R              # populate first value of matrix
            
        else: 
            # iterate over interval size up to 2**(n-1)
            for k in range(1, 2**s):
                if k % 2 != 0:                # only odds (k=1,3,5,..2**(n-1))
                    pSum = pSum + f(a + k*h)  # partial sum in integral routine
                  
            R = 0.5*Rnm[s-1][0] + h*pSum      # R(n,0) of Romberg method
            Rnm[s][0] = R                     # populate first column of matrix
    # ---------------------------------------------------
    
    # this part uses previous values in Romberg matrix to get new values
    # each new value is closer to the actual answer
    # --------------------------------------------------- 
    # iterate over rows and columns of Romberg matrix to get new values      
    for col in range(1,n):
        for row in range(col,n):
            romb1 = Rnm[row][col-1]
            romb0 = Rnm[row-1][col-1]
            coeff = 1/((4.0**col)-1)
                        
            # R(n,m) of Romberg method
            Rnm[row][col] = romb1 + coeff*(romb1 - romb0)
    # ---------------------------------------------------
    
    print 'Romberg Integral = ', Rnm[n-1][n-1]
    
     
#Romberg('sin(x)*exp(x)',-pi, pi, 10)
#Romberg('x*exp(x**2)',-1.0, 1.0, 10)        
#Romberg('1/(1+x**2)', -4.0, 4.0, 10)
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of interval steps up to 2**(n-1)
#-------------------------------------------------------------------    
#*******************************************************************         		
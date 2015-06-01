"""
Created on Thu May 30 21:01:03 2013
PHYS 510, Assignment 2 Functions
(Called by 'PHYS510 A2 Root Find.py' main program)
"""

# Root Finding
"""
The following functions are called for Assignment 2 in the root finding program.
1. Bisection method
2. Hybrid Bisection/Newton-Raphson method
3. Hybrid Bisection/Secant method
4. Hybrid Bisection/Muller-Brent method
"""

import math
import sympy

# 1
def rootBisection(g, f, xI, xF, Tol, nMax, Mthd):
# enter root finding algorithm by Bisection method
#*********************************************************************** 
        # initialize variables    
        error = 1
        n = 1
        found = 'no'        
        xiMid = 0  # initial midpoint value to store the n-1 value
        
        # loop until error is less than input tolerance
        while error > Tol:
            xMid = 0.5*(xI+xF)
            
            # set up main Bisection method:
            # make bracket interval smaller each iteration until root is found
            # check conditions and update bracket points     
            if f(xI)*f(xMid) > 0:
                xI = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
                
            elif f(xI)*f(xMid) < 0:
                xF = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint
                
            # break out of loop if root round    
            else:
                found = 'yes'
                break
        
        # output results to user
        if found == 'yes':
            print 'Exact Root Found =', xMid
            print 'Iterations = ', n 
            print '\n'          
        else:
            print 'Approximate Root =', xMid
            print 'Approximate Error =', error
            print 'Iterations = ', n-1
            print '\n'
            
# end rootBisection function
#***********************************************************************


# 2
def rootNewtonRaphson(g, f, xI, xF, Tol, nMax, Mthd):
# enter root finding algorithm by hybrid Bisection/Newton-Raphson method
#*********************************************************************** 
        # initialize variables    
        error = 1
        n = 1
        found = 'no'
        
        # get symbolic derivative of input function      
        x = sympy.Symbol('x')        # define x as symbolic variable
        sdf = sympy.diff(g, x)       # get symbolic derivative of input function g
        df = sympy.lambdify(x, sdf)  # turn symbolic derivative into numeric function
        
        # check condition for starting x value
        if f(xI) < f(xF): xn0 = xI
        else: xn0 = xF
        
        # loop until error is less than input tolerance
        while error > Tol:
            
            # set up main Newton-Raphson method:
            # use derivative of function as tangent line to get new x point 
            # calcuate new x value from f(x) and df(x)
            xn1 = xn0 - (f(xn0)/df(xn0))
            
            # if new x point is outside bracket interval then use midpoint instead
            if xn1 < xI or xn1 > xF:
                xn1 = 0.5*(xI+xF)
                
                # check conditions and update bracket points     
                if f(xI)*f(xn1) > 0:
                    xI = xn1
                    
                elif f(xI)*f(xn1) < 0:
                    xF = xn1
                
            error = abs(xn1 - xn0)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x
            
            # break out of loop if input max iterations met
            if n-1 >= nMax:
                found = 'yes'
                break
            
        # output results to user    
        if found == 'yes':
            print 'MAXIMUM ITERATIONS EXCEEDED'
            print 'Approximate Root =', xn1
            print 'Approximate Error =', error
            print 'Iterations =', n-1
            print '\n' 
        else: 
            print 'Approximate Root =', xn1
            print 'Approximate Error =', error
            print 'Iterations =', n-1
            print '\n' 
            
# end rootNewtonRaphson function
#***********************************************************************   
 
 
# 3 
def rootSecant(g, f, xI, xF, Tol, nMax, Mthd):
# enter root finding algorithm by hybrid Bisection/Secant method
#***********************************************************************        
        # initialize variables    
        error = 1
        n = 1
        found = 'no'
        
       # initialize starting values
        xn0 = xI
        xn1 = xF
        
        # loop until error is less than input tolerance
        while error > Tol:

            # set up main Secant method:
            # use secant line as linear approx to function to get new x point
            # calcuate new x value from secant line
            xn2 = xn1 - ((f(xn1)*(xn1 - xn0)) / (f(xn1) - f(xn0)))
            
            # if new x point is outside bracket  then use midpoint instead
            if xn2 < xI or xn2 > xF:
                xn2 = 0.5*(xI+xF)
                
                # check conditions and update bracket points     
                if f(xI)*f(xn2) > 0:
                    xI = xn2
                    
                elif f(xI)*f(xn2) < 0:
                    xF = xn2
                
            error = abs(xn2 - xn1)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x
            xn1 = xn2  # store the nth value for x
            
            # break out of loop if input max iterations met
            if n-1 >= nMax:
                found = 'yes'
                break
            
        # output results to user    
        if found == 'yes':
            print 'MAXIMUM ITERATIONS EXCEEDED'
            print 'Approximate Root =', xn2
            print 'Approximate Error =', error
            print 'Iterations =', n-1 
            print '\n' 
        else: 
            print 'Approximate Root =', xn2
            print 'Approximate Error =', error
            print 'Iterations =', n-1
            print '\n' 
            
# end rootSecant function            
#***********************************************************************


# 4
def rootMullerBrent(g, f, xI, xF, Tol, nMax, Mthd):
# enter root finding algorithm by hybrid Bisection/Muller-Brent method
#*********************************************************************** 
        # initialize variables    
        error = 1
        n = 1
        found = 'no'        
        
        # initialize starting values
        xn0 = xI
        xn1 = xF
        # use midpoint for 3rd point for use in quadratic approx to function
        xn2 = 0.5*(xn0+xn1)
        
        # loop until error is less than input tolerance
        while error > Tol:
            
            # set up main Muller-Brent method:
            # use 3 points as quadratic approx to function to get new x point
            # first calculate a, b, c from function values at 3 known points 
            c = f(xn2)        
            bNumer = (((xn0-xn2)**2)*(f(xn1)-f(xn2))) - (((xn1-xn2)**2)*(f(xn0)-f(xn2)))
            aNumer = ((xn1-xn2)*(f(xn0)-f(xn2))) - ((xn0-xn2)*(f(xn1)-f(xn2)))
            Denom = (xn0-xn1)*(xn0-xn2)*(xn1-xn2)
            b = bNumer / Denom  
            a = aNumer / Denom  
            
            # then calculate new x value from quadratic eq
            # check b and use alternate quadratic eq to avoid subtractive cancellations
            if b >= 0:
                xn3 = xn2 - (2*c) / (b + math.sqrt(b**2 - 4*a*c))
            else: 
                xn3 = xn2 + (2*c) / (-b + math.sqrt(b**2 - 4*a*c))
            
            # if new x point is outside bracket interval then use midpoint instead
            if xn3 < xI or xn3 > xF:
                xn3 = 0.5*(xI+xF)
                
                # check conditions and update bracket points     
                if f(xI)*f(xn3) > 0:
                    xI = xn3
                    
                elif f(xI)*f(xn3) < 0:
                    xF = xn3
                
            error = abs(xn3 - xn2)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x
            xn1 = xn2  # store the nth value for x
            xn2 = xn3  # store the n+1 value for x
            
            # break out of loop if input max iterations met
            if n-1 >= nMax:
                found = 'yes'
                break
            
        # output results to user    
        if found == 'yes':
            print 'MAXIMUM ITERATIONS EXCEEDED'
            print 'Approximate Root =', xn3
            print 'Approximate Error =', error
            print 'Iterations =', n-1
            print '\n'
        else: 
            print 'Approximate Root =', xn3
            print 'Approximate Error =', error
            print 'Iterations =', n-1
            print '\n'
        
# end rootMullerBrent function
#*********************************************************************** 
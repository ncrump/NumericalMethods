"""
Created on Thu May 30 21:01:03 2013
PHYS 510, Assignment 2 
(Calls 'rootFindFunctions.py' where algorithms are defined)
"""

# Root Finding
"""
Design a program that finds the roots of a given function so the user 
has the following options:
1. the user can enter root brackets between which a single root is sought
2. the user has the choice of the root-finding algorithm used
    * Bisection method
    * Hybrid Bisection / Newton-Raphson method
    * Hybrid Bisection / Secant method
    * Hybrid Bisection / Muller-Brent method
3. the user should be able to enter the function whose roots are to be found
4. the program should report the number of iterations
"""

import math
import rootFindFunctions

# define root finding function with choices for algorithm used
#--------------------------------------------------------------------------------------
# called as rootFind('function(x)', xInitial, xFinal, Tolerance, MaxIterations, Method)
# 'function(x)' = input function of x inside quotes (Ex: 'cos(x)-x')
# xInitial = initial x bracket value on one side of the root (entered as decimal)
# xFinal = initial x bracket value on the other side of the root (entered as decimal)
# Tolerance = the degree of accuracy in which to compute the root (Ex: 10e-5)
# MaxIterations = max number of iterations to compute if Tolerance has not been met
# Method = choice of root finding method to use, entered inside quotes
         # (options are: 'bisection', 'newton', 'secant', 'muller')
#-------------------------------------------------------------------------------------- 

def rootFind(g, xI, xF, Tol, nMax, Mthd):
    f = lambda x: eval(g)  # create function from input
    
    # check input bracket condition to ensure root is betwwen braket points
    # check to make sure function is defined at bracket points too
    # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
    notnum1 = math.isnan(f(xI)*f(xF))
    notnum2 = math.isinf(f(xI)*f(xF))
    if f(xI)*f(xF) >= 0 or notnum1 == True or notnum2 == True: 
        print '\nNO ROOTS BETWEEN THESE POINTS OR FUNCTION IS UNDEFINED.'
        print 'PLEASE ADJUST BRACKET VALUES.\n'
        return
    
    if Mthd == 'bisection':
        # enter root finding algorithm by Bisection method
        rootFindFunctions.rootBisection(g, f, xI, xF, Tol, nMax, Mthd)

    elif Mthd == 'newton':
       # enter root finding algorithm by hybrid Bisection/Newton-Raphson method
       rootFindFunctions.rootNewtonRaphson(g, f, xI, xF, Tol, nMax, Mthd)

    elif Mthd == 'secant':
        # enter root finding algorithm by hybrid Bisection/Secant method
        rootFindFunctions.rootSecant(g, f, xI, xF, Tol, nMax, Mthd)

    elif Mthd == 'muller':
        # enter root finding algorithm by hybrid Bisection/Muller-Brent method
        rootFindFunctions.rootMullerBrent(g, f, xI, xF, Tol, nMax, Mthd)


    else:
        print '\nINVALID ENTRY. PLEASE ADJUST INPUT AND TRY AGAIN.\n'
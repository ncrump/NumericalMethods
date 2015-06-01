"""
Created on Sun May 26 09:49:41 2013
PHYS 510, Assignment 1, Problem 3
"""

# The Quadratic Formula and Subtractive Cancellation
"""
1. Write a program that calculate the four solutions x1, x2, y1, and y2 
   for the quadratic equation with a set of arbitrary values of a, b, and c.
2. Take a=b=1, and c=10^-n, n=1,2,3,... and compute the four solutions for 
   increasing n. Make sure that your output have enough significant figures.
3. What is the largest n when machine precision causes Eq. (2) to fail?
4. Since b>0, subtractive cancellation will NOT occur in y1 and x2. Treat 
   these two values as the accurate solutions to the quadratic equation and 
   calcuate the errors in the other two solutions x1 and y2. Plot the fractional 
   relative errors of x1 and y2 as a function of 4ac in a log-log graph.
5. Comment on your results in relation to the machine precision.
"""

import math
from matplotlib import pyplot

# Part 1
#**************************************************************
# Define functions for quadratic eq and alternate quadratic eq
def quadraticEq(a,b,c):
    x1 = (-b + math.sqrt(b**2 - 4*a*c)) / (2*a)
    x2 = (-b - math.sqrt(b**2 - 4*a*c)) / (2*a)
    return x1, x2
    
def AltquadraticEq(a,b,c):
    y1 = (-2*c) / (b + math.sqrt(b**2 - 4*a*c))
    y2 = (-2*c) / (b - math.sqrt(b**2 - 4*a*c))
    return y1, y2
#**************************************************************
    
    
# Parts 2 & 3
a = 1
b = 1
n = 1 

x1Error = []
y2Error = []
array4ac = []
# By trail and error it was found that AltquadraticEq fails for n > 16
# This is due to division by zero as c becomes too small to represent 
while n < 17:
    print 'n =', n
    c = 10**-n
    x1, x2 = quadraticEq(a,b,c)
    y1, y2 = AltquadraticEq(a,b,c) 
    
    print 'x1 =', x1
    print 'x2 =', x2
    print 'y1 =', y1
    print 'y2 =', y2
    
    n = n + 1
    
   # Store the relative errors for x1 and y2 for plotting
    x1Error.append(abs((y1 - x1)/y1))  # relative error abs(x-x' / x)
    y2Error.append(abs((x2 - y2)/x2))  # relative error abs(x-x' / x)
    array4ac.append(4*a*c)
    
    
# Part 4
# Plot x1 and y2 relative errors against 4ac on a log-log plot
pyplot.plot(array4ac, x1Error, 'bo', array4ac, y2Error, 'ro')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.title('Relative Errors vs 4ac')
pyplot.xlabel('4ac')
pyplot.ylabel('Relative Error')  
pyplot.legend(('x1Error','y2Error'))     
    
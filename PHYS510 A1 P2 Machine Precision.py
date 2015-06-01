"""
Created on Sat May 25 23:38:45 2013
PHYS 510, Assignment 1, Problem 2
"""
 
# Determine your Machine Precision
"""
Write a program to determine the machine precision of your machine
(within a factor of 2) for double precision floating-point operations.
The output of your program should contain at least two columns showing
the value of x and eps as eps gets increasing small. 
What is the value of eps for your machine?
"""

eps = 1
x = 2
i = 0

# Iterate and make epsilon smaller by half each iteration
# Then add to 1 until epsilon become small enough where 1+eps = 1
while x > 1:
    eps = eps/2.0
    x = 1.0 + eps
    
    # use formatted print statement to get enough decimals of precision
    print '%2.0f' % i, '%0.16f' % eps, '%0.16f' % x  
    
    i = i + 1

print '\nIteration =', i-1
print 'eps =', eps
print 'x =', x


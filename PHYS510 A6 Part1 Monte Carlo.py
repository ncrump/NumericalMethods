"""
Created on Thu Jun 27 20:48:37 2013
PHYS 510, Assignment 6 Part 1
"""

# Monte Carlo Methods
"""
Part 1: Using the Monte Carlo approach, approximate the volume of a ball with radius 0.04 at 
        x,y,z = (1/3 ,1/3, 1/2). Use 10e6 three dimensional points with seed x0 = 1. How close 
        is your approximation to the correct answer?
"""


import math
import numpy as np
            
# Part 1
# Approximate Volume of Sphere using Monte Carlo Method
#*******************************************************************
def Sphere(r,x0,y0,z0,n):
    # r = radius of sphere, x0,y0,z0 = center of sphere, n = number of points
    # calculates volume of a sphere using the method of generating random points
    # within a cube and counting how many fall within the sphere bounded by the cube
    # Volume of Sphere = (Volume of cube)*[(# points in sphere)/(# points in cube)]
    
    # get bounding box of cube    
    xmin, xmax = x0 - r, x0 + r
    ymin, ymax = y0 - r, y0 + r
    zmin, zmax = z0 - r, z0 + r
    rSqr = r**2
    
    inSphere = 0  # start a counter for points landing in sphere
    
    # generate n random (x,y,z) points within cube 
    # use sphere equation to determine if within sphere
    for i in range(int(n)):
        x = np.random.uniform(xmin,xmax)
        y = np.random.uniform(ymin,ymax)
        z = np.random.uniform(zmin,zmax)
        
        if (x-x0)**2 + (y-y0)**2 + (z-z0)**2 < rSqr:
            inSphere = inSphere + 1
            
    AreaCube = (2*r)**3                 # area of cube
    AreaSphere = AreaCube*(inSphere/n)  # approx area of sphere
    Actual = (4.0/3.0)*math.pi*(r**3)   # actual area of sphere
    
    print 'Approximate Area = ', AreaSphere
    print 'Accurate to ', (abs(Actual-AreaSphere)/Actual)*100, '%'
    
    
Sphere(0.04,1/3,1/3,1/2,10e6)  
#*******************************************************************
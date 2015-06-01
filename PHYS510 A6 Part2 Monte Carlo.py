"""
Created on Thu Jun 27 20:48:37 2013
PHYS 510, Assignment 6 Part 2
"""

# Monte Carlo Methods
"""
Part 2: Simulate a beam of particles scattering off of an arbitrary target. We will assume the
        beam initially moves along the x-axis. At some negative value x = -L, far from the target
        (at x = 0), the beam has a circular cross section with radius R. Assume the particles are
        uniformly distributed throughout the circle. We'll assume the momentum of the particles
        form a Gaussian distribution with a mean and standard deviation. As they propagate, the
        target exerts a force F(r) on each particle at distance r. The goal is to find the trajectory
        of the particles and get the position distribution of the particles when they hit a large
        flat detector at x = L. 
"""


import numpy as np   
import matplotlib.pyplot as plt          
from mpl_toolkits.mplot3d import Axes3D  # used for 3D plotting

# Part 2
# Simulate a beam of particles using the Monte Carlo Method
#*******************************************************************
def ParticleSim(N,m,C,L,R,Pmu,Psig,dt):
    
    # initialize trajectory plot to create 3D scatter plot of particle positions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # create random uniform distribution of values for initial y and z locations
    # this creates random points within the box from -R to R of N points
    # initial x starts at -L
    y00 = np.random.uniform(-R,R,N)
    z00 = np.random.uniform(-R,R,N)
    
    # create arrays to store initial y,z points that fall within the cross section
    x0 = -L
    y0 = []
    z0 = []
    
    # loop over initial y,z points, check if they fall within the circle of radius R
    # the yz-axis is in the plane of the beam cross section
    # only keep y,z points within the beam circular cross section
    Rsqr = R**2
    for i in range(N):
        if y00[i]**2 + z00[i]**2 < Rsqr:
            y0.append(y00[i])
            z0.append(z00[i])
            
    # define new total number of particles after keeping only those within the beam            
    N = len(y0)
    
    # create random normal distribution of values for inital particle momentums
    # this uses a mean of Pmu and standard dev of Psig
    P0 = np.random.normal(Pmu,Psig,N)
    
    # create arrays to store final y,z values that reach the detector plate at x=L
    yFinal = []
    zFinal = []
    
    # loop over each of the N particles to generate new x,y,z positions of trajectory
    for i in range(N):
        # get initial x,y,z position of ith particle
        xp = x0
        yp = y0[i]
        zp = z0[i]
        
        # get x,y,z momentum of ith particle
        # initial y and z components of momentum are zero
        vpx = P0[i]/m
        vpy = 0
        vpz = 0
        
        # set inital time to zero 
        t = 0
        
        # create arrays for storing x,y,z positions of particles for trajectory plot
        xarr = [xp]
        yarr = [yp]
        zarr = [zp]
        
        # for a single ith particle this is where the time and position stepping is
        # x,y,z Force is calculated at time t to update position at time t+dt
        # x,y,z Position in then updated for time t+dt
        # x,y,z Velocity is then updated for time t+dt
        # then time is stepped to the next t+dt
        # the loop continues as long as xPosition < L and xVelocity is positive
        # this ensures that particles deflected back are not tracked
        while xp < L and vpx > 0:
            # calcuate x,y,z force using the Coulomb force: F=kqQ/r**2
            Fx = (C*xp)/(xp**2 + yp**2 + zp**2)**(3.0/2.0)
            Fy = (C*yp)/(xp**2 + yp**2 + zp**2)**(3.0/2.0)
            Fz = (C*zp)/(xp**2 + yp**2 + zp**2)**(3.0/2.0)
            
            # take a time step            
            t = t + dt
            
            # calculate x,y,z position using standard kinematic equations
            # x = x0 + v0t + 0.5at**2 (where a = Fx/m)
            # this assumes constant acceleration between t and t+dt (approximate)
            x = xp + vpx*t + (0.5/m)*Fx*(t**2)
            y = yp + vpy*t + (0.5/m)*Fy*(t**2)
            z = zp + vpz*t + (0.5/m)*Fz*(t**2)
            
            # calculate x,y,z velocity using v = v0 + at (where a = F/m)
            vx = vpx + (1.0/m)*Fx*t
            vy = vpy + (1.0/m)*Fy*t
            vz = vpz + (1.0/m)*Fz*t
            
            # store positions to arrays for trajectory plotting
            xarr.append(x)
            yarr.append(y)
            zarr.append(z)
                        
            # update values for next iteration (prev values = new values)
            xp = x
            yp = y
            zp = z
            vpx = vx
            vpy = vy
            vpz = vz
    
        # check if particle hits detector plate and store y,z position to array if so
        # I had to give the plate a thickness since no particle lands exactly at x=L
        # this is due to a fixed time step and particles arriving at different times
        if L <= x <= L+5 and -L <= y <= L and -L <= z <= L: 
            yFinal.append(y)
            zFinal.append(z)
            
        # plot 3D scatter of ith particle trajectory
        ax.scatter(xarr,yarr,zarr)
        
    # finish and show trajectory plot after all particles have been plotted
    # this is inside particle for loop but outside position stepping while loop
    plt.title('Particle Trajectories')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
        
    # plot 3D histogram after all particle trajectories have been calculated
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # create 2D histogram of final y,z positions (returns arrays for plotting)
    hist,xedges,yedges = np.histogram2d(yFinal,zFinal,bins=10)
        
    # set plot increments
    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1])
    
    # using 3D bar plot for making 3D histogram of x,y,z positions
    # each bar needs a x,y,z position and x,y,z width
    # had to flatten 2D hist array to make into 1D array for plotting
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 15*np.ones_like(zpos)
    dy = 15*np.ones_like(zpos)
    dz = hist.flatten()
    
    # make and show 3D bar plot
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
    ax.set_title('Particle Distribution at Detector Plate')
    ax.set_xlabel('y')
    ax.set_ylabel('z')
    plt.show()       
            

ParticleSim(1000,1.0,1.0,100,0.1,5.0,0.5,0.01)
#-------------------------------------------------------------------
# N = number of particles 
# m = mass of particle
# C = kqQ constant in Coulomb force
# L = length between origin and detector plate 
# R = radius of particle beam cross section
# Pmu = mean of Gaussian distribution describing initial particle momentum
# Psig = standard dev of Gaussian distribution describing initial momentum
# dt = time step
#------------------------------------------------------------------- 
#*******************************************************************
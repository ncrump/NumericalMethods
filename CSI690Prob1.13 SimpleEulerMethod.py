"""
Created on Tue Sep 01 21:44:01 2015
CSI 690, Assignment 1
Nick Crump
"""

# Problem 1.13: Euler method to solve simple 1st-order ODE

"""
From Numerical Methods for Engineers - Chapra 6th Ed
"""

from math import pi
import matplotlib.pyplot as plt

# set input parameters
k  = 0.1    # mm/min
r0 = 3.0    # mm
t0 = 0.0    # min
tf = 20.0   # min
dt = 0.25   # min

# get initial values
N  = int((tf-t0)/dt)
V0 = (4/3.0)*pi*r0**3
ti = t0

# print initial values
print '%4.2f %10.6f %10.6f' % (t0,V0,r0)

# loop to get numerical solution by simple Euler method
V = [V0]
t = [t0]
for i in range(N):
    Vi = V[i] - 4*pi*k*dt*((3*V[i])/(4*pi))**(2/3.0)
    ri = ((3*Vi)/(4*pi))**(1/3.0)
    ti = ti + dt
    V.append(Vi)
    t.append(ti)
    print '%4.2f %10.6f %10.6f' % (ti,Vi,ri)


# loop to get analytical solution done by hand
Vact = []
for i in t:
    Vt = (V0**(1/3.0) - (4*pi/3.0)**(1/3.0)*k*i)**3
    Vact.append(Vt)


# make plots
plt.plot(t,Vact,'r.-',label='Analytical Solution')
plt.plot(t,V,'b.-',label='Numerical Solution')
plt.xlabel('time (min)',fontsize=16)
plt.ylabel('volume (mm$^3$)',fontsize=16)
plt.legend()
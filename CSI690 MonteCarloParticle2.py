"""
Created on Fri Nov 27 22:50:06 2015
CSI 690, Assignment 11
Nick Crump
"""

# Monte Carlo simulation to compute average energy of single particle in 1D

import numpy as np
import matplotlib.pyplot as plt


# input parameters
# --------------------------------------
dmax = 10      # max step size
tmp0 = 0.1     # start temp
tmp1 = 1.0     # final temp
ntmp = 10      # number of temps
Etol = 0.1     # energy error tolerance
# --------------------------------------

# set temperatures
Tarr = np.linspace(tmp0,tmp1,ntmp)

# initialize arrays
Earr = []
Ethr = []

# loop over temps
for i in range(ntmp):
    # initialize for this temp
    T    = Tarr[i]
    x    = 100
    E    = 0.1*x*x
    aveE = E
    thrE = 0.5*T
    errE = 1
    acc  = 0
    step = 0
    # loop over atom
    while errE > Etol:
        # store old values
        xi = x
        Ei = E
        # take a step
        x += np.random.uniform(-dmax,dmax,1)[0]
        # calculate new energy
        E = 0.1*x*x
        d = E-Ei
        # accept move if new energy is lower
        if d < 0:
            acc += 1
        # otherwise check boltzmann condition
        else:
            prob = np.exp(-d/T)
            # accept move with this probability
            if np.random.uniform(0,1,1)[0] < prob:
                acc += 1
            # reject move otherwise
            else:
                x = xi
                E = Ei
        # increment and collect running average
        step  += 1
        aveEi = aveE
        aveE  = aveEi + (E-aveEi)/float(step)
        errE  = abs(thrE-aveE)/thrE

    # store values per temp
    accP = acc/float(step)
    Earr.append(aveE)
    Ethr.append(thrE)

    # print results
    if i == 0: print '\n T     <E>   %Err   %Acc    Steps'
    print '%4.2f %6.2f  %4.2f  %4.2f  %i' % (T,aveE,errE*100,accP*100,step)

# make plots
plt.figure()
plt.plot(Tarr,Earr,'bo-',label='Simulated')
plt.plot(Tarr,Ethr,'ro-',label='Theoretical')
plt.title('Monte Carlo Particle Sim')
plt.xlabel('T')
plt.ylabel('<E>')
plt.legend(loc=2)
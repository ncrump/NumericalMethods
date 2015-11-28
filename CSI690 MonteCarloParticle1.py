"""
Created on Fri Nov 27 22:50:06 2015
CSI 690, Assignment 11
Nick Crump
"""

# Monte Carlo simulation to compute average energy of single particle in 1D

import numpy as np
import matplotlib.pyplot as plt


# input parameters
# ----------------------------------
dmax = 10      # max step size
tmp0 = 0.1     # start temp
tmp1 = 1.0     # final temp
ntmp = 10      # number of temps
runs = 1000    # iterations per temp
# ----------------------------------

# set temperatures
Tarr = np.linspace(tmp0,tmp1,ntmp)

# get uniform random numbers
rand = np.random.uniform(-dmax,dmax,[ntmp,runs])
chek = np.random.uniform(0,1,[ntmp,runs])

# initialize arrays
Earr = []
Ethr = []
x    = 100
# loop over temps
for i in range(ntmp):
    # initialize for this temp
    T    = Tarr[i]

    E    = 0.1*x*x
    Esum = E
    acc  = 0
    # loop over atom
    for j in range(runs):
        # store old values
        xi = x
        Ei = E
        # take a step
        x += rand[i,j]
        # calculate new energy
        E = 0.1*x*x
        d = E-Ei
        # accept move if new energy is lower
        if d < 0:
            Esum += E
            acc  += 1
        # otherwise check boltzmann condition
        else:
            prob = np.exp(-d/T)
            # accept move with this probability
            if chek[i,j] < prob:
                Esum += E
                acc  += 1
            # reject move otherwise
            else:
                x     = xi
                E     = Ei
                Esum += E

    # store average energy and acceptance probability per temp
    aveE = Esum/float(runs+1)
    accP = (acc/float(runs+1))*100
    Earr.append(aveE)

    # theoretical result and error
    thrE = 0.5*T
    errE = (abs(thrE-aveE)/thrE)*100
    Ethr.append(thrE)

    # print results
    if i == 0: print '\n T     <E>   %Err   %Acc   Steps'
    print '%4.2f %6.2f  %4.2f  %4.2f  %i' % (T,aveE,errE,accP,runs)

# make plots
plt.figure()
plt.plot(Tarr,Earr,'bo-',label='Simulated')
plt.plot(Tarr,Ethr,'ro-',label='Theoretical')
plt.title('Monte Carlo Particle Sim')
plt.xlabel('T')
plt.ylabel('<E>')
plt.legend(loc=2)
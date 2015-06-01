"""
Example of how to use the matplotlib 3D module 
for making 3D histograms out of 3D bar plots.
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xpos=np.random.randint(1,10,10)
ypos=np.random.randint(1,10,10)
zpos=np.zeros(10)
dx=np.ones(10)
dy=np.ones(10)
dz=np.random.randint(1,10,10)

ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

plt.show()


"""

"""

## Import modules
from __future__ import print_function                  # make print function work in python 2.x
import numpy as np		                               # import numpy
import matplotlib.pyplot as plt 			           # plotting tools
from mpl_toolkits.mplot3d import Axes3D		           # plotting tools
import sys                                             # progress messages
import running as start                                # startup message
from new_nodal_positions import new_nodal_pos          # data export for physcial quantities

## Global settings
np.set_printoptions(threshold='nan')		# Do not truncate print

## User input
# Fix Python 2.x.
try: input = raw_input
except NameError: pass
sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.5
epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 1
T = input('Temperature, expressed in epsilon (default: 1): ') or 1
number_of_beads = input('Amount of beads per polymer: ') or 20
plot_data = input('Plot data? (y/n, default: y): ') or 'y'


## Fixed parameters
angle_dof = 10                              # Amount of different angles the polymer can move in
angles = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices

## Message at simulation start
start.message()

nodal_points = np.zeros((number_of_nodal_points,2),dtype=float)  # initialize all nodal point positions
posssible_nodal_points = np.zeros((len(angles),2),dtype=float)   # initialize list for all possible positions of the next nodal point


for N in range(0, number_of_nodal_points-1):



    possible_nodal_points = new_nodal_pos(nodal_points[N,:],angles)  # calculate all possible nodal points
    nodal_points[N+1,:] = possible_nodal_points[0,:]                              # choose first of all possible nodal points, so polymer is a random 2D walk at the moment. FOR TESTING ONLY!!
    plot_nodal_points = np.zeros((N+1,2),dtype=float)                             # this block is used to plot the polymer as it grows, only for tesing purposes
    plot_nodal_points = nodal_points[0:(N+1)]

    #energies = calculate_energies(nodal_points[N,:],    )
    #plt.scatter(plot_nodal_points[:,0],plot_nodal_points[:,1])
    #plt.plot(plot_nodal_points[:,0],plot_nodal_points[:,1])
    #plt.axis([-7, 7, -7, 7])
    #plt.show()

plt.scatter(plot_nodal_points[:,0],plot_nodal_points[:,1])
plt.plot(plot_nodal_points[:,0],plot_nodal_points[:,1])
plt.axis([-12, 12, -12, 12])
plt.show()
"""

"""

## Import modules
from __future__ import print_function       # make print function work in python 2.x
import numpy as np		                    # import numpy
import matplotlib.pyplot as plt 			# plotting tools
from mpl_toolkits.mplot3d import Axes3D		# plotting tools
import sys                                  # progress messages
import running as start                     # startup message

## Global settings
np.set_printoptions(threshold='nan')		# Do not truncate print

## User input
# Fix Python 2.x.
try: input = raw_input
except NameError: pass
sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.5
epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 1
T = input('Desired temperature (default: 1): ') or 1
plot_data = input('Plot data? (y/n, default: y): ') or 'y'
chain_length = input('Amount of nodal points: ') or 1000

## Message at simulation start
start.message()

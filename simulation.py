# Import libraries and modules
import numpy as np		                            # import numpy  
# import sys                                          # progress messages
from new_beads_positions import new_beads_pos       # calculate possible new bead positions
from calculate_energies import calculate_energies   # calculate energies for each new possible bead position
from calculate_energies2 import calculate_energies2   # calculate energies for each new possible bead position
from determine_new_bead import determine_new_bead   # function used to determine the final bead position by comparing the boltzmann factors
import list_tracking as track

def user_input():
    sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    T = input('Temperature, expressed in epsilon (default: 1): ') or 1
    number_of_beads = input('Amount of beads per polymer: ') or 150
    plot_data = input('Plot data? (y/n, default: y): ') or 'y'
    return float(sigma), float(epsilon), float(T), int(number_of_beads), plot_data

def start(number_of_beads,sigma,epsilon,T):
    ## Fixed parameters
    angle_dof = 6                              # Amount of different angles the polymer can move in
    angles1 = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    angles2 = np.linspace(0,2*np.pi,angle_dof)

    beads_pos = np.zeros((number_of_beads,3),dtype=float)  # initialize all bead positions
    possible_beads_pos = np.zeros((len(angles1)*len(angles2),3),dtype=float)   # initialize list for all possible positions of the next bead

    sigma_squared = sigma*sigma

    track.init(int(np.floor(number_of_beads/2)))

    for N in range(0, number_of_beads):
        possible_beads_pos = new_beads_pos(beads_pos[N-1,:],angles1,angles2)  # calculate all possible nodal points
        possible_beads_pos = possible_beads_pos.reshape(-1,3)
        energies = calculate_energies(possible_beads_pos,beads_pos[0:N],epsilon,sigma_squared)
        new_bead_index = determine_new_bead(energies,T)            # determine final new bead
        beads_pos[N,:] = possible_beads_pos[new_bead_index,:]    # add new final new bead to the polymer
        track.store(beads_pos[N,:],N)
    return beads_pos

def plot(beads_pos):
    import matplotlib.pyplot as plt 			           # plotting tools
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = pylab.figure()
    ax = Axes3D(fig)
    plt.plot(beads_pos[:,0], beads_pos[:,1],beads_pos[:,2])
    plt.show()
    return
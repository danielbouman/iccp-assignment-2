# Import libraries and modules
import numpy as np		                            # import numpy  
# import sys                                          # progress messages
from new_beads_positions import new_beads_pos       # calculate possible new bead positions
from determine_new_bead import determine_new_bead   # function used to determine the final bead position by comparing the boltzmann factors
import lj_energy
import matplotlib.pyplot as plt     # plotting tools

def user_input():
    # sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    # epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    # T = input('Temperature, expressed in epsilon (default: 1): ') or 1.0
    # number_of_beads = input('Amount of beads per polymer: ') or 150
    # plot_data = input('Plot data? (y/n, default: y): ') or 'y'
    sigma = 0.8
    epsilon = 0.25
    bending_energy = 10
    T = 1
    number_of_beads = 150
    plot_data = 'n'
    return float(sigma), float(epsilon), float(T), int(number_of_beads), plot_data

def start(number_of_beads,sigma,epsilon,T):
    ## Fixed parameters
    angle_dof = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    
    sigma_squared = sigma*sigma
    existing_pos = np.zeros((number_of_beads,2),dtype=float)    # initialize all bead positions
    candidate_pos = np.zeros((len(angles),2),dtype=float)       # initialize list for all possible positions of the next bead
    angle_last_bead = 0
    
    for N in range(1, number_of_beads):
        candidate_pos,angle_updated = new_beads_pos(existing_pos[N-1,:],angles)  # calculate all possible nodal points
        energies = lj_energy.func(existing_pos[0:N,:],candidate_pos,sigma_squared,epsilon,angle_dof,N,angle_updated,angle_last_bead,bending_energy) # calculate energies
        new_bead_index = determine_new_bead(energies,T)         # determine final new bead
        existing_pos[N,:] = candidate_pos[new_bead_index,:]    # add new final new bead to the polymer

        # plt.plot(existing_pos[0:N,0],existing_pos[0:N,1], 'b')
        #plt.plot(existing_pos[0:N,0],existing_pos[0:N,1], '.r')
        #plt.show()
    return existing_pos

def plot(beads_pos):
    import matplotlib.pyplot as plt     # plotting tools
    plt.plot(beads_pos[:,0],beads_pos[:,1], 'b')
    plt.plot(beads_pos[:,0],beads_pos[:,1], '.r')
    plt.show()

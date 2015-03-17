# Import libraries and modules
import numpy as np		                            # import numpy  
# import sys                                          # progress messages
from new_beads_positions import new_beads_pos       # calculate possible new bead positions
from calculate_energies import calculate_energies   # calculate energies for each new possible bead position
from determine_new_bead import determine_new_bead   # function used to determine the final bead position by comparing the boltzmann factors
import list_tracking as track

def user_input():
    sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    T = input('Temperature, expressed in epsilon (default: 1): ') or 0.01
    number_of_beads = input('Amount of beads per polymer: ') or 150
    plot_data = input('Plot data? (y/n, default: y): ') or 'y'
    return float(sigma), float(epsilon), float(T), int(number_of_beads), plot_data

def start(number_of_beads,sigma,epsilon,T):
    ## Fixed parameters
    angle_dof = 36                              # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    
    excisting_bead_pos = np.zeros((number_of_beads,2),dtype=float)  # initialize all bead positions
    candidate_bead_pos = np.zeros((len(angles),2),dtype=float)      # initialize list for all possible positions of the next bead

    sigma_squared = np.square(sigma)
    
    track.init(int(np.ceil(number_of_beads/2))) # initialize tracking
    
    for N in range(0, number_of_beads):
        candidate_bead_pos = new_beads_pos(excisting_bead_pos[N-1,:],angles)    # calculate all possible nodal points
        energies = calculate_energies(candidate_bead_pos,excisting_bead_pos[0:N],epsilon,sigma_squared) # calculate energies
        new_bead_index = determine_new_bead(energies,T)                         # determine final new bead
        excisting_bead_pos[N,:] = candidate_bead_pos[new_bead_index,:]          # add new final new bead to the polymer
        track.store(excisting_bead_pos[N,:],N)
    return excisting_bead_pos

def plot(beads_pos):
    import matplotlib.pyplot as plt 			           # plotting tools
    
    plt.plot(beads_pos[:,0],beads_pos[:,1], 'b')
    plt.plot(beads_pos[:,0],beads_pos[:,1], '.r')
    #plt.axis([-1*(number_of_beads), number_of_beads, -1*(number_of_beads), number_of_beads])
    
    plt.show()
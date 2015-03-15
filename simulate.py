# Import libraries and modules
import numpy as np		                               # import numpy
from new_beads_positions import new_beads_pos          # calculate possible new bead positions
from calculate_energies import calculate_energies      # calculate energies for each new possible bead position
from determine_new_bead import determine_new_bead      # function used to determine the final bead position by comparing the boltzmann factors

def start(number_of_beads,sigma,epsilon,T):
    ## Fixed parameters
    angle_dof = 18                              # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    
    beads_pos = np.zeros((number_of_beads,2),dtype=float)  # initialize all bead positions
    possible_beads_pos = np.zeros((len(angles),2),dtype=float)   # initialize list for all possible positions of the next bead
    
    for N in range(0, number_of_beads-1):
        possible_beads_pos = new_beads_pos(beads_pos[N,:],angles)  # calculate all possible nodal points
        energies = calculate_energies(possible_beads_pos,beads_pos[0:N],epsilon,sigma)
        #beads_pos[1,0]
        #print(energies[0])
        #print(energies[1])
        #print(energies[2])
        #print(energies[3])
        new_bead_index = determine_new_bead(energies,T)            # determine final new bead
        beads_pos[N+1,:] = possible_beads_pos[new_bead_index,:]    # add new final new bead to the polymer
        #plot_beads_pos = np.zeros((N+1,2),dtype=float)                             # this block is used to plot the polymer as it grows, only for tesing purposes
        #plot_beads_pos = beads_pos[0:(N+1)]
    return beads_pos
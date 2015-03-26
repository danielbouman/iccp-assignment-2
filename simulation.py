# Import libraries and modules
import numpy as np		                            # import numpy  
# import sys                                          # progress messages
import new_bead       # calculate possible new bead positions
import lj_energy

def user_input():
    # sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    # epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    # T = input('Temperature, expressed in epsilon (default: 1): ') or 1.0
    # number_of_beads = input('Amount of beads per polymer: ') or 150
    # plot_data = input('Plot data? (y/n, default: y): ') or 'y'
    sigma = 0.8
    epsilon = 0.25
    T = 1
    number_of_beads = 150
    plot_data = 'n'
    return float(sigma), float(epsilon), float(T), int(number_of_beads), plot_data

def start(number_of_beads,sigma,epsilon,T):
    ## Fixed parameters
    angle_dof = 6                                 # Amount of different angles the polymer can move in
    N_candidates = angle_dof**2
    angles1 = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    angles2 = np.linspace(0,2*np.pi,angle_dof)
    
    sigma_squared = sigma*sigma
    cutoff_length = int(2*np.ceil(1+2.5*sigma))
    existing_pos = np.zeros((number_of_beads,3),dtype=float)  # initialize all bead positions
    candidate_pos = np.zeros((len(angles1)*len(angles2),3),dtype=float)   # initialize list for all possible positions of the next bead
    angle_last_bead = [0,0]

    for N in range(1, number_of_beads):
        candidate_pos,angles1_updated,angles2_updated = new_bead.positions(existing_pos[N-1,:],angles1,angles2).reshape(-1,3)    # calculate candidate nodal points
        energies = lj_energy.func(existing_pos[0:N,:],candidate_pos,sigma_squared,epsilon,N_candidates,N,angles1_updated,angles2_updated)
        new_bead_index = new_bead.roulette(energies,T)            # determine final new bead

        angle_last_bead = [angles1_updated[ii], angles2_updated[ii]]
        existing_pos[N,:] = candidate_pos[new_bead_index,:]    # add new final new bead to the polymer
    return existing_pos

def plot(existing_pos):
    import matplotlib.pyplot as plt 			           # plotting tools
    import pylab
    from mpl_toolkits.mplot3d import Axes3D
    fig = pylab.figure()
    ax = Axes3D(fig)
    plt.plot(beads_pos[:,0], beads_pos[:,1],beads_pos[:,2])
    plt.show()
    return
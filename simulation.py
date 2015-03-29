# Import libraries and modules
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module

def user_input():
    # sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    # epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    # T = input('Temperature, expressed in epsilon (default: 1): ') or 1.0
    # number_of_beads = input('Amount of beads per polymer: ') or 150
    # plot_data = input('Plot data? (y/n, default: y): ') or 'y'
    sigma = 0.8
    epsilon = 0.25
    bending_energy = 0.0
    T = 0.1
    number_of_beads = 80
    plot_data = 'y'
    amount_of_polymers = 1000
    return float(sigma), float(epsilon), float(T), int(number_of_beads), plot_data, float(bending_energy), int(amount_of_polymers)

def start(number_of_beads,sigma,epsilon,T,bending_energy):
    ## Fixed parameters
    angle_dof = 72                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angle_dof)   # Split 2*pi radians up into angle_dof amount of slices
    
    sigma_squared = sigma*sigma
    existing_pos = np.zeros((number_of_beads,2),dtype=float)    # initialize all bead positions
    candidate_pos = np.zeros((len(angles),2),dtype=float)       # initialize list for all possible positions of the next bead
    angle_last_bead = 0
    total_weight_factor = 1

    for N in range(1, number_of_beads):
        candidate_pos,angles_updated = new_bead.positions(existing_pos[N-1,:],angles)  # calculate all possible nodal points
        energies = lj_energy.func(existing_pos[0:N,:],candidate_pos,sigma_squared,epsilon,bending_energy,angle_last_bead,angles_updated,angle_dof,N) # calculate energies
        new_bead_index,weight_factor = new_bead.roulette(energies,T)         # determine final new bead
        total_weight_factor = weight_factor*total_weight_factor
        existing_pos[N,:] = candidate_pos[new_bead_index,:]    # add new final new bead to the polymer
        angle_last_bead = angles_updated[new_bead_index]
        
    end_to_end_distance_squared = sum(np.square(existing_pos[0,:]-existing_pos[-1,:]))
    centre_of_mass = sum(existing_pos)/number_of_beads
    radius_of_gyration_squared = sum(sum(np.square(existing_pos[:,:]-centre_of_mass)))
    return existing_pos, total_weight_factor, end_to_end_distance_squared,radius_of_gyration_squared;

def calculate_expectation_value(weight_factors,quantities):
    #weight_factors[np.isnan(weight_factors)] = 0        # replace the weights that are too low for calculations by zero
    expectation_value_quantity = np.nansum(np.multiply(weight_factors,quantities))/np.nansum(weight_factors)
    return expectation_value_quantity;


def plot(beads_pos,end_to_end_distance):
    import matplotlib.pyplot as plt     # plotting tools
    plt.plot(beads_pos[:,0],beads_pos[:,1], 'b')
    plt.plot(beads_pos[:,0],beads_pos[:,1], '.r')
    #plt.plot(end_to_end_distance)
    plt.show()
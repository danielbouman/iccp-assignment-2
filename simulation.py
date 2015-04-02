# Import libraries and modulesa
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module
from datetime import datetime   # timer functions

def user_input():
    # sigma = input('Sigma value of L-J potential (default: 0.5): ') or 0.8
    # epsilon = input('Epsilon value of L-J potential (default: 0.5): ') or 0.25
    # T = input('Temperature, expressed in epsilon (default: 1): ') or 1.0
    # number_of_beads = input('Amount of beads per polymer: ') or 150
    # plotData = input('Plot data? (y/n, default: y): ') or 'y'
    #multi = input('Simulate more than one ensemble? [Y/n]  ') or 'n'
    multi = 'n'
    if multi.lower() == 'y':
        minBeads = input('Mininum number of beads: ') or 3
        maxBeads = input('Maximum number of beads: ') or 4
    else:
        minBeads = input('Number of beads: ') or 150
        maxBeads = minBeads        
    global nPolymers
    global sigma
    global epsilon
    global bendingEnergy
    global T
    global plotData

    sigma = 0.8
    epsilon = 0.25
    bendingEnergy = 0.0
    T = 0.1
    plotData = 'n'
    nPolymers = 100
    return float(sigma), float(epsilon), float(T), int(minBeads), int(maxBeads)+1, plotData, float(bendingEnergy), int(nPolymers)

def simulation(nBeads,multi,write_mode):    
    start_time = datetime.now()
    weight_factors = np.zeros((nPolymers),dtype=float)    # initialize all end_to_end distances
    end_to_end_distance_squared = np.zeros((nPolymers),dtype=float)    # initialize all end_to_end distances, squared
    radius_of_gyration_squared = np.zeros((nPolymers),dtype=float)    # initialize all end_to_end distances, squared

    # Fixed parameters
    angleDOF = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices
    
    sigmaSquared = sigma*sigma
    existingPos = np.zeros((nBeads,2),dtype=float)    # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)       # initialize list for all possible positions of the next bead
    angleLastBead = 0
    total_weight_factor = 1

    # Simulate polymers
    for ii in range(0, nPolymers):
        # beads_pos,weight_factors[ii],end_to_end_distance_squared[ii],radius_of_gyration_squared[ii] = start(nBeads,sigma,epsilon,T,bendingEnergy)            # Start simulation
        for N in range(1, nBeads):
            candidatePos,angles_updated = new_bead.positions(existingPos[N-1,:],angles)  # calculate all possible nodal points
            energies = lj_energy.func(existingPos[0:N,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,angles_updated,angleDOF,N) # calculate energies
            new_bead_index,weight_factor = new_bead.roulette(energies,T)         # determine final new bead
            weight_factors[ii] = weight_factor*total_weight_factor
            # upperWeight = 






            existingPos[N,:] = candidatePos[new_bead_index,:]    # add new final new bead to the polymer
            angleLastBead = angles_updated[new_bead_index]
        end_to_end_distance_squared[ii] = sum(np.square(existingPos[0,:]-existingPos[-1,:]))
        centre_of_mass = sum(existingPos)/nBeads
        radius_of_gyration_squared[ii] = sum(sum(np.square(existingPos[:,:]-centre_of_mass)))

    print(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    
    # Collect properties
    end_to_end_distance = np.sqrt(end_to_end_distance_squared)
    radius_of_gyration = np.sqrt(radius_of_gyration_squared)
    exp_end_to_end_distance = calculate_expectation_value(weight_factors,end_to_end_distance)
    exp_end_to_end_distance_squared = calculate_expectation_value(weight_factors,end_to_end_distance_squared)
    exp_radius_of_gyration = calculate_expectation_value(weight_factors,radius_of_gyration)
    
    if multi:
        # Save data to file
        save.save(exp_end_to_end_distance_squared,"R_squared",header="",write_mode=write_mode)
        save.save(exp_radius_of_gyration,"exp_radius_of_gyration",header="",write_mode=write_mode)
    else:
        print("End to end distance: " + str(exp_end_to_end_distance))
        print("Radius of gyration: " + str(exp_radius_of_gyration))

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
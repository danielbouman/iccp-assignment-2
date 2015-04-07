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
        #minBeads = input('Number of beads: ') or 150
        minBeads = 150
        maxBeads = minBeads

    global nBeads       
    global nPolymers
    global sigma
    global epsilon
    global bendingEnergy
    global T
    global plotData
    global upLimInit
    global lowLimInit

    upLimInit = 1.2
    lowLimInit = 1.0
    sigma = 0.8
    epsilon = 0.25
    bendingEnergy = 0.0
    T = 1
    plotData = 'n'
    nPolymers = 5
    nBeads = 35
    return float(sigma), float(epsilon), float(T), int(minBeads), int(maxBeads)+1, plotData, float(bendingEnergy), int(nPolymers)

def addBead(L,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim):
    candidatePos,angles_updated = new_bead.positions(existingPos[L-1,:],angles)  # calculate all possible nodal points
    energies = lj_energy.func(existingPos[0:L,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,angles_updated,angleDOF,L) # calculate energies
    new_bead_index,weight_factor = new_bead.roulette(energies,T)         # determine final new bead
    total_weight_factors = weight_factor*total_weight_factors
    print("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(total_weight_factors)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim))
    existingPos[L,:] = candidatePos[new_bead_index,:]    # add new final new bead to the polymer
    angleLastBead = angles_updated[new_bead_index]

    
    upLim = 1.1*upLim
    lowLim = 0.9*lowLim

    if L < nBeads-1:
        if total_weight_factors>upLim:
            print('upLim')
            # Clone current polymer
            addBead(L+1,N+1,existingPos,candidatePos,angles,angleLastBead,0.5*total_weight_factors,upLim,lowLim)
            # Multiply weight by 0.5 and continue growing current polymer
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,0.5*total_weight_factors,upLim,lowLim)
            
            # nPolymers = nPolymers + 1
            
        elif total_weight_factors<lowLim:
            print('lowLim')
            RNG = np.random.random()
            if RNG<0.5:
                print('Polymer weight doubled')
                return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,2*total_weight_factors,upLim,lowLim)
            else:
                print('Polymer removed')
                return addBead(1,N,existingPos,candidatePos,angles,angleLastBead,1,upLimInit,lowLimInit)
        else:
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim)

    N = N + 1
    return existingPos,total_weight_factors,N


def simulation(nBeads,multi,write_mode):    
    start_time = datetime.now()
    weight_factors = np.ones((nPolymers),dtype=float)    # initialize all end_to_end distances
    end_to_end_distance_squared = np.zeros((nPolymers),dtype=float)    # initialize all end_to_end distances, squared
    radius_of_gyration_squared = np.zeros((nPolymers),dtype=float)    # initialize all end_to_end distances, squared
    

    # Fixed parameters
    global sigmaSquared
    global angleDOF
    global angles
    angleDOF = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices

    sigmaSquared = sigma*sigma
    existingPos = np.zeros((nBeads,2),dtype=float)    # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)       # initialize list for all possible positions of the next bead
    angleLastBead = 0
    total_weight_factors = 1

    N = 0
    L = 1
    
    # Simulate polymers
    while N < nPolymers:
        
        print("nPolymers: "+str(nPolymers))
        # Create polymer
        total_weight_factors = 1
        upLim = upLimInit
        lowLim = lowLimInit
        L = 1
        existingPos,weight_factor_single,N = addBead(L,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim)
        
        end_to_end_distance_squared[N-1] = sum(np.square(existingPos[0,:]-existingPos[-1,:]))
        centre_of_mass = sum(existingPos)/nBeads
        radius_of_gyration_squared[N-1] = sum(sum(np.square(existingPos[:,:]-centre_of_mass)))
        print("nPolymers: "+str(nPolymers))

    print(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    
    # Collect properties
    end_to_end_distance = np.sqrt(end_to_end_distance_squared)
    radius_of_gyration = np.sqrt(radius_of_gyration_squared)
    exp_end_to_end_distance = calculate_expectation_value(weight_factors,end_to_end_distance)
    exp_end_to_end_distance_squared = calculate_expectation_value(weight_factors,end_to_end_distance_squared)
    exp_radius_of_gyration = calculate_expectation_value(weight_factors,radius_of_gyration)
    # print(weight_factors)
    # print(exp_end_to_end_distance)
    if multi:
        # Save data to file
        save.save(exp_end_to_end_distance_squared,"R_squared",header="",write_mode=write_mode)
        save.save(exp_radius_of_gyration,"exp_radius_of_gyration",header="",write_mode=write_mode)
    else:
        print("End to end distance: " + str(exp_end_to_end_distance))
        print("Radius of gyration: " + str(exp_radius_of_gyration))
    return existingPos;


def calculate_expectation_value(weight_factors,quantities):
    #weight_factors[np.isnan(weight_factors)] = 0        # replace the weights that are too low for calculations by zero
    expectation_value_quantity = sum(np.multiply(weight_factors,quantities))/sum(weight_factors)
    return expectation_value_quantity;


def plot(existingPos):
    import matplotlib.pyplot as plt     # plotting tools
    plt.plot(existingPos[:,0],existingPos[:,1], 'b')
    plt.plot(existingPos[:,0],existingPos[:,1], '.r')
    #plt.plot(end_to_end_distance)
    plt.show()
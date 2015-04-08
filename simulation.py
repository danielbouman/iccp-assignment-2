# Import libraries and modulesa
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module
from datetime import datetime   # timer functions

def user_input():
    # Call global variables
    global nBeads       
    global nPolymers
    global sigma
    global epsilon
    global bendingEnergy
    global T
    global plotData
    global upLimInit
    global lowLimInit
    
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

    upLimInit = 1.2
    lowLimInit = 1.0
    sigma = 0.8
    epsilon = 0.25
    bendingEnergy = 0.0
    T = 1
    plotData = 'n'
    nPolymers = 300
    nBeads = 150
    return float(sigma), float(epsilon), float(T), int(minBeads), int(maxBeads)+1, plotData, float(bendingEnergy), int(nPolymers)

def addBead(L,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim):
    # Determine candidate positions for new bead
    candidatePos,angles_updated = new_bead.positions(existingPos[L-1,:],angles)  
    # Calculate total energy for each candidate position
    energies = lj_energy.func(existingPos[0:L,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,angles_updated,angleDOF,L)
    new_bead_index,weight_factor = new_bead.roulette(energies,T)         # determine final new bead
    total_weight_factors = weight_factor*total_weight_factors
    existingPos[L-1,:] = candidatePos[new_bead_index,:]    # add new final new bead to the polymer
    angleLastBead = angles_updated[new_bead_index]
    
    
    print("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(total_weight_factors)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim))

    upLim = 1.1*upLim
    lowLim = 0.11*lowLim

    if L < nBeads:
        global nPolymers
        if total_weight_factors>upLim:
            # Enrichment
            nPolymers = nPolymers + 1
            print('upLim reached, clone current polymer. nPolymers: '+str(nPolymers))
            # Clone current polymer and grow
            addBead(L+1,nPolymers,existingPos,candidatePos,angles,angleLastBead,0.5*total_weight_factors,upLim,lowLim)
            print('Done with cloned polymer, continue growing original')
            # When done with cloned polymer, multiply weight original polymer by 0.5 and continue growing
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,0.5*total_weight_factors,upLim,lowLim)
            
        if total_weight_factors<lowLim:
            # Pruning
            RNG = np.random.random()
            if RNG<0.5:
                print('lowLim reached, polymer weight doubled')
                # Double current polymer weight and continue growing
                return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,2*total_weight_factors,upLim,lowLim)
            else:
                # Remove current bead
                nPolymers = nPolymers - 1
                print('lowLim reached, polymer removed. nPolymers: '+str(nPolymers))
                return existingPos,2,N
        else:
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim)
            
    N = N + 1
    return existingPos,total_weight_factors,N


def simulation(nBeads,multi,write_mode):    
    start_time = datetime.now()
    weight_factors = np.ones((nPolymers*2),dtype=float)                   # initialize all end_to_end distances
    end_to_end_distance_squared = np.zeros((nPolymers*2),dtype=float)     # initialize all end_to_end distances, squared
    radius_of_gyration_squared = np.zeros((nPolymers*2),dtype=float)      # initialize all end_to_end distances, squared
    
    # Fixed parameters
    global sigmaSquared
    global angleDOF
    global angles
    angleDOF = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices

    sigmaSquared = sigma*sigma
    existingPos = np.zeros((nBeads,2),dtype=float)          # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)    # initialize list for all possible positions of the next bead
    angleLastBead = 0
    total_weight_factors = 1

    N = 0
    L = 1
    
    # Simulate polymers
    while N < nPolymers:
        # Create polymer
        total_weight_factors = 1    # Starting polymer weight factor
        upLim = upLimInit           # When reached enrich
        lowLim = lowLimInit         # When reached prune
        L = 1                       # Starting bead
        # Grow beads
        existingPos,weight_factor_single,N = addBead(L,N,existingPos,candidatePos,angles,angleLastBead,total_weight_factors,upLim,lowLim)
        # Collect data
        end_to_end_distance_squared[N-1] = sum(np.square(existingPos[0,:]-existingPos[-1,:]))
        centre_of_mass = sum(existingPos)/nBeads
        radius_of_gyration_squared[N-1] = sum(sum(np.square(existingPos[:,:]-centre_of_mass)))
        print("nPolymers: "+str(nPolymers))

    print(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    
    # Collect data
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
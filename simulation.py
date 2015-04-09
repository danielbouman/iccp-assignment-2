# Import libraries and modulesa
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module
from datetime import datetime   # timer functions

def user_input():
    # Get global variables
    global nBeads       
    global nPolymers
    global sigma
    global epsilon
    global bendingEnergy
    global T
    global plotData
    global upLimInit
    global lowLimInit
    global startPolymers
    
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
        minBeads = 30
        maxBeads = minBeads

    upLimInit = 1.2     # Initial upper polymerWeight limit for enrichment
    lowLimInit = 1    # Initial lower polymerWeight limit for pruning
    sigma = 0.8         # Sigma for the Lennard Jones potential
    epsilon = 0.25      # epsilon for the Lennard Jones potential
    bendingEnergy = 0.0
    T = 1               # Temperature
    plotData = 'n'      # Plot data boolean
    startPolymers = 30  # Ensemble size
    nBeads = 30         # Numer of beads per polymer
    return int(minBeads), int(maxBeads)+1, plotData
    
def partitionFunction(L,oldZ):
    completedPolymers
    currentZ = oldZ(L)
    return currentZ

def addBead(L,N,existingPos,candidatePos,baseAngles,angleLastBead,currentPolymerWeight,upLim,lowLim):
    # Variables
    global diagFile
    
    # Determine candidate positions for new bead
    candidatePos,candidateAngles = new_bead.positions(existingPos[L-1,:],baseAngles)  
    # Calculate energies for each candidate position
    energies = lj_energy.func(existingPos[0:L,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,candidateAngles,angleDOF,L)
    # Choose bead position from candidates and determine bead weight
    chosenBeadIndex,beadWeight = new_bead.roulette(energies,T)
    # Add chosen bead position to existing positions
    existingPos[L-1,:] = candidatePos[chosenBeadIndex,:]
    # Update current polymer weight
    currentPolymerWeight = beadWeight*currentPolymerWeight
    # Save chosen angle of current bead for bending energy with next bead
    angleLastBead = candidateAngles[chosenBeadIndex]
    # Update upper polymerWeight limit for enrichment and lower polymerWeight limit for pruning
    # partitionFunction = 
    # Write to diagnostics file
    diagFile.write("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(currentPolymerWeight)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim)+"\n")
    upLim = 1.11*upLim
    lowLim = 0.11*lowLim

    # Grow new bead if final length has not been reached
    if L < nBeads:
        global nPolymers
        # Enrichment
        if currentPolymerWeight>upLim:
            # Increase population with one polymer
            diagFile.write('upLim reached, clone current polymer. nPolymers: '+str(nPolymers)+'\n')                 # Write to diagnostics file
            global beadPos
            global polymerWeights
            global polymerAdded
            global nPolymersEnrich
            nPolymersEnrich = nPolymersEnrich + 1
            beadPos, polymerWeights[nPolymersEnrich], Ntemp = addBead(L+1,nPolymersEnrich,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight,upLim,lowLim) # Clone current polymer and grow clone 
            diagFile.write('Done with cloned polymer, continue growing original polymer ('+str(N)+')\n')                              # Write to diagnostics file
            if beadPos[-1,0] != 0:
                diagFile.write('Cloned polymer '+str(nPolymersEnrich)+' was fully grown'+'\n')
                print('Cloned polymer '+str(nPolymersEnrich)+' was fully grown'+'\n')
            else:
                diagFile.write('Cloned polymer was removed'+'\n')
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight,upLim,lowLim)      # When finished growing cloned polymer, multiply weight original polymer by 0.5 and continue growing
            
        # Pruning
        if currentPolymerWeight<lowLim:
            RNG = np.random.random()
            if RNG<0.5:
                # Double current polymer weight and continue growing
                diagFile.write('lowLim reached, polymer weight doubled\n')
                return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,2*currentPolymerWeight,upLim,lowLim)
            else:
                # Remove current bead
                nPolymers = nPolymers - 1
                diagFile.write('lowLim reached, polymer removed. nPolymers: '+str(nPolymers)+'\n')
                return existingPos, currentPolymerWeight, N
        else:
            # Continue with next bead
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,currentPolymerWeight,upLim,lowLim)
    N = N + 1  # Select the next polymer
    return existingPos, currentPolymerWeight, N


def simulation(nBeads,multi,write_mode):    
    start_time = datetime.now()
    
    # Variables
    global sigmaSquared
    global angleDOF
    global angles
    global diagFile
    global nPolymersEnrich
    global startPolymers
    global nPolymers
    global beadPos
    global polymerWeights
    nPolymersEnrich = startPolymers
    angleDOF = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices
    sigmaSquared = sigma*sigma
    beadPos = np.zeros((nBeads,2),dtype=float)          # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)    # initialize list for all possible positions of the next bead
    angleLastBead = 0
    total_weight_factors = 1
    nPolymers = startPolymers
    polymerWeights = np.ones((startPolymers*2),dtype=float)                   # initialize all end_to_end distances
    end_to_end_distance_squared = np.zeros((startPolymers*2),dtype=float)     # initialize all end_to_end distances, squared
    radius_of_gyration_squared = np.zeros((startPolymers*2),dtype=float)      # initialize all end_to_end distances, squared
    N = 0   # First polymer
    
    # Open file for diagnostic data
    diagFile = open("diag.txt","w")
    
    # Simulate polymers
    while N < nPolymers:
        
        # Variables
        polymerWeightInit = 1   # Initial polymer weight factor
        upLim = upLimInit       # When reached enrich
        lowLim = lowLimInit     # When reached prune
        L = 1                   # Starting bead
        
        # Grow beads recursively
        beadPos, polymerWeights[N], N = addBead(L,N,beadPos,candidatePos,angles,angleLastBead,polymerWeightInit,upLim,lowLim)
        
        # Collect data
        end_to_end_distance_squared[N-1] = sum(np.square(beadPos[0,:]-beadPos[-1,:]))
        centre_of_mass = sum(beadPos)/nBeads
        radius_of_gyration_squared[N-1] = sum(sum(np.square(beadPos[:,:]-centre_of_mass)))

    # Write diagnostics
    diagFile.write(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    diagFile.close()
    
    print(str(nPolymers) + " polymers with " + str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    
    # Collect data
    #print(polymerWeights)
    #polymerWeights = polymerWeights[0:nPolymers]
    #print(polymerWeights)
    # polymerWeights = [0 if x==1 else x for x in polymerWeights]

    print(polymerWeights)

    end_to_end_distance = np.sqrt(end_to_end_distance_squared)
    radius_of_gyration = np.sqrt(radius_of_gyration_squared)
    exp_end_to_end_distance = calculate_expectation_value(polymerWeights,end_to_end_distance)
    exp_end_to_end_distance_squared = calculate_expectation_value(polymerWeights,end_to_end_distance_squared)
    exp_radius_of_gyration = calculate_expectation_value(polymerWeights,radius_of_gyration)
    if multi:
        # Save data to file
        save.save(exp_end_to_end_distance_squared,"R_squared",header="",write_mode=write_mode)
        save.save(exp_radius_of_gyration,"exp_radius_of_gyration",header="",write_mode=write_mode)
    else:
        print("End to end distance: " + str(exp_end_to_end_distance))
        print("Radius of gyration: " + str(exp_radius_of_gyration))
    return beadPos;


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
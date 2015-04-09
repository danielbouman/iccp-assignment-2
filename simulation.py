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
    sigma = 0.8         # Sigma for the Lennard Jones potential
    epsilon = 0.25      # epsilon for the Lennard Jones potential 
    bendingEnergy = 0.0
    T = 0.5               # Temperature
    plotData = 'n'      # Plot data boolean
    startPolymers = 30  # Ensemble size
    nBeads = 150         # Numer of beads per polymer
    return int(minBeads), int(maxBeads)+1, plotData
    
# Partition function
def partitionFunction(oldZ,completedPolymers,currentPolymerWeight):
    return (oldZ*completedPolymers+currentPolymerWeight)/(completedPolymers+1)
    
# Add bead function
def addBead(L,N,existingPos,candidatePos,baseAngles,angleLastBead,currentPolymerWeight):
    # Declare global variables
    global diagFile
    global completedPolymers
    global Z
    global oldZ
    
    candidatePos,candidateAngles = new_bead.positions(existingPos[L-1,:],baseAngles)    # Determine candidate positions for new bead)
    energies = lj_energy.func(existingPos[0:L-1,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,candidateAngles,angleDOF,L) # Calculate energies for each candidate position
    chosenBeadIndex,beadWeight = new_bead.roulette(energies,T) # Choose bead position from candidates and determine bead weight
    existingPos[L,:] = candidatePos[chosenBeadIndex,:] # Add chosen bead position to existing positions
    currentPolymerWeight = beadWeight*currentPolymerWeight # Update current polymer weight
    angleLastBead = candidateAngles[chosenBeadIndex] # Save chosen angle of current bead for bending energy with next bead 
    
    Z[L] = partitionFunction(oldZ[L],completedPolymers,currentPolymerWeight)
    
    if N != 0:
        upLim = 2*Z[L]
        lowLim = 1.2*Z[L]
    else:
        upLim = float('inf')
        lowLim = 0
        
    # Write to diagnostics information to file
    diagFile.write("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(currentPolymerWeight)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim)+"\n")
    
    # Grow new bead if final length has not been reached
    if L < nBeads-1:
        global nPolymers
        global polymerPruned
        # print(currentPolymerWeight)
        # plt.plot(existingPos[:,0],existingPos[:,1], 'b')
        # plt.show
        if currentPolymerWeight>upLim:  # Enrich
            # Declare global variables
            global beadPos
            global polymerWeights
            global polymersAdded
            # global nPolymersEnrich
            diagFile.write('upLim reached, clone current polymer. nPolymers: '+str(nPolymers)+'\n') # Write diagnostics to file
            # nPolymersEnrich = nPolymersEnrich + 1 # Increase population
            beadPos, polymerWeights[nPolymersEnrich], Ntemp = addBead(L+1,N+1,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight) # Clone current polymer and grow clone 
            diagFile.write('Done with cloned polymer, continue growing original polymer ('+str(N)+')\n') # Write diagnostics to file
            # Write diagnostics to file
            diagFile.write(str(beadPos))
            if polymerPruned:
                diagFile.write('Cloned polymer was removed'+'\n')
                polymerPruned = False
            else:
                diagFile.write('Cloned polymer '+str(Ntemp)+' was fully grown'+'\n')
                polymersAdded = polymersAdded + 1
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight) # When finished growing cloned polymer, multiply weight original polymer by 0.5 and continue growing
            
        if currentPolymerWeight<lowLim: # Prune
            RNG = np.random.random()
            if RNG<0.5: # Double current polymer weight and continue growing
                diagFile.write('lowLim reached, polymer weight doubled\n')
                return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,2*currentPolymerWeight) # Write diagnostics to file
            else:   # Remove current bead
                nPolymers = nPolymers - 1 # Decrease population
                polymerPruned = True
                diagFile.write('lowLim reached, polymer removed. nPolymers: '+str(nPolymers)+'\n') # Write diagnostics to file
                return existingPos, 0, N
        else:
            # Continue with next bead
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,currentPolymerWeight)
    completedPolymers = completedPolymers + 1   # Count completed polymers
    oldZ = Z
    diagFile.write('Polymers completed:'+str(completedPolymers)+'\n')
    N = N + 1 + polymersAdded                                 # Select the next polymer
    return existingPos, currentPolymerWeight, N


def simulation(multi,write_mode):    
    start_time = datetime.now()
    # Declare global variables
    global sigmaSquared
    global angleDOF
    global angles
    global diagFile
    global nPolymersEnrich
    global startPolymers
    global nPolymers
    global beadPos
    global polymerWeights
    global nBeads
    global completedPolymers
    global oldZ
    global Z
    global polymersAdded
    polymersAdded = 0
    completedPolymers = 0
    nPolymersEnrich = startPolymers
    angleDOF = 6                               # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices
    sigmaSquared = sigma*sigma
    beadPos = np.zeros((nBeads,2),dtype=float)          # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)    # initialize list for all possible positions of the next bead
    angleLastBead = 0
    nPolymers = startPolymers
    polymerWeights = np.ones((startPolymers*3),dtype=float)                   # initialize all end_to_end distances
    end_to_end_distance_squared = np.zeros((startPolymers*3),dtype=float)     # initialize all end_to_end distances, squared
    radius_of_gyration_squared = np.zeros((startPolymers*3),dtype=float)      # initialize all end_to_end distances, squared
    Z = np.zeros((nBeads),dtype=float)
    oldZ = np.zeros((nBeads),dtype=float)
    N = 0   # First polymer
    diagFile = open("diag.txt","w") # Open file for diagnostic data
    
    while N < nPolymers: # Simulate polymers
        # Variables
        polymerWeightInit = 1   # Initial polymer weight factor
        L = 2                   # Starting bead
        beadPos[1,:] = [1,0]
        # Grow beads recursively
        beadPos, polymerWeights[N], N = addBead(L,N,beadPos,candidatePos,angles,angleLastBead,polymerWeightInit)
        # Collect data
        end_to_end_distance_squared[N-1] = sum(np.square(beadPos[0,:]-beadPos[-1,:]))
        centre_of_mass = sum(beadPos)/nBeads
        radius_of_gyration_squared[N-1] = sum(sum(np.square(beadPos[:,:]-centre_of_mass)))

    # Write diagnostics
    diagFile.write(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    diagFile.close()
    
    print(str(completedPolymers) + " polymers with " + str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    # Collect data
    #print(polymerWeights)
    #polymerWeights = polymerWeights[0:nPolymers]
    #print(polymerWeights)
    polymerWeights = [0 if x==1 else x for x in polymerWeights]

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
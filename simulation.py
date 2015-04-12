# Import libraries and modulesa
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module
from datetime import datetime   # timer functions
import sys
import save_data as data

"""
This function defines all the global parameters for the simulation.
When multiple ensembles are run these are also defined here.
"""
def variables():
    # Get global variables
    global nBeads       
    global nPolymers
    global sigma
    global epsilon
    global bendingEnergy
    global T
    global plotData
    global startPolymers
    global multiplePolymerLengths
    global multipleTemperatures
    global multipleBendingenergies
    
    sigma = 0.8             # Siagma for the Lennard Jones potential
    epsilon = 0.25          # epsilon for the Lennard Jones potential 
    bendingEnergy = 0.0
    plotData = False        # Plot data boolean
    startPolymers = 10000   # Ensemble size

    # Single ensemble
    nBeads = 30
    T = 1

    # Here multiple ensembles with varying polymer length, temperature or bending energy are defined
    multiplePolymerLengths = False
    minBeads = 148
    maxBeads = 150
    
    multipleTemperatures = False
    minT = 0.05
    maxT = 10
    stepT = 0.1

    multipleBendingenergies = False
    minBending = 0.20
    maxBending = 1.5
    stepBending = 0.05
    
    if multiplePolymerLengths:
        multipleTemperatures = False
        multipleBendingenergies = False
        return int(minBeads), int(maxBeads+1), T, T+1, stepT, bendingEnergy, bendingEnergy+1, 1
    elif multipleTemperatures:
        return int(nBeads), int(nBeads+1), minT, maxT+1, stepT, bendingEnergy, bendingEnergy+1, 1
    elif multipleBendingenergies:
        return int(nBeads), int(nBeads+1), T, T+1, 1, minBending, maxBending, stepBending
    else:
        return int(nBeads), int(nBeads+1), T, T+1, 1, bendingEnergy, bendingEnergy+1, 1
    
"""
Recursive function for growing beads, this is the main part of the simulation.
Main features: energy calculations, weights calculations, pruning and enrichment.
"""
def addBead(L,N,existingPos,candidatePos,baseAngles,angleLastBead,currentPolymerWeight):
    # Declare global variables
    # global diagFile
    global completedPolymers
    global Z
    global oldZ
    global end_to_end_distance_squared
    global radius_of_gyration_squared
    global nBeads
    # Here the new bead position, bead weight and angle with the previous bead are determined
    candidatePos,candidateAngles = new_bead.positions(existingPos[L-1,:],baseAngles)    # Determine candidate positions for new bead)
    energies = lj_energy.func(existingPos[0:L-1,:],candidatePos,sigmaSquared,epsilon,bendingEnergy,angleLastBead,candidateAngles,angleDOF,L) # Calculate energies for each candidate position
    chosenBeadIndex,beadWeight = new_bead.roulette(energies,T,L) # Choose bead position from candidates and determine bead weight
    existingPos[L,:] = candidatePos[chosenBeadIndex,:] # Add chosen bead position to existing positions
    currentPolymerWeight = beadWeight*currentPolymerWeight # Update current polymer weight
    angleLastBead = candidateAngles[chosenBeadIndex] # Save chosen angle of current bead for bending energy with next bead 
    Z[L] = (oldZ[L]*completedPolymers+currentPolymerWeight)/(completedPolymers+1) # Partition function
    # This makes sure that the first polymer is grown using RR algorithm and subsequent polymers with PERM
    if N != 0:
        upLim = 2.43*oldZ[L]/oldZ[2]
        lowLim = 1.46*oldZ[L]/oldZ[2]
    else:
        if currentPolymerWeight == 0:
            sys.exit("First polymer (RR algorithm) has weight 0, restart simulation.")
        upLim = float('inf')
        lowLim = 0
    # Write to diagnostics information to file
    # diagFile.write("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(currentPolymerWeight)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim)+"\n")
    
    # Grow new bead if final length has not been reached
    if L < nBeads-1:
        global nPolymers
        # When the upper limit is reached enrich the population
        if currentPolymerWeight>upLim:
            sys.stdout.write("Polymer: %d of %d \r" % (N,nPolymers) )
            sys.stdout.flush()
            # Declare global variables
            global beadPos
            global polymerWeights
            # Write diagnostic inforation to file
            # diagFile.write('upLim reached, clone current polymer. nPolymers: '+str(nPolymers)+'\n')
            # diagFile.write(str(existingPos)+'\n'+str(L)+'\n')
            # Look for available locations for a new polymer, and split current polymer
            tempN = np.nonzero(polymerWeights == 1)[0][0]

            polymerWeights[tempN] = 0
            nPolymers = nPolymers + 1
            beadPos, polymerWeights[tempN], _ = addBead(L+1,tempN,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight)
            # diagFile.write('Done with cloned polymer, continue growing original polymer ('+str(N)+')\n')
            # diagFile.write(str(existingPos))
            # When done with split polymer, continue growing original polymer
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight)
        # When the lower limit is reached, prune the polymer
        if currentPolymerWeight<lowLim:
            RNG = np.random.random()
            if RNG<0.5: 
                # The current polymer weight is doubled and growing continues
                # diagFile.write('lowLim reached, polymer weight doubled\n')
                return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,2*currentPolymerWeight) # Write diagnostics to file
            else:   
                # The current bead is removed
                nPolymers = nPolymers - 1
                # diagFile.write('lowLim reached, polymer removed. nPolymers: '+str(nPolymers)+'\n') # Write diagnostics to file
                return existingPos, 1, N
        # When neither limit is reached, continue growing the polymer
        else:
            return addBead(L+1,N,existingPos,candidatePos,angles,angleLastBead,currentPolymerWeight)
    # A polymer is now fully grown
    # Partition function is updated and quantities are determined
    completedPolymers = completedPolymers + 1
    centre_of_mass = sum(existingPos)/nBeads
    radius_of_gyration_squared[N] = sum(sum(np.square(existingPos[:,:]-centre_of_mass)))/nBeads
    end_to_end_distance_squared[N] = sum(np.square(existingPos[0,:]-existingPos[-1,:]))
    oldZ[:] = Z[:]
    if plotData:
        plot(existingPos)
    # diagFile.write('Polymers completed:'+str(completedPolymers)+'\n')
    # Next polymer
    N = np.nonzero(polymerWeights == 1)[0][0]
    return existingPos, currentPolymerWeight, N

"""
Function to start the simuluation. 
From here the polymers are grown and the quantities of interest are determined and collected.
"""
def start(nBeadsVar,TVar,bendVar):    
    start_time = datetime.now()
    # Declare global variables
    global sigmaSquared
    global angleDOF
    global angles
    # global diagFile
    global startPolymers
    global nPolymers
    global beadPos
    global polymerWeights
    global nBeads
    global completedPolymers
    global oldZ
    global Z
    global end_to_end_distance_squared
    global radius_of_gyration_squared
    global multipleBendingenergies
    global T
    
    # Switch to variable properties when simulating mulitple ensembes
    if multiplePolymerLengths:
        nBeads = nBeadsVar
    elif multipleTemperatures:
        T = TVar
    elif multipleBendingenergies:
        bendingEnergy = bendVar
        
    # Value assignment 
    completedPolymers = 0
    angleDOF = 6                              # Amount of different angles the polymer can move in
    angles = np.linspace(0,2*np.pi,angleDOF)   # Split 2*pi radians up into angleDOF amount of slices
    sigmaSquared = sigma*sigma
    beadPos = np.zeros((nBeads,2),dtype=float)          # initialize all bead positions
    candidatePos = np.zeros((len(angles),2),dtype=float)    # initialize list for all possible positions of the next bead
    angleLastBead = 0
    nPolymers = startPolymers
    # Pre-allocation of arrays
    polymerWeights = np.ones((startPolymers*3),dtype=float)
    end_to_end_distance_squared = np.zeros((startPolymers*3),dtype=float)
    radius_of_gyration_squared = np.zeros((startPolymers*3),dtype=float)
    Z = np.zeros((nBeads),dtype=float)
    oldZ = np.zeros((nBeads),dtype=float)
    N = 0   # First polymer
    # diagFile = open("diag.txt","w") # Open file for diagnostic data
    
    while N < nPolymers: # Simulate polymers
        # Variables
        polymerWeightInit = 1   # Initial polymer weight factor
        L = 2                   # Starting bead
        beadPos[1,:] = [1,0]
        polymerWeights[N] = 0
        # Grow beads recursively
        beadPos, polymerWeights[N], N = addBead(L,N,beadPos,candidatePos,angles,angleLastBead,polymerWeightInit)
        sys.stdout.write("Polymer: %d of %d \r" % (N,nPolymers) )
        sys.stdout.flush()

    # Diagnost;ic information
    # diagFile.write(str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    # diagFile.close()
    print(str(completedPolymers) + " polymers with " + str(nBeads) + " beads, done in: " + str(datetime.now() - start_time))
    polymerWeights = [0 if x==1 else x for x in polymerWeights]

    # Quantities are determined
    end_to_end_distance = np.sqrt(end_to_end_distance_squared)
    radius_of_gyration = np.sqrt(radius_of_gyration_squared)
    exp_end_to_end_distance = calculate_expectation_value(polymerWeights,end_to_end_distance)
    exp_end_to_end_distance_squared = calculate_expectation_value(polymerWeights,end_to_end_distance_squared)
    exp_radius_of_gyration = calculate_expectation_value(polymerWeights,radius_of_gyration)
    exp_radius_of_gyration_squared = calculate_expectation_value(polymerWeights,radius_of_gyration_squared)
    sd_end_to_end_distance = np.sqrt(exp_end_to_end_distance_squared-np.square(exp_end_to_end_distance))
    sd_radius_of_gyration = np.sqrt(exp_radius_of_gyration_squared-np.square(exp_radius_of_gyration))
    
    # Save quantities to external files
    if multiplePolymerLengths:
        opt_data_R = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_end_to_end_distance)
        data.save(exp_end_to_end_distance_squared,"R_squared_polVar",header="",write_mode="a",optional_data=opt_data_R)
        opt_data_G = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_radius_of_gyration)
        data.save(exp_radius_of_gyration_squared,"radius_of_gyration_polVar",header="",write_mode="a",optional_data=opt_data_G)
    elif multipleTemperatures:
        opt_data_R = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_end_to_end_distance)+" "+str(T)
        data.save(exp_end_to_end_distance_squared,"R_squared_tempVar",header="",write_mode="a",optional_data=opt_data_R)
        opt_data_G = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_radius_of_gyration)+" "+str(T)
        data.save(exp_radius_of_gyration_squared,"radius_of_gyration_tempVar",header="",write_mode="a",optional_data=opt_data_G)
    elif multipleBendingenergies:
        opt_data_R = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_end_to_end_distance)+" "+str(bendingEnergy)
        data.save(exp_end_to_end_distance_squared,"R_squared_bendVar",header="",write_mode="a",optional_data=opt_data_R)
        opt_data_G = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_radius_of_gyration)+" "+str(bendingEnergy)
        data.save(exp_radius_of_gyration_squared,"radius_of_gyration_bendVar",header="",write_mode="a",optional_data=opt_data_G)
    
    return beadPos;

"""
Function to determine weighted average quantities
"""
def calculate_expectation_value(weight_factors,quantities):
    return sum(np.multiply(weight_factors,quantities))/sum(weight_factors);

"""
Function to plot bead positions
"""
def plot(existingPos):
    import matplotlib.pyplot as plt
    plt.plot(existingPos[:,0],existingPos[:,1], 'b')
    plt.plot(existingPos[:,0],existingPos[:,1], '.r')
    plt.show()

"""
Function returns a range with float steps. Input:
start   : start Value
stop    : end value
step    : step size
"""
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step
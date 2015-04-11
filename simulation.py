# Import libraries and modulesa
import numpy as np	# import numpy  
import new_bead     # determine new bead positions
import lj_energy    # fortran lj_energy module
from datetime import datetime   # timer functions
import sys
import save_data as data

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
    
    multiplePolymerLengths = True
    minBeads = 108
    maxBeads = 150 
    sigma = 0.8             # Sigma for the Lennard Jones potential
    epsilon = 0.25          # epsilon for the Lennard Jones potential 
    bendingEnergy = 0.0
    T = 350                 # Temperature
    plotData = 'n'          # Plot data boolean
    startPolymers = 10000   # Ensemble size
    nBeads = 3              # Numer of beads per polymer
    return int(minBeads), int(maxBeads)+1, plotData
    
# Add bead function
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
    # This makes sure that only the first polymer is grown with the Rosenbluth method
    if N != 0:
        upLim = 8.53*oldZ[L]/oldZ[2]
        lowLim = 1.16*oldZ[L]/oldZ[2]
    else:
        if currentPolymerWeight == 0:
            sys.exit("First polymer (RR) has weight 0")
        upLim = float('inf')
        lowLim = 0
    # Write to diagnostics information to file
    # diagFile.write("Polymer "+str(N)+", bead "+str(L)+", weight: "+ str(currentPolymerWeight)+" lowLim: "+str(lowLim)+", upLim: "+str(upLim)+"\n")
    
    # Grow new bead if final length has not been reached
    if L < nBeads-1:
        global nPolymers
        # When the upper limit is reached enrich the population
        if currentPolymerWeight>upLim:
            # sys.stdout.write("Polymer: %d of %d \r" % (N,nPolymers) )
            # sys.stdout.flush()
            # Declare global variables
            global beadPos
            global polymerWeights
            # Write diagnostic inforation to file
            # diagFile.write('upLim reached, clone current polymer. nPolymers: '+str(nPolymers)+'\n')
            # diagFile.write(str(existingPos)+'\n'+str(L)+'\n')
            # Look for free locations for a new polymer, and split current polymer
            tempN = np.nonzero(polymerWeights == 1)[0][0]

            polymerWeights[tempN] = 0
            nPolymers = nPolymers + 1
            beadPos, polymerWeights[tempN], _ = addBead(L+1,tempN,existingPos,candidatePos,angles,angleLastBead,0.5*currentPolymerWeight)
            # diagFile.write('Done with cloned polymer, continue growing original polymer ('+str(N)+')\n')
            # diagFile.write(str(existingPos))
            # When done with split polymer, continue growing original polymer
            oldZ[:] = Z[:]
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
                oldZ[:] = Z[:]
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
    # diagFile.write('Polymers completed:'+str(completedPolymers)+'\n')
    # Next polymer
    N = np.nonzero(polymerWeights == 1)[0][0]
    return existingPos, currentPolymerWeight, N


def start(nBeadsVar):    
    start_time = datetime.now()
    # Declare global variables
    global sigmaSquared
    global angleDOF
    global angles
    # global diagFile
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
    global end_to_end_distance_squared
    global radius_of_gyration_squared
    if multiplePolymerLengths:
        nBeads = nBeadsVar
    polymersAdded = 0
    completedPolymers = 0
    nPolymersEnrich = startPolymers
    angleDOF = 6                              # Amount of different angles the polymer can move in
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
    # diagFile = open("diag.txt","w") # Open file for diagnostic data
    
    while N < nPolymers: # Simulate polymers
        # Variables
        polymerWeightInit = 1   # Initial polymer weight factor
        L = 2                   # Starting bead
        beadPos[1,:] = [1,0]
        polymerWeights[N] = 0
        # Grow beads recursively
        beadPos, polymerWeights[N], N = addBead(L,N,beadPos,candidatePos,angles,angleLastBead,polymerWeightInit)
        # sys.stdout.write("Polymer: %d of %d \r" % (N,nPolymers) )
        # sys.stdout.flush()

    # Diagnostic information
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
    opt_data_R = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_end_to_end_distance)
    data.save(exp_end_to_end_distance_squared,"R_squared",header="",write_mode="a",optional_data=opt_data_R)
    opt_data_G = str(completedPolymers)+" "+str(nBeads)+" "+str(sd_radius_of_gyration)
    data.save(exp_radius_of_gyration_squared,"radius_of_gyration",header="",write_mode="a",optional_data=opt_data_G)
    
    
    
    # print("End to end distance: " + str(exp_end_to_end_distance))
    # print("Radius of gyration: " + str(exp_radius_of_gyration))
    # print("End to end distance squared: " + str(exp_end_to_end_distance_squared))
    # print("Radius of gyration squared: " + str(exp_radius_of_gyration_squared))
    # print("End to end distance SD: " + str(np.sqrt(exp_end_to_end_distance_squared-np.square(exp_end_to_end_distance))))
    # print("Radius of gyration SD: " + str(np.sqrt(exp_radius_of_gyration_squared-np.square(exp_radius_of_gyration))))
    # print("Ratio: " + str(exp_end_to_end_distance_squared/exp_radius_of_gyration_squared))
    # polymerWeights[:] = (value for value in polymerWeights if value != 0)
    # print(polymerWeights)
    # print("number of finished POLYMERS: " + str(np.count_nonzero(polymerWeights)))
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
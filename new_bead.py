"""

"""

## Import libraries
import numpy as np
def positions(previous_beads_pos, angles1,angles2):
    offset1 = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    offset2 = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    angles1 = np.add(angles1,offset1)                      # add offset to all possible angles
    angles2 = np.add(angles2,offset2)                      # add offset to all possible angles

    possible_beads_pos = np.zeros((len(angles1),(len(angles2)),3),dtype=float)      # initialize the possible bead possitions, wrt the origin
    possible_rel_prev_beads_pos = np.zeros((len(angles1),(len(angles2)),3),dtype=float) # initialize the possible bead positions, relative to the previous bead

    for ii in range(0,len(angles1)):
        for iii in range(0,len(angles2)):
            possible_rel_prev_beads_pos[ii,iii,0] = np.sin(angles1[ii])*np.cos(angles2[iii])
            possible_rel_prev_beads_pos[ii,iii,1] = np.sin(angles1[ii])*np.sin(angles2[iii])
            possible_rel_prev_beads_pos[ii,iii,2] = np.cos(angles1[ii])

    return previous_beads_pos[:] + possible_rel_prev_beads_pos[:,:],angles1,angles2;
    
def roulette(energies,T):
    
    probabilities = np.exp(np.divide(energies,-T))   # unnormalized probabilities
    probabilities = np.divide(probabilities,sum(probabilities)) # normalize probabilities by dividing by their sum
    cumsum_probabilities = np.cumsum(probabilities)

    RNG = np.random.random() # rng is a random number chosen from a uniform distribution between 0 and 1. It is used to select one of the possible bead positions with corresponding probability

    for ii in range(0, len(probabilities)):
        if cumsum_probabilities[ii] > RNG:
            break

    return ii;
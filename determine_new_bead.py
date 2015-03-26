"""

"""

## Import libraries
import numpy as np
def determine_new_bead(energies,T):
    probabilities = np.exp(np.divide(energies,-T))   # unnormalized probabilities
    probabilities = np.divide(probabilities,sum(probabilities)) # normalize probabilities by dividing by their sum
    
    cumsum_probabilities = np.cumsum(probabilities)
    RNG = np.random.random() # rng is a random number chosen from a uniform distribution between 0 and 1. It is used to select one of the possible bead positions with corresponding probability

    for ii in range(0, len(probabilities)):
        if cumsum_probabilities[ii] > RNG:
            break

    new_bead_index = ii
    return new_bead_index;
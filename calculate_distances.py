"""
calculate_distances calculates for all points in positions_to_check what the distance is to all elements in all_positions, and returns this as the variable distances

"""

## Import libraries
import numpy as np
def calculate_distances(positions_to_check, all_positions):

    distances = np.zeros(positions_to_check.shape[0],all_positions.shape[0]),dtype=float)

    for ii in range(0, positions_to_check.shape[0]):                        # runs over each point in positions_to_check
        distances[ii,:] = np.subtract(positions_to_check,all_positions)


    return distances;
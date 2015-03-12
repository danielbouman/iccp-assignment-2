"""

"""

## Import libraries
import numpy as np
def new_bead_pos(previous_bead_pos, angles):
    offset = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    angles = np.add(angles,offset)                      # add offset to all possible angles

    possible_bead_pos = np.zeros((len(angles),2),dtype=float)
    possible_rel_prev_bead_pos = np.zeros((len(angles),2),dtype=float)

    possible_rel_prev_bead_pos[:,0] = np.cos(angles)         # calcate x-position relative to previous point
    possible_rel_prev_bead_pos[:,1] = np.sin(angles)         # calcate y-position relative to previous point

    possible_bead_pos[:,0] =  previous_bead_pos[0] + possible_rel_prev_bead_pos[:,0] # calculate absolute x - positions
    possible_bead_pos[:,1] =  previous_bead_pos[1] + possible_rel_prev_bead_pos[:,1] # calculate absolute y - positions

    return possible_bead_pos;
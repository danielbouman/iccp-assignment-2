"""

"""

## Import libraries
import numpy as np
def new_beads_pos(previous_beads_pos, angles):
    offset = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    angles = np.add(angles,offset)                      # add offset to all possible angles

    possible_beads_pos = np.zeros((len(angles),2),dtype=float)      # initialize the possible bead possitions, wrt the origin
    possible_rel_prev_beads_pos = np.zeros((len(angles),2),dtype=float) # initialize the possible bead positions, relative to the previous bead

    possible_rel_prev_beads_pos[:,0] = np.cos(angles)         # calcate x-position relative to previous point
    possible_rel_prev_beads_pos[:,1] = np.sin(angles)         # calcate y-position relative to previous point

    #possible_beads_pos[:,0] =  previous_beads_pos[0] + possible_rel_prev_beads_pos[:,0] # calculate absolute x - positions
    #possible_beads_pos[:,1] =  previous_beads_pos[1] + possible_rel_prev_beads_pos[:,1] # calculate absolute y - positions

    possible_beads_pos[:,:] =  previous_beads_pos[:] + possible_rel_prev_beads_pos[:,:] # calculate absolute position

    return possible_beads_pos;
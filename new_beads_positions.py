"""

"""

## Import libraries
import numpy as np
def new_beads_pos(previous_beads_pos, angles1,angles2):
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
    
    possible_beads_pos[:,:] =  previous_beads_pos[:] + possible_rel_prev_beads_pos[:,:] # calculate absolute position

    return possible_beads_pos;
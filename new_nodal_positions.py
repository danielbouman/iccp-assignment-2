"""

"""

## Import libraries
import numpy as np
def new_nodal_pos(previous_nodal_point, angles):
    offset = 2*np.pi*(np.random.random())               # create random offset of each angle to avoid preferential angles
    angles = np.add(angles,offset)                      # add offset to all possible angles

    pos_nodal_point = np.zeros((len(angles),2),dtype=float)
    pos_rel_prev_nodal_point = np.zeros((len(angles),2),dtype=float)

    pos_rel_prev_nodal_point[:,0] = np.cos(angles)         # calcate x-position relative to previous point
    pos_rel_prev_nodal_point[:,1] = np.sin(angles)         # calcate y-position relative to previous point

    pos_nodal_point[:,0] =  previous_nodal_point[0] + pos_rel_prev_nodal_point[:,0] # calculate absolute x - position
    pos_nodal_point[:,1] =  previous_nodal_point[1] + pos_rel_prev_nodal_point[:,1] # calculate absolute y - position

    return pos_nodal_point;
"""
List tracking

init(size)
Initialises the grid for the bead positions
size    size of the grid

store(pos,n)
Stores the location of a newly created bead in the grid.
pos     position of created bead [x,y]
n       index of created bead

get(pos,cutoff)
Gets the indexes of all beads within specified cutoff range
pos     position of current bead [x,y]
cutoff  range within beads should be considered
"""
## Import libraries
import numpy as np

def init():
    global grid
    grid = {}
    # for i in range(-size,size):       # Create grid
    #     for j in range(-size,size):
    #         for k in range(-size,size):
    #             grid['g' + str(i) + str(j) + str(k)] = []

def store(pos,n):
    grid_pos = 'g'+''.join(str(e) for e in np.floor_divide(pos,0.5).astype(int))
    grid[grid_pos] = []
    grid[grid_pos].append(n)
    
def get(pos,cutoff):
    relevant_beads = []
    current_grid_pos = np.floor_divide(pos,0.5).astype(int)
    for i in range(current_grid_pos[0]-cutoff,current_grid_pos[0]+cutoff):
        for j in range(current_grid_pos[1]-cutoff,current_grid_pos[1]+cutoff):
            for k in range(current_grid_pos[2]-cutoff,current_grid_pos[2]+cutoff):
                if 'g'+str(i)+str(j)+str(k) in grid:
                    found_beads = grid['g'+str(i)+str(j)+str(k)]
                # if found_beads:
                    relevant_beads.extend(found_beads)
    return relevant_beads
            
def show():
    print(grid)


"""
Testing:
"""
# init(10)

# store([0.1,0.5],1)
# store([0.5,1.1],2)
# store([0.4,2],3)
# store([2.2,4.3],4)
# store([2.4,4],5)
# store([2.6,4],6)
# store([2.4,3.6],7)
# store([1.4,0],8)
# store([2.4,1],9)

# rel_beads = get([2.4,4],1)
# print(rel_beads)
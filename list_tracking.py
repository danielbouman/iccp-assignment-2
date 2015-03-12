"""
List tracking

store
pos     Position [x,y]
"""
## Import libraries
import numpy as np

def init(size):
    global grid
    grid = {}
    for n in range(size):
        for m in range(size):
            grid['g' + str(n) + str(m)] = []

def store(pos,n):
    global grid_pos
    grid_pos = np.floor_divide(pos,0.5).astype(int)
    grid_pos = 'g'+''.join(str(e) for e in grid_pos)
    grid[grid_pos].append(n)
    
def get(pos,cutoff):
    relevant_beads = []
    current_grid_pos = np.floor_divide(pos,0.5).astype(int)
    for n in [current_grid_pos[0]-cutoff,current_grid_pos[0],current_grid_pos[0]+cutoff]:
        for m in [current_grid_pos[1]-cutoff,current_grid_pos[1],current_grid_pos[1]+cutoff]:
            # print(str(n)+str(m))
            found_beads = grid['g'+str(n)+str(m)]
            if found_beads:
                relevant_beads.extend(found_beads)
    return relevant_beads
            
def show():
    print(grid)


"""
Testing:
"""
init(10)

store([0.1,0.5],1)
store([0.5,1.1],2)
store([0.4,2],3)
store([2.2,4.3],4)
store([2.4,4],5)
store([2.6,4],6)
store([2.4,3.6],7)
store([1.4,0],8)
store([2.4,1],9)

rel_beads = get([2.4,4],1)
print(rel_beads)
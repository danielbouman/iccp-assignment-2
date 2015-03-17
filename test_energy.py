import energies
import numpy as np
# f = open('fortran_doc.txt', 'w')
# f.write( energies.lj_energy.__doc__ + '\n')
# f.close()
sigma = 0.8
sigma_squared = sigma*sigma
epsilon = 0.25
n_possible_pos = 6
n_beads = 10
possible_pos = np.arange(0,n_possible_pos)
pos = np.arange(0,n_beads)

for n in range(n_possible_pos):
    possible_pos[n] = 1+n
for m in range(n_beads):
    pos[m] = 1+m
    
energies1 = energies.lj_energy(possible_pos,pos,epsilon,sigma_squared,[n_possible_pos,n_beads])

# print(possible_energies)
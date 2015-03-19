import test
import numpy as np

sigma_squared = 0.64
epsilon = 0.25
N_existing = 3
N_candidates = 5
pos = np.ones((N_existing,2))
candidate_pos = np.ones((N_candidates,2))
pos[:,0] = range(3)
pos[:,1] = range(3,6)
candidate_pos[:,0] = range(5)
candidate_pos[:,1] = range(5,10)

# square = test.func(my_array,my_array_size)
energies = test.func(pos,candidate_pos,sigma_squared,epsilon,N_candidates,N_existing)
print(energies)
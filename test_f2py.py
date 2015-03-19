import test
import numpy as np

my_array = np.ones((3,2))
my_array[:,0] = range(3)
my_array[:,1] = range(3,6)
my_array_size = 3

square = test.func(my_array,my_array_size)
print(square)
import numpy as np

# # Loop over all connections once
# i = 3
# j = 3
# for n in range(i):
#     for m in range(j):
#         if n != m and n < m:
#             print(str(n) + str(m))

x = 3
y = 4

for n in [x-1,x,x+1]:
    for m in [y-1,y,y+1]:
        print(str(n)+str(m))
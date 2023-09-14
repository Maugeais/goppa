import sys
sys.path.insert(0, '..')

import numpy as np
from goppa import error_control


G = np.matrix([[1, 1, 0, 1], 
               [1, 0, 1, 1], 
               [1, 0, 0, 0], 
               [0, 1, 1, 1], 
               [0, 1, 0, 0], 
               [0, 0, 1, 0], 
               [0, 0, 0, 1]], dtype = int)


# Creation of the code defined by the above generator matrix
E = error_control.code(G)

# Random code with random error
m = np.random.randint(0, 2, G.shape[1])
c = E.code(m)
e = np.zeros_like(c, dtype = int)
e[np.random.choice([i for i in range(E.G.shape[1])], E.t)] = 1
cp = (c+e) % 2

# Decoding
print(E.decode(cp), m)



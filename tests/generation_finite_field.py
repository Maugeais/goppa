#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field
import time
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

N = 10

I = []
T = []
M = []
m = []

def par_generate(a) :
    
    t0 = time.time()
    
    F = finite_field.field(2, a[0], seed = a[1])

    t = time.time()-t0
    
    return(t)

for i in range(2, 40) :
    
    print("Generation of the field GF(2**{})...".format(i))

    with multiprocessing.Pool() as pool:
        
        t = pool.map(par_generate, [(i, n) for n in range(N*i)])
    
    T.append(np.average(t))
    M.append(max(t))
    m.append(min(t))
    
    I.append(i)
        
    
plt.scatter(I, np.array(T), label="average")
# plt.plot(I, np.array(m), label="best")
# plt.plot(I, np.array(M), label="worst")


z = np.polyfit(I, T, 4)
p = np.poly1d(z)
plt.plot(I, p(I), "--")

plt.grid(True)
plt.xlabel("Dimension over the prime subfield")
plt.ylabel("Time (s)")
plt.legend()
plt.savefig("generation_finite_field.svg")
# plt.show()

 #!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field, polynomials_finite_field
import time
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing


def par_generate(a) :
    
    t0 = time.time()
    
    P = polynomials_finite_field.genIrred(a[0], a[1], seed = a[2])

    t = time.time()-t0
        
    return(t)

N = 40

for j in range(7, 8) :
    
    I = []
    T = []
    
    print("Irreducible polynomial of degree {}...".format(j))
    
    for i in range(1, 13) :

        print("       generation of the field GF(2**{})...".format(i))
    
        F = finite_field.field(2, i)
        
        
        with multiprocessing.Pool() as pool:
            
            t = pool.map(par_generate, [(F, j, n) for n in range(N)])

     
        T.append(np.average(t))
        I.append(i)

    plt.scatter(I, T, label="Degree {}".format(j))
    z = np.polyfit(I, T, 1)
    p = np.poly1d(z)
    plt.plot(I, p(I), "--")

plt.grid(True)
plt.legend()
plt.xlabel("Dimension over the prime subfield")
plt.ylabel("Time (s)")
plt.savefig("generation_irreducible.svg")

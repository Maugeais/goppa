#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field, polynomials_finite_field
import time
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

N = 10

dimFieldMax = 2
degMax = 30

def par_irred(a) :
    
    P = X**a[0]
 
    P.coef[:-1] =  F.rand(a[0], seed = a[1])
        
    t0 = time.time()

    P.isIrred()
    
    return(time.time()-t0)

for i in range(1, dimFieldMax+1) :

    print("Generation of the field GF(2**{})...".format(i))

    F = finite_field.field(2, i)
    X = polynomials_finite_field.pol([0, 1], F)


    I = []
    T = []

    for j in range(1, degMax+1) :

        print("       irreducible polynomial of degree {}...".format(j))
                    
        with multiprocessing.Pool() as pool:
    
            t = pool.map(par_irred, [(j, n) for n in range(N*j)])
            
           
        T.append(np.average(t))
        I.append(j)

    plt.scatter(I, T, label="field cardinal = {}".format(2**i))
    z = np.polyfit(I, T, 3)
    p = np.poly1d(z)
    plt.plot(I, p(I), "--")

plt.grid(True)
plt.legend()
plt.xlabel("Degree of the polynomial")
plt.ylabel("Time (s)")
plt.savefig("test_irreducibility.svg")

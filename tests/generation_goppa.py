#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field, polynomials_finite_field, goppa
import time
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np


plt.figure("Irreducible polynomials")

N = 10

maxDim = 11
maxCor = 8

def par_generate(a) :
    
    t0 = time.time()
    
    G = goppa.goppa(a[0].field, a[0], a[0].field.card, seed = a[1])

    t = time.time()-t0
        
    return(t)



for j in range(2, 5) :

    Itemp = []
    Ttemp = []
    I = []
    T = []
    
    gen = (i for i in range(3, 8) if 2**i > j*i)
    for i in gen :

        print("Generation of the field GF(2**{})...".format(i))

        F = finite_field.field(2, i)
        
        
        print("       irreducible polynomial of degree {}...".format(j))
        
        P = polynomials_finite_field.genIrred(F, j)
        
        for a in F.elements() :
            if P(a) == 0 :
                print(P, a, P(a))
        
        

        with multiprocessing.Pool() as pool:
            
            t = pool.map(par_generate, [(P, n) for n in range(N)])
        
        T.append(np.average(t))
        I.append(2**i)
        
        # Ttemp.append(np.average(t))
        # Itemp.append(i)
        
    plt.scatter(I, T, label="Correction capacity = {}".format(j))
    
        
    z = np.polyfit(I, (T), 2)
    p = np.poly1d(z)
    I = np.linspace(min(I), max(I), 100)
    plt.plot(I, (p(I)), "--")

plt.grid(True)
plt.legend()
plt.xlabel("Cardinality of L")
plt.ylabel("Time (s)")
plt.savefig("generation_goppa.svg")

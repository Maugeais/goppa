#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field, polynomials_prime_field
import time
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing

N = 10

dimFieldMax = 2
degMax = 50

def par_irred(a) :
    
    P = X**a[0]
    
    np.random.seed((time.time_ns() + a[1]) % 2**32)
    
    P.coef[:-1] = np.random.randint(0, 2, j)    

    t0 = time.time()

    P.isIrred()
    
    return(time.time()-t0)

for i in range(2, dimFieldMax+2) :

    print("Generation of the field Z/{}Z...".format(i))

    X = polynomials_prime_field.polZnZ([0, 1], i)


    I = []
    T = []

    for j in range(1, degMax+1) :

        print("       irreducible polynomial of degree {}...".format(j))
                    
        with multiprocessing.Pool() as pool:
    
            t = pool.map(par_irred, [(j, n) for n in range(N*j)])
            
           
        T.append(np.average(t))
        I.append(j)

    plt.scatter(I, T, label="field cardinal = {}".format(i))
    z = np.polyfit(I, T, 3)
    p = np.poly1d(z)
    plt.plot(I, p(I), "--")

plt.grid(True)
plt.legend()
plt.xlabel("Degree of the polynomial")
plt.ylabel("Time (s)")
plt.savefig("test_irreducibility_prime.svg")

#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import goppa, error_control
import numpy as np
import time

print("Generation of the binary Goppa code...")


t0 = time.time()


G = goppa.goppa(12, 8, 2**11)

print(G)
print("Generation time : {:.2f}s".format(time.time()-t0))
print("Size of the syndromes table: {:.2f} Go".format(error_control.ballSize(G.H.shape[0], G.g.deg*2+1)/(2**30)))

t0 = time.time()


for i in range(1) :
    # # Création d'un vecteur 
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
        
    error = np.zeros(G.H.shape[1], dtype = int)
    I = np.random.choice(range(G.H.shape[1]), G.g.deg, replace = False)
    error[I] = 1
        
    t0 = time.time()
    
    mp = G.decode((c+error))
    
    print("Decoding is = ", sum(abs(mp-m)) == 0)
    print("    Decoding time {:.2f}s".format(time.time()-t0 ))
    
   



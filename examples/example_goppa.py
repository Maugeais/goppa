#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import goppa, error_control
import numpy as np
import time

print("Génération du code de Goppa...")

G = goppa.goppa(7, 8, 2**6)
# G = goppa.goppa(5, 4, 2**5)
print(G)

print("Table des syndromes...")
print("     Taille de la table {} Go".format(G.H.shape[0],error_control.ballSize(G.H.shape[0], G.g.deg*2+1)/(2**30)))

t0 = time.time()


for i in range(1) :
    # # Création d'un vecteur 
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
    
    # s = G.H.dot(X) % 2
    
    error = np.zeros(G.H.shape[1], dtype = int)
    I = np.random.choice(range(G.H.shape[1]), G.g.deg, replace = False)
    error[I] = 1
        
    t0 = time.time()
    
    mp = G.decode((c+error))
    
    print("Difference = ", sum(mp-m), "Temps de calcul", time.time()-t0 )
    
   



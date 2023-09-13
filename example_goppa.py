from goppa import goppa, error_control
import numpy as np
import time

print("Génération du code de Goppa...")

G = goppa.goppa(7, 8, 2**6)
# G = goppa.goppa(5, 4, 2**5)
print(G)

print("Table des syndromes...")
print("     Taille de la table {} Go".format(G.dimSyndrome,error_control.ballSize(G.dimSyndrome, G.g.deg*2+1)/(2**30)))

t0 = time.time()

# syndTable = error_control.syndrom(G.dimSyndrome, G.g.deg)
# print("     Temps de calcul", time.time()-t0 )



for i in range(1) :
    # # Création d'un vecteur 
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
    
    # s = G.H.dot(X) % 2
    
    error = np.zeros(G.H.shape[1], dtype = int)
    I = np.random.choice(range(G.H.shape[1]), G.g.deg, replace = False)
    error[I] = 1
        
    t0 = time.time()
    
    cp = G._patterson((c+error))
    
    print("Différence = ", sum(cp-c), "Temps de calcul", time.time()-t0 )
    
    t0 = time.time()
    
    synd = G.H.dot(c+error)
    
    # for s in syndTable :
    #     if all(s == synd) :
    #         print("trouvé !")
    #         break
    
    print("Différence = ", sum(cp-c), "Temps de calcul", time.time()-t0 )


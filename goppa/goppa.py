#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 10:19:30 2023

@author: maugeais
"""

from . import finite_field
from . import polynomials_finite_field
import numpy as np
from . import error_control


def traceMatrix(H) :
    
    
    m, n = H.shape
    dim = H[0, 0].F.dim
            
    Hp = np.zeros((m*dim, n), dtype = int)
    
    for i in range(m) :
        for j in range(n) :
                       
            Hp[i*dim:(i+1)*dim, j] = np.pad(H[i][j].val.coef, (0, dim-1-max(H[i][j].val.deg, 0)), 'constant')
            
    return(Hp)

class goppa :
    def __init__(self, m, t, n, pol = None) :
        """
        Initialisation de la classe goppa

        Parameters
        ----------
        m : int
            dimension du corps de base sur F_2, le cardinal est donc 2**m
        t : int
            degré du polynôme générateur
        n : int
            nombre de points considérés
        pol : TYPE, optional
            DESCRIPTION. The default is None.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
                 
        if (n-t*m <= 0) :
            raise Exception("La dimension doit être positive")
            
        self.F = finite_field.field(2, m)
        
        # Tirace de n éléments SANS remise
        self.L = np.random.choice(self.F.elements(), n, replace = False)
        self.X = polynomials_finite_field.pol([0, 1], self.F)
        self.t = t
        
        self.dimCode = n-t*m
        self.dimSyndrome = t*m
        self.length = n

        if pol != None :
            self.g = polynomials_finite_field.pol(pol, self.F)
        else :
            # Trouve un polynôme unitaire qui ne s'annule sur aucun des elements
            P = self.X**t
                    
            while not P.isIrred() :
                            
                for j in range(t) :
                    P.coef[j] = self.F.rand()
                                        
                
                # if P.isIrred :
                #     break
                # dP = P.der()
                # Q = polynomials_finite_field.gcd(P, dP)
                
                # vals = [P(gamma) == 0 for gamma in self.L]   
                # print("Q irred", P.isIrred())

                # # Si le polynôme n'a pas de racines multiples
                # if Q.deg == 0 and not any(vals):        
                #     break
                
            self.g = P
        
        h = [1/self.g(gamma) for gamma in self.L]
                
                            
        self.Hf = np.array([[h[i]*self.L[i]**j for i in range(n)] for j in range(t)])
        
        self.H = traceMatrix(self.Hf)
        
        A, B = error_control.inv(self.H.T)
        
        # Matrice génératrice
        
        self.G = B.T
        
        # print("Distance de Hamming", error_control.d_ham(np.matrix(self.G)))
        
    def __repr__(self) :
        s  = "Code de goppa \n"
        s += "   Corps GF({})\n".format(self.F.card)
        s += "   polynôme : {}\n".format(self.g)
        s += "   dimension : {}\n".format(self.dimCode)
        s += "   correction : {}\n".format(self.g.deg)
        s += "   type : ({}, {}, {})\n".format(self.length, self.dimCode, self.g.deg)
        
        return(s)
        
        
        
#     def _vector2field(self, x) :
#
#         if (len(x) != self.dimCode*self.F.card//2) :
#             print("problem de dimension: len(x) = ", self.dimCode*self.F.card//2)
#
#         res = [prime_field.GF(x[i*self.F.dim:(i+1)*self.F.dim], self.F) for i in range(len(x)//self.F.dim)]
#
#         return(res)
    
    def _syndrome(self, Y) :
        """
        Calcul du syndromze 

        Parameters
        ----------
        Y : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
             
        s = 0*self.X
        
        for i, gamma in enumerate(self.L) :
            
            if Y[i] != 0 :
                _, Q, _ = polynomials_finite_field.bezout(self.X-gamma, self.g)

                s += Q
                
        return(s)
    
    def _patterson(self, Y) :
        """
        Algorithme de Patterson pour décoder un syndrome

        Parameters
        ----------
        Y : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        s = self._syndrome(Y % 2)
                                
        if s == 0 :
            return(Y)
        
                     
        _, sm1, _ =  polynomials_finite_field.bezout(s, self.g)
        
        sqrt = sqrtModg(sm1-self.X, self.g)
        
        a, b = LBR((sqrt, polynomials_finite_field.pol([1], self.F)),
                    (self.g, polynomials_finite_field.pol([0], self.F)), length = lambda x : (x[0]**2+self.X*x[1]**2).deg)
                        
        p = (a**2+self.X*b**2) 
        
        err = p._zeros(self.L)
        
        cor = np.zeros_like(Y, dtype = int)
        cor[err] = 1
                
        return((Y-cor) % 2)
         
    def decode(self, s) :         
                
        e = G._patterson(s)
        
        return(e)
        
        
        
                
def sqrtModg(a, g) :
    """
    Calcul de la racine sm1-x mod g : on est dans une extension de coprs fini en car 2

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    g : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # 
       
    
    n = a.field.card**g.deg//4
    
    ap = a
    
      
    while n > 0 :
           
         ap = (ap*ap) % g
         
         n //= 2
  
    
    # print('Vérif du carré = ', (ap**2-a) % g)
    return(ap)

def LBR(u, v, length) :
    """
    Lattice basis reduction

    Parameters
    ----------
    a : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    if length(u) < length(v) :
        u, v = v, u
        
        
    while length(v) < length(u) :
             
        if length(u) % 2 == 0 :
            q = u[0]//v[0]
        else :
            q = u[1]//v[1]
                
        r = (u[0]-q*v[0], u[1]-q*v[1])
        u = v
        v = r
                
    return(u)

def testGoppaStandard(G) :
    
    # code de goppa standard (8, 2, 5)
    G = goppa(3, 2, 8, pol = [1, 1, 1])
    # G = goppa(4, 3, 2**4)
    
    
    # # Création d'un vecteur 
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
    
    # s = G.H.dot(X) % 2
    
    error = [0, 0, 0, 0, 0, 1, 1, 0] 
    
    cp = G._patterson((c+error))
    cp = G._patterson(0*c+error)
    
    print(c, cp)
    
def test(G) :

    # # Création d'un vecteur 
    m = np.random.randint(0, 2, G.G.shape[1])
    c = G.G.dot(m) % 2
    
    # s = G.H.dot(X) % 2
    
    error = np.zeros(G.H.shape[1], dtype = int)
    I = np.random.choice(range(G.H.shape[1]), G.g.deg, replace = False)
    error[I] = 1
    
    print("Nombre d'erreur = ", sum(error))
    
    cp = G._patterson((c+error))
    
    print("Différence = ", sum(cp-c))
    
if __name__ == "__main__" :
    
    G = goppa(6, 4, 2**5)
    
    print(G)
    
    
    test(G)
    

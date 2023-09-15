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
            the base field is GF(2**m)
        t : int
            degree of the polynomial
        n : int
            number of points in the list L
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
            raise Exception("The dimension must be positive")
            
        self.F = finite_field.field(2, m)
        
        # Tirace de n éléments SANS remise
        self.L = np.random.choice(self.F.elements(), n, replace = False)
        self.X = polynomials_finite_field.pol([0, 1], self.F)
        self.t = t
        
        self.k = n-t*m
        self.length = n

        if pol != None :
            self.g = polynomials_finite_field.pol(pol, self.F)
        else :
            # Find random irreducible polynomial of degree t
            P = self.X**t
                    
            while not P.isIrred() :
                            
                for j in range(t) :
                    P.coef[j] = self.F.rand()
                                        
                
            self.g = P
        
        h = [1/self.g(gamma) for gamma in self.L]
                
                            
        self.Hf = np.array([[h[i]*self.L[i]**j for i in range(n)] for j in range(t)])
        
        self.H = traceMatrix(self.Hf)
        
        _, B = error_control.inv(self.H.T)
        
        # Generator matrix
        self.G = B.T
        
        # Decoding matrix
        self.D, _ = error_control.inv(self.G)
        
        # print("Distance de Hamming", error_control.d_ham(np.matrix(self.G)))
        
    def __repr__(self) :
        s  = "Goppa code \n"
        s += "   Corps GF({})\n".format(self.F.card)
        s += "   polynôme : {}\n".format(self.g)
        s += "   dimension : {}\n".format(self.k)
        s += "   correction : {}\n".format(self.g.deg)
        s += "   type : ({}, {}, {})\n".format(self.length, self.k, 2*(self.g.deg)+1)
        
        return(s)
        
        
        
#     def _vector2field(self, x) :
#
#         if (len(x) != self.k*self.F.card//2) :
#             print("problem de dimension: len(x) = ", self.k*self.F.card//2)
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
                
        e = self._patterson(s)
        
        return(self.D.dot(e) % 2)
        
        
        
                
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


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 13:12:08 2023

@author: maugeais
"""

from .polynomials_prime_field import polZnZ, gcd, bezout, genIrred
import numpy as np
import time
    
    
class field:
    def __init__(self, car = 2, n = 4, varName = 'α', seed = 0):
        """
        Definition of the field

        Parameters
        ----------
        car : int
            characteristic of the field
        n : int, or 
            degree as an extension over the prime field
            car**n is the cardinal of the field
        n : list of int, or polynomials_prime_field.polZnZ
            defines a polynomial over Z/car Z
            a test of irreducibility is performed
        varName : character, optional
            name of the generator of the field. The default is 'α'.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        
        self.car = car   
        self.varName = varName
        
        if isinstance(n, list) :
            n = polZnZ(n, self.car)
        
        # If n is a polynomial
        if isinstance(n, polZnZ) :
            
                self.gen = n 
                self.var = GF([0, 1], self)
                
                self.dim = self.gen.deg
                self.card = car**self.dim
                
                if not self.gen.isIrred() :
                    raise Exception("The polynomial t is not irrecubible")
                    

                return              
            
        self.card = car**n
        self.dim = n
    
        P = genIrred(self.car, self.dim, seed = seed)
        
        self.gen = P   
        self.var = GF([0, 1], self)
        
    def rand(self, n = 1, seed = 0) :
        """
        Generate a random element of the current field

        Returns
        -------
        Element of the current field

        """
        
        if seed != 0 :
            np.random.seed((time.time_ns() + seed) % 2**32)
            
        if n == 1 :
            return(GF(np.random.randint(0, self.car, self.dim), self))
        return([GF(np.random.randint(0, self.car, self.dim), self) for _ in range(n)])
    
    def elements(self) :
        """
        Produces a list of all the elements of the current field

        Returns
        -------
        List

        """
        
        _elements = [GF([0], self)]
        
        for n in range(1, self.card) :
            
            L = []
            m = n
            # Construit 
            while (m != 0) :
                L.append(m % self.car)
                m = m//self.car
                
            _elements.append(GF(L, self))
            
        return(_elements)
    

        
class GF:
     """ Definition dese éléments du corps """
    
     def __init__(self, coef, F):
         """
         

         Parameters
         ----------
         coef : TYPE
             DESCRIPTION.
         F : TYPE
             DESCRIPTION.

         Raises
         ------
         Exception
             DESCRIPTION.

         Returns
         -------
         None.

         """
        
         self.F = F
                
         if isinstance(coef, list) :
             self.val = polZnZ(coef, F.car) % F.gen
            
         elif isinstance(coef, int) :
             self.val = polZnZ(coef, F.car) 
            
         elif isinstance(coef, np.ndarray) :
             self.val = polZnZ(coef, F.car) % F.gen            
            
         elif isinstance(coef, polZnZ) :
             self.val = coef % F.gen
            
         else :
             raise Exception("Type non reconnu pour la définition d'un élément GF")
        
        
     def __repr__(self):
         """
         Display the element

         Returns
         -------
         None.

         """
         
         return(self.val.__repr__(var = self.F.varName))


     # Addition de polynomes
     def __add__(self, P):
         """
         

         Parameters
         ----------
         P : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
         
         if isinstance(P, int) :
             P = GF([P], self.F)
         
         return(GF((self.val+P.val)%self.F.gen, self.F))
     
     def __radd__(self, P) :
         """
         

         Parameters
         ----------
         P : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
         return(GF((P+self.val)%self.F.gen, self.F))
     

     def __sub__(self, P):
         """
         

         Parameters
         ----------
         P : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
         if isinstance(P, int) :
             P = GF([P], self.F)
         return(GF((self.val-P.val)%self.F.gen, self.F))
     
     def __rsub__(self, P):
         """
         

         Parameters
         ----------
         P : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
        
         return(-P+self)
     
     def __neg__(self): 
         """
         Opposite of the current element

         Returns
         -------
         TYPE
             DESCRIPTION.

         """
         return GF(-self.coef, self.F)    

     
     def __mul__(self, P):
         """
         Multiplication of the current element by another

         Parameters
         ----------
         P : GF
             DESCRIPTION.

         Returns
         -------
         self*P

         """
         if isinstance(P, int) :
             P = GF([P], self.F)
         elif isinstance(P, GF) :
             pass
         else :
             return(P*self)
             
         return(GF((self.val*P.val)%self.F.gen, self.F))
     
     def __rmul__(self, l):
         """
         Multiplication of an interger by the current element

         Parameters
         ----------
         l : int
             DESCRIPTION.

         Returns
         -------
         l*self

         """ 
         return(self*l)
        

     
     def __pow__(self, n) :
         """ Fast computation of the power of teh current element by n
         Parameters
         ----------
         n : int
             DESCRIPTION.

         Returns
         -------
         self**n

         """
         if n >= 0 :
             
             
             a = 1
             b = self
            
             for i in reversed(bin(n)[2:]) :
                 if i == '1' :
                     a = a*b
                     
                 b = b*b
            
             return(a)
             
        
         else :
             return((1/self)**(-n))



     def __truediv__(self, B):
         """
         Division of the current element by another one

         Parameters
         ----------
         B : int or GF
             DESCRIPTION.

         Raises
         ------
         Exception
             in case B is zero

         Returns
         -------
         GF

         """
         
         if isinstance(B, int) :
             B = GF([B], self.F)
         
         if B.val.deg < 0 :
             raise Exception("Division par zéro") 
             
         _, u, _ = bezout(B.val, self.F.gen)          
     
         return(GF((self.val*u) % self.F.gen, self.F))
  
     
     def __rtruediv__(self, B) :
         """
         Division of an element by the current one

         Parameters
         ----------
         B : int
             DESCRIPTION.

         Returns
         -------
         GF

         """
         
         B = GF(B, self.F)
         
         return(B/self)
     
     
     def __eq__(self, P) :
         """
         Test the equality between the current element and P

         Parameters
         ----------
         P : int or GF
             DESCRIPTION.

         Returns
         -------
         boolean

         """
         
         
         
         if isinstance(P, int) :
            
             P = GF([P], self.F)
         
         return(self.val == P.val)
         
     def __ne__(self, P):
         """
         Test the difference between the current element to P

         Parameters
         ----------
         other : TYPE
             DESCRIPTION.

         Returns
         -------
         boolean

         """
         
         return not self.__eq__(P)
     
    
        

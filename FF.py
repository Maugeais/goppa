#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 13:12:08 2023

@author: maugeais
"""

from polZnZ import polZnZ, gcd, bezout
import numpy as np
    
    
class field:
    def __init__(self, card, car = 2, varName = 'α'):
        """
        

        Parameters
        ----------
        card : TYPE
            DESCRIPTION.
        car : TYPE, optional
            DESCRIPTION. The default is 2.
        varName : TYPE, optional
            DESCRIPTION. The default is 'α'.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        
        self.card = card
        self.car = car 
        self.n = int(np.log(card)/np.log(car))
        
        self.varName = varName
        X = polZnZ([0, 1], self.car)
                
        # Génération aléatoire de polynôme
        refs = []
        for k in range(self.n-2, self.n-1) :
            
            refs.append(X**(car**k)-X)
                        
        def isIrred(P) :
            
            for ref in refs :
                a = gcd(ref, P)
                
                if a.deg > 0 :
                    return False
            return(True)
                
        
        while True :
            
            coef = np.random.randint(0, self.car, self.n+1)
            coef[0] = coef[-1] = 1
            
            P = polZnZ(coef, self.car)
            
            # Et on vérifie que la divisibilité est bonne

            if isIrred(P) :
                break
        
        self.gen = P        
        self.var = GF([0, 1], self)
        
    def rand(self) :
        
        return(GF(np.random.randint(0, self.car, self.n), self))

        
class GF:
     """ Definition dese éléments du corps """
     def __init__(self, coef, F):
        
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
        
        
     # Affichage
     def __repr__(self):
         
         return(self.val.__repr__(var = self.F.varName))


     # Addition de polynomes
     def __add__(self, P):
         
         if isinstance(P, int) :
             P = GF([P], self.F)
         
         return(GF((self.val+P.val)%self.F.gen, self.F))
     
     def __radd__(self, P) :
         return(GF((P+self.val)%self.F.gen, self.F))
     

     def __sub__(self, P):
         if isinstance(P, int) :
             P = GF([P], self.F)
         return(GF((self.val-P.val)%self.F.gen, self.F))
     
     def __rsub__(self, P):
        
         return(-P+self)
     
     def __neg__(self): 
         """
         Négation de l'élément courant

         Returns
         -------
         TYPE
             DESCRIPTION.

         """
                  
         return GF(-self.coef, self.car)    

     
     def __mul__(self, P):
         """
         Multiplication de l'élément courant par un autre

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
         elif isinstance(P, GF) :
             pass
         else :
             return(P*self)
             
         return(GF((self.val*P.val)%self.F.gen, self.F))
     
     def __rmul__(self, l):
         """
         Multiplication d'un entier par l'élément courant

         Parameters
         ----------
         l : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """ 
         return(self*l)
        

     
     def __pow__(self, n) :
         """ À refaire avec une exponentiation rapide

         Parameters
         ----------
         n : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
         if n >= 0 :
             return(GF((self.val**n)%self.F.gen, self.F))
         else :
             return(GF(((1/self).val**(-n))%self.F.gen, self.F))


     def __truediv__(self, B):
         """
         Division de l'élément courant par un autre

         Parameters
         ----------
         B : TYPE
             DESCRIPTION.

         Raises
         ------
         Exception
             DESCRIPTION.

         Returns
         -------
         None.

         """
         
         if isinstance(B, int) :
             B = GF([B], self.F)
         
         if B.val.deg < 0 :
             raise Exception("Division par zéro") 
             
         a, u, v = bezout(B.val, self.F.gen)
         
         if (a.deg > 0) :
             print('toto', B.val, self.F.gen, gcd(B.val, self.F.gen))
                  
     
         return(GF((self.val*u) % self.F.gen, self.F))
  
     
     def __rtruediv__(self, B) :
         """
         Division d'un entier par l'élément courant

         Parameters
         ----------
         B : TYPE
             DESCRIPTION.

         Returns
         -------
         None.

         """
         
         B = GF(B, self.F)
         
         return(B/self)
     
     
     def __eq__(self, P) :
         """
         Test d'égalité

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
         
         return(self.val == P.val)
         
     def __ne__(self, other):
         """
         Test de différence

         Parameters
         ----------
         other : TYPE
             DESCRIPTION.

         Returns
         -------
         TYPE
             DESCRIPTION.

         """
         
         return not self.__eq__(other)
     
     
         
        
if __name__ == '__main__' : 

    
    F = field(5**5, 5)
    alpha = F.var
    P = (2+alpha+2*alpha**3+1*alpha**4)
    
    Q = 1/P
    print('(1/P)*P = ', Q*P)
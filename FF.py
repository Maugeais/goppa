#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 13:12:08 2023

@author: maugeais
"""

import numpy as np

class pol:
    """ Class des polynomes, reprise de rfp2 """
    def __init__(self, coef, car = 2):
        
        if isinstance(coef, str) :
            """ Ce serait bien de faire ça ..."""
            pass
        elif isinstance(coef, int) :
            coef = [coef]
        
        self.car = car            
        self.coef = np.mod(np.array(coef, dtype = int), self.car)
        I = np.where(self.coef != 0)[0]
        
        if len(I) == 0 :
            self.deg = -1
            self.coef = np.array([0], dtype = int)
        else :
            self.deg=I[-1]
            self.coef = self.coef[:self.deg+1]
        


    # Evaluation a l'aide du schema de Horner
    def eval(self, x):
        y = self.coef[self.deg]+0*x
        for i in range(self.deg-1, -1, -1):
            y = y*x+self.coef[i]

        return(y)

    # Affichage
    def __repr__(self, var = 'X'):
        
        if self.deg == -1 :
            return('0')

        if (self.deg == 0):
            poly = str(self.coef[0])
            return(poly)

        poly = ''

        for i in range(0, self.deg+1):
            if self.coef[i] != 0 :
                if self.coef[i] != 1 :
                    poly += str(self.coef[i])
                    
                elif i == 0 :
                    poly += '1'
                    
                if i == 1 :
                    poly += var
                    
                if i > 1 :
                    poly = poly+var+'**'+str(i)

                poly += '+'

        return(poly[:-1])

    # Addition de polynomes
    def __add__(self, P):
        
        if isinstance(P, int) :
            P = pol([P], self.car)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype = int)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] += P.coef[i]

        return(pol(Q, self.car))
    
    def __radd__(self, P) :
        return(self+pol(P, self.car))
    

    def __sub__(self, P):
        
        if isinstance(P, int) :
            P = pol([P], self.car)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype = int)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] -= P.coef[i]

        return(pol(Q, self.car))
    
    def __rsub__(self, P) :
        return((-1)*self+pol(P, self.car))
    
    def __mul__(self, P):
        
        if self.deg + P.deg < 0 :
            return(pol([0], self.car))
        Q = np.zeros(self.deg+P.deg+1, dtype = int)
        
        for i in range(self.deg+1) :
            for j in range(P.deg+1) :
                Q[i+j] += self.coef[i]*P.coef[j]
        
        return(pol(Q, self.car))
    
    def __pow__(self, n) :
        
        if n == 0 :
            return(pol([1], self.car))
        
        Q = self
        
        for i in range(n-1) :
            Q = Q*self
            
        return(Q)
    
    def __truediv__(self, B):
        
        
        if isinstance(B, int) :
            B = pol([B], self.car)
        
        if (B.deg < 0) :
            raise Exception("Division par zéro")
        
        
        
        R = pol(self.coef, self.car)
        Q = pol([0], self.car)
        
        while R.deg >= B.deg :
                        
            # Calcule l'inverse de B.coef[-1]
            _, u, v = bezoutint(B.coef[-1], self.car)
            
            Q = Q + (R.coef[-1]*u)*pol([0, 1], self.car)**(R.deg-B.deg)
            R = self-Q*B
 
        return(Q, R)
    
    def __mod__(self, P) :
        
        Q, R = self/P
        return(R)
        
    # Multiplication par une constante
    def __rmul__(self, l):
        return(self*pol(l, self.car))
    
    
    def __eq__(self, P) :
        
        if isinstance(P, int) :
            P = pol([P], self.car)
        
        if (self.deg != P.deg) :
            return(False)
                
        return(all(self.coef == P.coef))
    
    def __ne__(self, other):
       return not self.__eq__(other)
        
        
    
def gcd(a, b) :
    
             
    while (b.deg >= 0) :

        r = a % b
                
        a = b
        b = r
        
    return(a)

def bezoutint(a, b) :
    u0, u1, v0, v1 = 1, 0, 0, 1
    
    while (b != 0) :
        
        q = a//b
        r = a % b
        
        a = b
        b = r
        
        u0, u1 = u1, u0-u1*q
        v0, v1 = v1, v0-v1*q
        
    return(a, u0, v0)

def bezout(a, b) :
    u0, u1, v0, v1 = pol(1, a.car), pol(0, a.car), pol(0, a.car), pol(1, a.car)
    
    while (b.deg >= 0) :
        
        q, r = a/b
                
        a = b
        b = r
        
        u0, u1 = u1, u0-u1*q
        v0, v1 = v1, v0-v1*q
        
    return(a, u0, v0)
    
    
class field:
    """ Class des polynomes, reprise de rfp2 """
    def __init__(self, card, car = 2, varName = 'α'):
        
        self.card = card
        self.car = car 
        self.n = int(np.log(card)/np.log(car))
        
        self.varName = varName
        X = pol([0, 1], self.car)
                
        # Génération aléatoire de polynôme
        refs = []
        for k in range(1, self.n-1) :
            
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
            
            P = pol(coef, self.car)
            
            # Et on vérifie que la divisibilité est bonne

            if isIrred(P) :
                break
        
        self.gen = P        
        self.var = GF([0, 1], self)
        
    def rand(self) :
        
        return(GF(np.random.randint(0, self.car, self.n), self))

        
class GF:
     """ Class des polynomes, reprise de rfp2 """
     def __init__(self, coef, F):
        
        self.F = F
                
        if isinstance(coef, list) :
            self.val = pol(coef, F.car) % F.gen
            
        elif isinstance(coef, int) :
            self.val = pol(coef, F.car) 
            
        elif isinstance(coef, np.ndarray) :
            self.val = pol(coef, F.car) % F.gen            
            
        elif isinstance(coef, pol) :
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

     
     def __mul__(self, P):
         if isinstance(P, int) :
             P = GF([P], self.F)
             
         return(GF((self.val*P.val)%self.F.gen, self.F))

     
     def __pow__(self, n) :
         if n >= 0 :
             return(GF((self.val**n)%self.F.gen, self.F))
         else :
             return(GF(((1/self).val**(-n))%self.F.gen, self.F))


     def __truediv__(self, B):
         
         if isinstance(B, int) :
             B = GF([B], self.F)
         
         if B.val.deg < 0 :
             raise Exception("Division par zéro") 
             
         a, u, v = bezout(B.val, self.F.gen)
         
         if (a.deg > 0) :
             print('toto', B.val, self.F.gen, gcd(B.val, self.F.gen))
                  
         _, r, s = bezoutint(a.coef[0], self.F.car)   
                 
         for i in range(u.deg+1) :
             u.coef[i] = (u.coef[i]*r)%self.F.car
     
         return(GF((self.val*u) % self.F.gen, self.F))
  
     
     def __rtruediv__(self, B) :
         
         B = GF(B, self.F)
         
         return(B/self)
     
     # Multiplication par une constante
     def __rmul__(self, l):
         return(self*l)
     
     def __eq__(self, P) :
         
         if isinstance(P, int) :
             P = GF([P], self.F)
         
         return(self.val == P.val)
         
     def __ne__(self, other):
         
         return not self.__eq__(other)
         
        
if __name__ == '__main__' : 

    X = pol([0, 1], car = 11)
    
    P = (2+X+3*X**3+2*X**8)
    
    Q = (1+X+3*X**2)
    
    print(P+P+3*P)
    
    print(P-2*Q**4)
    
    A, B = P / Q
    
    print(P-(A*Q+B))
    
    
    F = field(5**3, 5)
    alpha = F.var
    P = (2+alpha+2*alpha**3)
    
    Q = 1/P
    print('hola', Q*P)
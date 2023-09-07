#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:15:04 2023

@author: maugeais
"""

import FF
import numpy as np

class pol:
    """ Class des polynomes, reprise de rfp2 """
    def __init__(self, coef, field):
        
        print(coef)
        
        self.coef = []
        if isinstance(coef, str) :
            """ Ce serait bien de faire ça ..."""
            pass
        elif isinstance(coef, int) :
            coef = [FF.GF(coef, field)]
        elif isinstance(coef, list) :
            
            for c in coef : 
                if isinstance(c, int) :
                    self.coef += [FF.GF(c, field)]
                if isinstance(c, FF.GF) :
                    self.coef += [c]
                    
            self.coef = np.array(self.coef, dtype=FF.GF)
            
        self.field = field
                
        
        # self.coef = np.array([FF.GF(c.val, field) for c in coef], dtype = FF.GF)
        I = np.where(self.coef != 0)[0]
        
        print(I, self.coef)
        
        if len(self.coef) == 0 or len(I) == 0 :
            self.deg = -1
            self.coef = np.array([FF.GF(0, field)])
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
            P = pol([P], self.field)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype=FF.GF)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] += P.coef[i]
            
        print("----", [q for q in Q])

        return(pol(Q, self.field))
    
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


if __name__ == '__main__' : 
    field = FF.field(2**6, 2)
    
    X = pol([0, 1], field)
    
    alpha = field.var
    
    X = pol([0, 1], field)
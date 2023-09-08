#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 12:15:04 2023

@author: maugeais
"""

import FF
import numpy as np

class pol:
    def __init__(self, coef, field):
        """

        Parameters
        ----------
        coef : TYPE
            DESCRIPTION.
        field : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
                
        self.coef = []
        if isinstance(coef, str) :
            """ Ce serait bien de faire ça ..."""
            pass
        elif isinstance(coef, int) :
            self.coef = [FF.GF(coef, field)]

        elif isinstance(coef, list) or isinstance(coef, np.ndarray):
            
            for c in coef : 
                if isinstance(c, int) :
                    self.coef += [FF.GF(c, field)]
                if isinstance(c, FF.GF) :
                    self.coef += [c]
                    
            self.coef = np.array(self.coef, dtype=FF.GF)
            
        self.field = field
                
        
        # self.coef = np.array([FF.GF(c.val, field) for c in coef], dtype = FF.GF)
        I = np.where(self.coef != 0)[0]
                
        if len(self.coef) == 0 or len(I) == 0 :
            self.deg = -1
            self.coef = np.array([FF.GF(0, field)])
        else :
            self.deg=I[-1]
            self.coef = self.coef[:self.deg+1]
        


    # Evaluation a l'aide du schema de Horner
    def __call__(self, x):
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
                    I = np.where(self.coef[i].val.coef != 0)[0]
                    if (len(I) > 1) :
                        
                        poly += '('+str(self.coef[i])+')'
                        
                    else :
                        poly += str(self.coef[i])
                        
                    poly += '.'
                    
                elif i == 0 :
                    poly += '1'
                    
                if i == 1 :
                    poly += var
                    
                if i > 1 :
                    poly = poly+var+'**'+str(i)

                poly += '+'

        return(poly[:-1])

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
            P = pol([P], self.field)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype=FF.GF)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] += P.coef[i]
            
        return(pol(Q, self.field))
    
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
        return(self+pol(P, self.field))
    

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

        return(self+(-1)*P)
    
    def __rsub__(self, P) :
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        return(-self+P)
    
    def __mul__(self, P):
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if isinstance(P, FF.GF) :
            Q = np.copy(self.coef)
            for i in range(self.deg+1) :
                Q[i] *= P
            return(pol(Q, self.field))
        
        if self.deg + P.deg < 0 :
            return(pol([0], self.field))
        Q = np.zeros(self.deg+P.deg+1, dtype = FF.GF)
        
        for i in range(self.deg+1) :
            for j in range(P.deg+1) :
                Q[i+j] += self.coef[i]*P.coef[j]
                        
        return(pol(Q, self.field))
    
    # Multiplication par une constante
    def __rmul__(self, l):
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        return(self*pol(l, self.field))
    
    def __pow__(self, n) :
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if n == 0 :
            return(pol([1], self.field))
        
        Q = self
        
        for i in range(n-1) :
            Q = Q*self
            
        return(Q)
    
    def __truediv__(self, P):
        
        
        if isinstance(P, int) or isinstance(P, FF.GF) :
            P = pol([P], self.car)
        
        if (P.deg < 0) :
            raise Exception("Division par zéro")
        
        
        R = pol(self.coef, self.field)
        Q = pol([0], self.field)
        X = pol([0, 1], self.field)
        while R.deg >= P.deg :
                        
            # Calcule l'inverse de B.coef[-1]
            u = 1/P.coef[-1]
            
            Q = Q + (R.coef[-1]*u)*X**(R.deg-P.deg)
            R = self-Q*P
 
        return(Q, R)
    
    def __mod__(self, P) :
        
        Q, R = self/P
        return(R)
    
    def __floordiv__(self, P) :
        
        Q, R = self/P
        return(Q)
        
    
    
def gcd(a, b) :
    
             
    while (b.deg >= 0) :

        r = a % b
                
        a = b
        b = r
        
    return(a)


def bezout(a, b) :
    u0, u1, v0, v1 = pol(1, a.field), pol(0, a.field), pol(0, a.field), pol(1, a.field)
    
    while (b.deg >= 0) :
        
        q, r = a/b
        
        a = b
        b = r
        
        u0, u1 = u1, u0-u1*q
        v0, v1 = v1, v0-v1*q
        
    # # Et on normalise pour avoir un pgcd unitaire
    if a.deg >= 0 :
        
        r = 1/a.coef[-1]
                
        for i in range(u0.deg+1) :
            u0.coef[i] *= r
            
        for i in range(v0.deg+1) :
            v0.coef[i] *= r
            
        for i in range(a.deg+1) :
            a.coef[i] *= r
                
        
    return(a, u0, v0)


if __name__ == '__main__' : 
    field = FF.field(5**4, 5)
    
    X = pol([0, 1], field)
    
    alpha = field.var
    
    P = (1+(1+alpha)*X+1*X**3+1*X**7)
    
    Q = (1+X+alpha*X**3)
    
    # print(P+P+3*P)
    
    # print(P-1*Q**4)
    
    A, B = P / Q
    print("P-(A*Q+B) = ", P-(A*Q+B))
    
    a, u, v = bezout(P, Q)
    print(a, '=', u*P+v*Q)

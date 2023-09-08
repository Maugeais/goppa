#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 08:08:05 2023

@author: maugeais
"""

import numpy as np

class polZnZ:
    def __init__(self, coef, car = 2):
        """
        

        Parameters
        ----------
        coef : TYPE
            DESCRIPTION.
        car : TYPE, optional
            DESCRIPTION. The default is 2.

        Returns
        -------
        None.

        """
        
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
        

    def __call__(self, x):
        """
        Evaluation a l'aide du schema de Horner  

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        y = self.coef[self.deg]+0*x
        for i in range(self.deg-1, -1, -1):
            y = (y*x+self.coef[i]) % self.car

        return(y)

    # Affichage
    def __repr__(self, var = 'X'):
        """
        

        Parameters
        ----------
        var : TYPE, optional
            DESCRIPTION. The default is 'X'.

        Returns
        -------
        None.

        """
        
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

    
    def __add__(self, P):
        """
        # Addition de polynomes
        
        Parameters
        ----------
        var : TYPE, optional
            DESCRIPTION. The default is 'X'.

        Returns
        -------
        None.

        """
        
        if isinstance(P, int) :
            P = polZnZ([P], self.car)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype = int)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] += P.coef[i]

        return(polZnZ(Q, self.car))
    
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
        return(self+polZnZ(P, self.car))
    

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
        return((-1)*self+polZnZ(P, self.car))
    
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
        
        if self.deg + P.deg < 0 :
            return(polZnZ([0], self.car))
        Q = np.zeros(self.deg+P.deg+1, dtype = int)
        
        for i in range(self.deg+1) :
            for j in range(P.deg+1) :
                Q[i+j] += self.coef[i]*P.coef[j]
        
        return(polZnZ(Q, self.car))
    
    def __rmul__(self, l):
        """
        

        Parameters
        ----------
        l : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        return(self*polZnZ(l, self.car))
    
    def __pow__(self, n) :
        """
        Il faudrait refaire avec exponentiation rapide

        Parameters
        ----------
        n : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        
        if n == 0 :
            return(polZnZ([1], self.car))
        
        Q = self
        
        for i in range(n-1) :
            Q = Q*self
            
        return(Q)
    
    def __truediv__(self, P):
        
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
            P = polZnZ([P], self.car)
        
        if (P.deg < 0) :
            raise Exception("Division par zéro")
        
        
        
        R = polZnZ(self.coef, self.car)
        Q = polZnZ([0], self.car)
        
        while R.deg >= P.deg :
                        
            # Calcule l'inverse de B.coef[-1]
            _, u, v = bezoutint(P.coef[-1], self.car)
            
            Q = Q + (R.coef[-1]*u)*polZnZ([0, 1], self.car)**(R.deg-P.deg)
            R = self-Q*P
 
        return(Q, R)
    
    def __mod__(self, P) :
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        Q, R = self/P
        return(R)
    
    def __floordiv__(self, P) :
        
        Q, R = self/P
        return(Q)
    
    
    def __eq__(self, P) :
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
            P = polZnZ([P], self.car)
        
        if (self.deg != P.deg) :
            return(False)
                
        return(all(self.coef == P.coef))
    
    def __ne__(self, P):
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        return not self.__eq__(P)
        
    
    def __neg__(self):
        """
        

        Parameters
        ----------
        None

        Returns
        -------
        None.

        """
        return polZnZ(-self.coef, self.car)
        
    
def gcd(a, b) :
    """
    

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
    
             
    while (b.deg >= 0) :

        r = a % b
                
        a = b
        b = r
        
    return(a)

def bezoutint(a, b) :
    """
    

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
    """
    

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
    u0, u1, v0, v1 = polZnZ(1, a.car), polZnZ(0, a.car), polZnZ(0, a.car), polZnZ(1, a.car)
    
    while (b.deg >= 0) :
        
        q, r = a/b
                
        a = b
        b = r
        
        u0, u1 = u1, u0-u1*q
        v0, v1 = v1, v0-v1*q
        
    # Et on normalise pour avoir un pgcd unitaire
    if a.deg >= 0 :
        
        _, r, s = bezoutint(a.coef[-1], a.car)   
                
        for i in range(u0.deg+1) :
            u0.coef[i] = (u0.coef[i]*r)%a.car
            
        for i in range(v0.deg+1) :
            v0.coef[i] = (v0.coef[i]*r)%a.car
            
        for i in range(a.deg+1) :
            a.coef[i] = (a.coef[i]*r)%a.car
        
    return(a, u0, v0)

if __name__ == '__main__' : 

    X = polZnZ([0, 1], car = 11)
    
    P = (2+X+3*X**3+2*X**8)
    
    Q = (1+X+3*X**2)
    
    print(P+P+3*P)
    
    print(P-2*Q**4)
    
    A, B = P / Q
    
    print(P-(A*Q+B))
    
    a, u, v = bezout(P, Q)
    print(a, '=', u*P+v*Q)
    
    
    
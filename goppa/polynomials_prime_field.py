#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 08:08:05 2023

@author: maugeais
"""

import numpy as np
import sympy

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
            coef = np.array([coef], dtype = int)
        elif isinstance(coef, int) :
            coef = np.array(coef, dtype = int)
        
        self.car = car            
        self.coef = np.mod(coef, self.car)
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
        coef = np.copy(self.coef)
        coef[0] = (coef[0]+P) % self.car
        
        return(polZnZ(coef, self.car))
    

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
            P = polZnZ([P], self.car)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype = int)
        for i in range(self.deg+1):
            Q[i] = self.coef[i]

        for i in range(P.deg+1):
            Q[i] -= P.coef[i]

        return(polZnZ(Q, self.car))

    
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
        coef = np.copy(-self.coef)
        coef[0] = (coef[0]+P) % self.car
        
        return(polZnZ(coef, self.car))
        # return((-1)*self+polZnZ(P, self.car))
    
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
            
            if R.deg > P.deg :
                S = np.zeros(R.deg-P.deg+1, dtype = int)
                S[-1] = R.coef[-1]*u
                S = polZnZ(S, self.car)
            else :
                S = polZnZ([R.coef[-1]*u], self.car)
            
            Q = Q + S
        
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
        
        _, R = self/P
        return(R)
    
    def __floordiv__(self, P) :
        """
        

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        Q, _ = self/P
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
    
    def powMod(self, n, g) :
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
            return(polZnZ([1], self.car))
        
        a = 1
        b = self
       
        for i in reversed(bin(n)[2:]) :
            if i == '1' :
                a = (a*b) %g
                
            b = (b*b) % g
       
        return(a) 
    
    def isIrred(self) :
        """
        Test d'irréductibilité de Rabin

        Returns
        -------
        bool
            DESCRIPTION.

        """
        N = sympy.primefactors(self.deg)
        X = polZnZ([0, 1], self.car)
        
        for n in N :
            
            h = X.powMod(self.car**(self.deg//n), self)-X
    
            if gcd(self, h).deg > 0 :
                return False
                
        g = X.powMod(self.car**(self.deg), self)-X
    
        
        if (g % self).deg < 0 :
            return(True)
        
        return(False)
        
    
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

    X = polZnZ([0, 1], car = 2)
    
    P = (2+X+3*X**3+2*X**8)
    P = 1+X+X**2
    
    # Q = (1+X+3*X**2)
    
    # print(P+P+3*P)
    
    # print(P-2*Q**4)
    
    # A, B = P / Q
    
    # print(P-(A*Q+B))
    
    # a, u, v = bezout(P, Q)
    # print(a, '=', u*P+v*Q)
    
    
    print(P.isIrred())
    
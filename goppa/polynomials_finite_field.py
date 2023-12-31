from . import finite_field
import numpy as np
from . import arithmetic_tools
import time 


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
            self.coef = [finite_field.GF(coef, field)]

        elif isinstance(coef, list) or isinstance(coef, np.ndarray):
            for c in coef : 
                if isinstance(c, int) or isinstance(c, np.int64) :            
                    self.coef += [finite_field.GF(int(c), field)]
                if isinstance(c, finite_field.GF) :
                    self.coef += [c]
                    
            self.coef = np.array(self.coef, dtype=finite_field.GF)
            
        self.field = field
                
        # print(self.coef)
        # self.coef = np.array([finite_field.GF(c.val, field) for c in coef], dtype = finite_field.GF)
        I = np.where(self.coef != 0)[0]
                
        if len(self.coef) == 0 or len(I) == 0 :
            self.deg = -1
            self.coef = np.array([finite_field.GF(0, field)])
        else :
            self.deg = int(I[-1])
            self.coef = self.coef[:self.deg+1]
        


    
    def __call__(self, x):
        """
        Evaluation du polynôme a l'aide du schema de Horner

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
            y = y*x+self.coef[i]

        return(y)

    def __repr__(self, var = 'X'):
        """
        Affichage

        Parameters
        ----------
        x : TYPE
            DESCRIPTION.

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
        Addition du polynôme courant avec un autre

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if isinstance(P, int) or isinstance(P, finite_field.GF):
            P = pol([P], self.field)
            
        Q = np.zeros(max(self.deg, P.deg)+1, dtype=finite_field.GF)
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
        
        return((-1)*self+P)
    
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
        if isinstance(P, finite_field.GF) :
            Q = np.copy(self.coef)
            for i in range(self.deg+1) :
                Q[i] *= P
            return(pol(Q, self.field))
        
        if (self.deg < 0) or (P.deg < 0) :
            return(pol([0], self.field))
        
        # If P is a monommial
        if len(np.where(P.coef != 0)[0]) == 1 :
            Q = np.zeros(self.deg+P.deg+1, dtype = finite_field.GF)
                        
            Q[P.deg:] = (P.coef[-1]*self.coef) 
            return(pol(Q, self.field))
        
        # If self a monommial
        if len(np.where(self.coef != 0)[0]) == 1 :
            Q = np.zeros(self.deg+P.deg+1, dtype = finite_field.GF)
            
            Q[self.deg:] = (P.coef*self.coef[-1]) 
            return(pol(Q, self.field))
        
        
        
        Q = np.zeros(self.deg+P.deg+1, dtype = finite_field.GF)
        
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
        
        a = 1
        b = self
       
        for i in reversed(bin(n)[2:]) :
            if i == '1' :
                a = a*b
                
            b = b*b
       
        return(a)        
       
    
    def __truediv__(self, P):
        """
        Division : calcul du quotient et du reste

        Parameters
        ----------
        P : TYPE
            DESCRIPTION.

        Raises
        ------
        Exception
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if isinstance(P, int) or isinstance(P, finite_field.GF) :
            P = pol([P], self.car)
        
        if (P.deg < 0) :
            raise Exception("Division par zéro")
        
        
        R = pol(self.coef, self.field)
        Q = pol([0], self.field)
        # X = pol([0, 1], self.field)
        while R.deg >= P.deg :     
                                    
            # Calcule l'inverse de B.coef[-1]
            u = 1/(P.coef[-1])
            
            
            if R.deg > P.deg :
                S = np.zeros(R.deg-P.deg+1, dtype = finite_field.GF)
                S[-1] = R.coef[-1]*u
                S = pol(S, self.field)
            else :
                S = pol([R.coef[-1]*u], self.field)
                
            Q += S
            
            R -= S*P
            
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
        return(Q)
    
    def der(self) :
        """
        

        Returns
        -------
        None.

        """
        a = [i*self.coef[i] for i in range(1, self.deg+1)]
                
        return(pol(a, self.field))
    
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
            P = pol([P], self.field)
        
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
    
    def __neg___(self) :
        """
        

        Returns
        -------
        None.

        """
        
        return((-1)*self)
    
    def _zeros(self, L) :
        """
        Parameters
        ----------
        L : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        z = []
        
        
        for i, gamma in enumerate(L) :
            
            if self(gamma) == 0 :
                z.append(i)
                
        return(z)

        
    
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
            return(pol([1], self.field))
        
        a = 1
        b = self
       
        for i in reversed(bin(n)[2:]) :
            if i == '1' :
                a = (a*b) %g
                
            b = (b*b) % g
       
        return(a) 
    
    def isIrred(self) :
        """
        Rabin test of irreducibility
        https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Rabin.27s_test_of_irreducibility

        Returns
        -------
        None.

        """
        
        N = arithmetic_tools.prime_divisors(self.deg)
        X = pol([0, 1], self.field)
        
        for n in N :
                        
            h = X.powMod(self.field.card**(self.deg//n), self)-X

            if gcd(self, h).deg > 0 :
                
                return False
                
        g = (X.powMod(self.field.card**(self.deg), self)-X)
        
        if (g % self).deg < 0 :
            return(True)
        
        else :
            return(False)
        
    
def gcd(a, b) :
    """
    Calcul du pgcd de a et b, non normalisé !

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
             
    while (b.deg >= 0) :
        
        r = a % b
                        
        a = b
        b = r
        
    if a.deg >= 0 :
        
        r = 1/a.coef[-1]
    
        for i in range(a.deg+1) :
            a.coef[i] *= r
        
    return(a)


def bezout(a, b) :
    """
    Calcul une décomposition de Bézout normalisée de a et b

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
    u0, u1, v0, v1 = pol([1], a.field), pol([0], a.field), pol([0], a.field), pol([1], a.field)

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


def genIrred(baseField, degree, verbose = False, seed = 0) :

    X = pol([0, 1], baseField)
    P = X**degree
    
    n = 0
    
    if seed != 0 :
        np.random.seed((time.time_ns() + seed) % 2**32)
                        
    while not P.isIrred() :
        
                    
        for j in range(degree) :
            P.coef[j] = baseField.rand()
            
        n += 1
        
    if verbose : 
        print("Number of trials = ", n)
            
    return(P)

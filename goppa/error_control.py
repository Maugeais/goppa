import numpy as np
import scipy

def inv(A):
    """
    Calcul un inverse à gauche  d'une matrice A injective

    Parameters
    ----------
    A : TYPE
        DESCRIPTION.

    Raises
    ------
    Exception
        DESCRIPTION.

    Returns
    -------
    None.

    """
    m, n = A.shape
    
    Ap = np.matrix.copy(A)
    Xp = np.eye(m, dtype = int)
        
    # Triangularisation
        
    for j in range(n) :
        
        # si le pivot est nul, on doit échanger des lignes
        if (Ap[j, j] == 0 ) :
                        
            # on prend le plus grand pivot 
            
            mx = np.max(np.abs(Ap[j+1:, j]))
            if mx == 0 :
                                
                raise Exception("Matrice non injective")
            
            ell = np.argmax(np.abs(Ap[j+1:, j]))+j+1
            
            # echange la colonne j et la colonne ell
            
            tmp = np.copy(Ap[ell, :])
            Ap[ell, :] = Ap[j, :]
            Ap[j, :] = tmp
            
            tmp = np.copy(Xp[ell, :])
            Xp[ell, :] = Xp[j, :]
            Xp[j, :] = tmp
            
            # L[j], L[ell] = L[ell], L[j]
                    
        
        for i in range(j+1, m) :
            
            
            piv = Ap[i, j]//Ap[j, j]
         
            Ap[i, :] = (Ap[i, :] - piv*Ap[j, :]) % 2
            Xp[i, :] = (Xp[i, :] - piv*Xp[j, :]) % 2
                                    
    # Remontée
                
    for j in reversed(range(n)) :
           
        for i in range(j) :
            
            piv = Ap[i, j]
            
            Ap[i, :] = (Ap[i, :] - piv*Ap[j, :])%2
            Xp[i, :] = (Xp[i, :] - piv*Xp[j, :])%2

                        
    return(Xp[:n,:], Xp[n:, :])
    

class code :
    def __init__(self, G) :
        """
        Build the error control for a given generator matrix

        Parameters
        ----------
        G : TYPE
            DESCRIPTION.

        Returns
        -------
        self.P = parity check matrix
        self.D = decoding matrix

        """
        self.G = G
        self.D, self.P = inv(G)
        
        self.d = self.d_ham(G)
        self.t = int((self.d-1)/2)
        
        L = self._syndrom_(G.shape[0], self.t)
                
        self.syndromes = [[ell, self.P.dot(ell)] for ell in L]
        
        
    def code(self, m) :
        
        y = np.array(self.G.dot(np.array(m)))[0] %2
        
        return(y)
        
    def syndrom(self, c) :
        
        y = np.array(self.P.dot(np.array(c))) % 2
        
        return(y)


    def d_ham(self, G) :
        """
        Calcule la distance de Hamming du code engendré par G

        Parameters
        ----------
        G : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        n, k = G.shape
        
        D = []
        # Génére tous les vecteurs de F_2^p
        
        for i in range(1, 2**k) :
            b = "{:0{}b}".format(i, k)
            X = np.reshape([int(a) for a in b], [k, 1])
            D.append((G*X%2).sum())
                    
        return(min(D))
    
    def _syndrom_(self, n, t) :
        """
        Calcule les syndromes

        Parameters
        ----------
        G : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if (n == 0) :
            return([[]])
            
        if (t == 0) :
            return([[0 for a in range(n)]])
            
        # Maintenant, n > 0 et t > 0
        
        L = []
    
        # On rajoute les syndrome avec premier bit = 0    
        L0 = self._syndrom_(n-1, t)
        
        for S in L0 :
            L += [[0]+S]
            
        # Puis on rajoute les syndrome commençant par 1
            
        L0 = self._syndrom_(n-1, t-1)
        
        for S in L0 :
            L += [[1]+S]
            
        return(L)
    
    def error(self, c) :
        
        s = np.array(self.P.dot(np.array(c))) % 2
        
        for S in self.syndromes :
            
            if all(s == S[1]) :
                return(S[0])
        
        raise Exception("Syndrom not found")   
        
    
    def decode(self, c) :
        
        error = self.error(c)
        
        return(self.D.dot(c+error) % 2)
    
    

def ballSize(ell, d) :
    
    t = np.floor((d-1)/2)
    
    return(int(sum(scipy.special.comb(ell, np.arange(0, t+1)))))
            
        
# trouve une matrice de taille ell, n avec une distance hamming d, maximum N étape        
def findMat(ell, n, d, N=100) :
    
    i = 0
    
    while (i < N) :
        
        G = np.matrix(np.random.randint(0, 2, (ell, n)))
        
        dt = d_ham(G)
                
        if dt == d :
            return(G)
        
        i += 1
        
    return(None)
        
    
# Test si le type d'un codez est parfait    
def isPerfect(ell, n, d) :
    
    return(2**ell == ballSize(ell, d)*2**n)

# Pour une capacité de correction donnée, cherche les longueur potentielle des codes parfaits
def potentialPerfect(t, L = 100) :
    ell = np.arange(0, L, dtype= int)
    
    if t == 1 :
    
        a = 1+ell
        
    elif t == 2 :
        a = 1+ell+ell*(ell-1)//2
        
    elif t == 3 :
        a = a = 1+ell+ell*(ell-1)//2 + ell*(ell-1)*(ell-2)//6
    
    b = 2**(np.array((np.floor(np.log2(a))), dtype=int))
    
    res = np.where(a==b)[0]

    
    n = ell[res]-np.log2(a[res]).astype(int)
    
    
    return(res, n)
    
        

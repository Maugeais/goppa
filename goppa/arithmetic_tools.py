def comb(n, k):
    
    if (k >  n) :
        raise Exception("k must be small than n")
    
    if(k == 0) :
        return 1
    
    if(k > n/2):
        return comb(n,n-k)
    
    return n * comb(n-1,k-1) // k

def prime_divisors(n) :
    """
    Compute the prime divisors of n

    Parameters
    ----------
    n : int
        DESCRIPTION.

    Returns
    -------
    List of the prime divisors

    """

    p = 2
    divisors = []
    
    while p**2 < n :
        
        if n % p == 0 :
            
            divisors.append(p)
            
            while n % p == 0 :
                n //= p
        
        if p == 2 :
            p = 3
        else :
            p += 2
           
    if n != 1 :
        divisors.append(n)
            
    return(divisors)

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
    u0, u1, v0, v1 = 1, 0, 0, 1
    
    while (b != 0) :
        
        q = a//b
        r = a % b
        
        a = b
        b = r
        
        u0, u1 = u1, u0-u1*q
        v0, v1 = v1, v0-v1*q
        
    return(a, u0, v0)

         
if __name__ == "__main__" : 
    # a = prime_divisors(2**7-1)     
    # print(a)
    
    print(comb(10, 7))
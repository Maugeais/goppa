#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import polynomials_prime_field

X = polynomials_prime_field.polZnZ([0, 1], car = 2)
    
P = (2+X+3*X**3+2*X**8)
    
Q = (1+X+3*X**2)
    
print(P-2*Q**4)
    
A, B = P / Q
    
print(P-(A*Q+B))
    
R = 1+X+X**2
    
print(R.isIrred())
    

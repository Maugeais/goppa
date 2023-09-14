#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

import sys
sys.path.insert(0, '..')

from goppa import finite_field, polynomials_finite_field

field = finite_field.field(5, 4)

X = polynomials_finite_field.pol([0, 1], field)

alpha = field.var

P = (1+(1+alpha)*X+1*X**3+1*X**7)

Q = (1+X+alpha*X**2)

# print(P+P+3*P)

# print(P-1*Q**4)

A, B = P / Q
print("P-(A*Q+B) = ", P-(A*Q+B))


print(Q.isIrred())

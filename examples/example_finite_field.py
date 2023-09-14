#!/usr/bin/env python3
import sys
sys.path.insert(0, '..')

from goppa import finite_field

F = finite_field.field(3, 5)
alpha = F.var
P = (2+alpha+2*alpha**3+1*alpha**4)

Q = 1/P
print('(1/P)*P = ', Q*P)

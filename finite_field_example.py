#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 09:43:16 2023

@author: maugeais
"""

import FF

F = FF.field(5**6, 5)
alpha = F.var


P = 1+2*alpha**2+alpha**4
Q = alpha+2*2*alpha**3

b = alpha

print((P/Q)*Q, P)

# for i in range(2**6):
    
#     if b == 1 :
#         break
#     b = b*alpha

# print(i)

# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:59:53 2024

@author: kietb
"""

import math
def simp(a, b, N):
    # Use Simpson's Rule to calculate integral
    v = a
    dv = (b-a)/N
    ans = 0
    for _ in range (N):
        v1 = v 
        v2 = v + dv
        v3 = v + dv + dv
        
        f1 = math.exp(-v1**2)* (v1*2)
        f2 = math.exp(-v2**2)* (v2**2)
        f3 = math. exp(-v3**2)* (v3**2)
        
        ans += dv / 6 * (f1 + 4*f2 +f3)
        v += dv
    return ans

print (simp(1, 10, 2000000)/ math.pi**(1/2))

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 21:07:06 2019

@author: ADMIN
"""

match = 1
mismatch = -1
gap = -2

def hallarScore(secuencia1,  secuencia2):
    score = 0
    x = 0
    for i in secuencia1:
        if i == secuencia2[x]:  
            score += match
        elif i=='-' or secuencia2[x]=='-':
            score += gap
        else:
            score += mismatch
        x = x + 1
    return score

secuencia1 = ('GG-A-TC')
secuencia2 = ('GGAAATC')
print(secuencia1)
print(secuencia2)
print('score =',hallarScore(secuencia1, secuencia2),"\n")


secuencia1 = ('GGA--TC')
secuencia2 = ('GGAAATC')
print(secuencia1)
print(secuencia2)
print('score =',hallarScore(secuencia1, secuencia2))

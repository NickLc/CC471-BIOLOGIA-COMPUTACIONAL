# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:26:37 2019

@author: ADMIN
"""

secuencia1 = ('ACTGG')
secuencia2 = ('ACTAA')
interseccion = []
x = 0
for i in secuencia1:
    if i == secuencia2[x]:
        interseccion.append('|')
    else:
        interseccion.append(' ')
    x = x + 1
    
print(''.join(secuencia1))
print(''.join(interseccion))
print(''.join(secuencia2))
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:26:37 2019

@author: ADMIN
"""

secuencia1 = ('ACTGGTCAACTGGTCA')
secuencia2 = ('ACTAACTGGTCAATCA')
interseccion = []

print(' ',' '.join(secuencia2))      #Secuencia 2 en horizontal

for i in secuencia1:
    interseccion.append(i)          #Secuencia 1 en vertical
    for j in secuencia2:
        if i == j:
            interseccion.append(' •')
        else:
            interseccion.append('  ')
    print(''.join(interseccion))
    interseccion = []




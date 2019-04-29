# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:23:28 2019

@author: ADMIN
"""

secuencia1 = ('TGAGAGTTGACAGTTGACAGT')
secuencia2 = secuencia1
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




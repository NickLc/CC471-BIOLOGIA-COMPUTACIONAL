# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 23:02:33 2019

@author: ADMIN
"""
# Leer la primera secuencia
direccion = 'secuencia1.txt'
f = open (direccion,'r')
secuencia1 = f.read()
print('Read', direccion)
f.close()

# Leer la segunda secuencia
direccion = 'secuencia2.txt'
f = open (direccion,'r')
secuencia2 = f.read()
print('Read', direccion)
f.close()

# Guardar la dot-plot
direccion = "salida.txt"
file = open(direccion, "w")
interseccion = []

file.write(' '+' '.join(secuencia2)+'\n')

for i in secuencia1:
    interseccion.append(i)          #Secuencia 1 en vertical
    for j in secuencia2:
        if i == j:
            interseccion.append(' •')
        else:
            interseccion.append('  ')
    file.write(''.join(interseccion)+'\n')
    interseccion = []

print('Dot plot creado '+ direccion)
file.close()
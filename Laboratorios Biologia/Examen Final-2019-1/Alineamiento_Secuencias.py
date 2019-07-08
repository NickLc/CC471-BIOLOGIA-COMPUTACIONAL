
# Examen Final CC471-2019-1 
# Problema N° 1b 
# Autor: Lazaro Camasca Edson
# Nombre del Programa: Alineamiento_Secuencias.py

#========================================================================
# 1.b) Alineamiento de secuencias - Generar el archivo clustal
# Utilizando Clustal Omega para Windows

#./clustalo -i globin.fa -o globin.aln --outfmt=clu --force

import os

def alinear_Secuencias(filename = 'Sec_Unidas.fasta'):
    dir_inicio='data_gen'
    dir_final = 'Clustal_Omega'    
    salida_ali = 'Sec_Alineadas.clustal'  #Nombre el archivo que contiene el alineamiento
    print('Alineando las secuencias......')    
    mover_archivo(dir_inicio, dir_final, filename)
    #Realizar el alineamiento
    comando = 'cd Clustal_Omega & clustalo.exe -i {} -o {} --outfmt=clu --force'.format(filename, salida_ali)
    os.system(comando)
    print('Se genero el archivo: {}'.format(salida_ali))

    mover_archivo(dir_final, dir_inicio, salida_ali)
    mover_archivo(dir_final, dir_inicio, filename)

	
#----------------------------------------------------------------------------

def mover_archivo(dir_inicio, dir_final, filename):
    comando = 'move {}\{} {}\{}'.format(dir_inicio,filename, dir_final, filename)
    os.system(comando)

if __name__ == '__main__':
    alinear_Secuencias('protsec.fasta')

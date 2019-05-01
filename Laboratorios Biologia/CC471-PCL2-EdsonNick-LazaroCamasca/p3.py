#==================================================================
# Alumno: Edson Lazaro Camasca
# Codigo: 20160468E
# Curso: Biologia computacional
# 30/04/2019
# p3.py
#==================================================================

from Bio import SeqIO
#Creamos una lista de sequencuas
lista_sequences = [] # Setup an empty list

#Leemos el archivo gb con un parse
for record in SeqIO.parse("PCL2.gb", "genbank"):
    lista_sequences.append(record)

#==========================================================

print("Pregunta 3.1\n")
#Convertir a formato fasta
print("Convetiendo a formato fasta")
SeqIO.write(lista_sequences, "PCL2.fasta", "fasta")


#==========================================================
print("Pregunta 3.2\n")
#Archivo de salida
output_filename = "PCLFILTRO.fasta"
#Secuencia terminal - Tomando 3 registros del ejercicio 2
sec_ter = 'FSL' 
#Lista de las secuencias que cumplan con el filtro
list_filtro = []
count = 0 
total = 0
for record in SeqIO.parse("PCL2.fasta", "fasta"):
    # Longitud de las secuencias
    lon = len(record.seq) 
    #Los ultimos cuatro requistros
    ult_reg = str(record.seq[lon-3:lon])
    #Contar las secuencias que cumplas con la condicion
    total +=1  
    # Aplicando el filtro de que terminen el 'TFS'
    if ult_reg == sec_ter:
        count += 1
        print(ult_reg)
        #Recoremos la lista imprimiento el id y la longitud de la sequencia
        list_filtro.append(record)
        
SeqIO.write(list_filtro, output_filename, "fasta")
print(str(count) + " records selected out of " + str(total))
    
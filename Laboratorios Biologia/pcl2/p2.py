#==================================================================
# Alumno: Edson Lazaro Camasca
# Codigo: 20160468E
# Curso: Biologia computacional
# 30/04/2019
# p2.py
#==================================================================

from Bio import SeqIO
#Creamos una lista de sequencuas
lista_sequences = [] # Setup an empty list

#Leemos el archivo gb con un parse
for record in SeqIO.parse("PCL2.gb", "genbank"):
    lista_sequences.append(record)
#Fichero para escribir
f = open ("texto p2.txt", "w")

f.write("Pregunta 2\n")
f.write("Nombre: Edson Lazaro Camasca\nCodigo: 20160468")

#Secuencia terminal
sec_ter = 'TFSL'
count = 0
total = 0
for sequence in lista_sequences:
    # Longitud de las secuencias
    lon = len(sequence.seq) 
    #Los ultimos cuatro requistros
    ult_reg = str(sequence.seq[lon-4:lon])
    total +=1
    #Contar las secuencias que cumplas con la condicion
    if ult_reg == sec_ter:
        count += 1
        #Recoremos la lista imprimiento el id y la longitud de la sequencia
        f.write("Id: {}\nLongitud: {}\n".format(sequence.id,len(sequence.seq)))
        print("Id: {}\nLongitud: {}\n".format(sequence.id,len(sequence.seq)))
    
    
print(str(count) + " records selected out of " + str(total))

f.close()
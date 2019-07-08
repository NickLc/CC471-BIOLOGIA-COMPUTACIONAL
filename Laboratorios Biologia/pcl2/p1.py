#==================================================================
# Alumno: Edson Lazaro Camasca
# Codigo: 20160468E
# Curso: Biologia computacional
# 30/04/2019
# p1.py
#==================================================================

from Bio import SeqIO
#Creamos una lista de sequencuas
lista_sequences = [] # Setup an empty list

#Leemos el archivo gb con un parse
for record in SeqIO.parse("PCL2.gb", "genbank"):
    lista_sequences.append(record)
#Fichero para escribir
f = open ("texto p1.txt", "w")

f.write("Pregunta 1.1")
f.write("Nombre: Edson Lazaro Camasca\nCodigo: 20160468")

f.write("\nPregunta 1.2")
for sequence in lista_sequences:
	#Recoremos la lista imprimiento el id y la longitud de la sequencia
    f.write("Id: {}\nLongitud: {}\n".format(sequence.id,len(sequence.seq)))

f.write("\nPregunta 1.3")
f.write("Secuencias con longitud mayor a 233\n")
for seq in lista_sequences:
	#imprimimos solo las secuencias con longitud mayor a 233
    longitud = len(seq)
    if longitud > 233:
        f.write("Id: {}\nLongitud: {}\n".format(seq.id,longitud))
        
f.close()
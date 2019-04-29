#importamos el paquete AlignIO
from Bio import AlignIO

#creamos el objeto del alineamiento
alignment = AlignIO.read(open("LHBs_variants.aln"), "clustal")

#Imprimimos la longitud del alineamiento
print "Aligment lenght: ",alignment.get_alignment_length()

print alignment[0]
#print alignment.get_column(38)
#print alignment.get_column(39)
#print alignment.get_column(40)

#Columna 38
print alignment[:, 38]
print alignment[:, 39]
print alignment[:, 40]


print "RESULTADO DE LAS 10 PRIMERAS COLUMNAS"
for i in range(10):
	print "columna ",i,": ",alignment[:,i]

#==================================================================
# Alumno: Edson Lazaro Camasca
# Codigo: 20160468E
# Curso: Biologia computacional
# 30/04/2019
# p4.py
#==================================================================

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
#Biopython ahora sabrá la diferencia entre una 
#base de ADN A de adenina y un residuo de proteina A de alanina.

#Archivo de salida
output_filename = "mRNA.fasta"

#Lista de las secuencias codificantes
list_my_rna = []
for record in SeqIO.parse("PCLFILTRO.fasta", "fasta"):
    secuencia = str(record.seq)
    list_my_rna.append(Seq(secuencia, generic_rna))

for i in list_my_rna:
    print(i)
    
#SeqIO.write(list_my_rna, output_filename, "fasta")
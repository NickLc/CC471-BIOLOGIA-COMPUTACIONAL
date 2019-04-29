from Bio.Align import AlignInfo
from Bio import AlignIO

#creamos el objeto del alineamiento
alignment = AlignIO.read(open("LHBs_variants.aln"), "clustal")

#Imprimimos la longitud del alineamiento
print ("Aligment lenght: ",alignment.get_alignment_length())


#objeto para informacion sumarizada
summary_align = AlignInfo.SummaryInfo(alignment)

#calcular una secuencia de concenso simple

consensus = summary_align.dumb_consensus()

print (consensus)

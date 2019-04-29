from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SubsMat #Diccionario de remplazos

#creamos el objeto del alineamiento
alignment = AlignIO.read(open("LHBs_variants.aln"), "clustal")

#Imprimimos la longitud del alineamiento
print "Aligment lenght: ",alignment.get_alignment_length()


#objeto para informacion sumarizada
summary_align = AlignInfo.SummaryInfo(alignment)

#calcular una secuencia de concenso simple

consensus = summary_align.dumb_consensus()

print consensus

#Encontrar la PSSM 
my_pssm = summary_align.pos_specific_score_matrix(consensus)

print my_pssm[0]["M"]

#Matrix de Sustitucion
replace_info = summary_align.replacement_dictionary()
print replace_info[ ("A", "G")]
print replace_info[ ("A", "K")]

my_arm = SubsMat.SeqMat(replace_info)
my_lom = SubsMat.make_log_odds_matrix(my_arm)
my_lom.print_full_mat()



#======================================================
# Unir secuencias en formato fasta 
from Bio import SeqIO
import os

def join_Sqs(address_dir, address_output):
	"""Junta todas las secuencias que se encuentran en address_dir y los retorna en address_output"""
	records = []  
	for filename in os.listdir(address_dir):
		handle = open(address_dir + "/" + filename)
		record = SeqIO.read( handle, "swiss")
		records.append(record)

	print ("Numero de registros:", len(records))
	SeqIO.write(records, address_output, "fasta")

#======================================================
# Leer una secuencia sea fasta o gebank
from Bio import SeqIO
def read_Sq(address_sq, tipo = 'fasta'):
	"""Leer secuencias sea fasta o gebank"""
	
	if tipo == 'fasta':
		record = SeqIO.read(address_sq, "fasta")
	if tipo == 'genbank':
		record = SeqIO.read(address_sq, "genbank")

	print("Secuencia {} leida con codigo {}".format(address_sq, record.id))
	#print ("Codigo: ",record.id)
	#print ("Descricion:", record.description)
	#print ("Nombre: ",record.name)
	#print ("Secuencia:", record.seq)
	return record
	# file/P17102.fasta
	# file/P17102.gbk

#======================================================
# Leer secuencia alineada

from Bio import AlignIO
def read_SeqAli(address_Seq):
	#creamos el objeto del alineamiento
	alignment = AlignIO.read(open(address_Seq), "clustal")
	print("Secuencias alineadas leidas")
	#Imprimimos la longitud del alineamiento
	#print("Aligment lenght: ",alignment.get_alignment_length())
	return alignment

#======================================================

def get_columnSecAli(alignment,pos = -1, pos_init = 0, pos_end=0):
	"""Optiene las columnas del Alineamiento de secuencias
	pos si se quiere una columna en especifico y [init,end] si se quiere un bloque"""
	# Una sola columna 
	if pos != -1:
		column = alignment[:, pos]
		return column
	#Varias columnas
	else:
		column = [alignment[:,i] for i in range(pos_init, pos_end)]
		return column

#======================================================

from Bio.Align import AlignInfo
def get_sumatyInfo(alignment):
	#objeto para informacion sumarizada
	summary_align = AlignInfo.SummaryInfo(alignment)
	return summary_align

#======================================================
def get_consesus(summary_align):
	"""Una secuencia de consenso, es una secuencia ideal que representa los nucleótidos o 
	aminoácidos que se encuentran con mayor frecuencia en cada posición de un fragmento de 
	DNA o de una proteína, respectivamente"""

	#calcular una secuencia de concenso simple
	consensus = summary_align.dumb_consensus()
	#dumb_consensus() que sirve para calcular una secuencia de consenso simple.
	return(consensus)

#======================================================

from Bio import SubsMat #Diccionario de remplazos

# Matrices de Score de Posiciones Especificas (PSSMs)
def get_PSSM(summary_align):
	"""Un PSSM, es un tipo de matriz de puntuación utilizada en las búsquedas de proteínas BLAST en
	las que las puntuaciones de sustitución de aminoácidos se dan por separado para cada posición en
	una alineación de secuencia múltiple de proteínas"""
	#calcular una secuencia de concenso simple
	consensus = summary_align.dumb_consensus()
	#Encontrar la PSSM 
	pssm = summary_align.pos_specific_score_matrix(consensus)
	return pssm
#======================================================

def get_element_of_PSSM(pssm, num_sq, name_res):
	"""Se accede a un elemento de la PSSM con el numero de la secuencia y el nombre del residuo"""
	return pssm[num_sq][name_res]

#======================================================

from Bio import SubsMat
def create_matriz_Sustitucion(summary_align):
	"""Esta información nos da nuestro número aceptado de reemplazos, o con qué frecuencia esperamos
	que diferentes residuos se sustituyan entre sí"""
	replace_info = summary_align.replacement_dictionary()
	#print replace_info[ ("A", "G")]
	#print replace_info[ ("A", "K")]
	#La funcion SeqMat() toma como parametro el diccionario de reemplazos
	my_arm = SubsMat.SeqMat(replace_info)
	#crear una matriz de reemplazo aceptada (Accepted Replacement Matrix - ARM).
	my_lom = SubsMat.make_log_odds_matrix(my_arm)
	my_lom.print_full_mat()
	return my_lom
#======================================================

if __name__ == 'Alineamiento_Secuencias':
	print("Alineamiento_Secuencicas se ha importado correctamente.")

if __name__ == '__main__':
	join_Sqs("data/LHBs", "data/LHBs_variants.fasta")
	seq = read_Sq('data/P17102.fasta', 'fasta')
	print(seq)
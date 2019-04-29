# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:55:33 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO
wanted_ids = ["PGSC0003DMP400019313", "PGSC0003DMP400020381", "PGSC0003DMP400020972"]
fileName = 'PGSC_DM_v3.4_pep_representative.fasta'
output = 'wanted_potato_proteins_in_order.fasta'
fasta_index = SeqIO.index(fileName, 'fasta')
total = len(fasta_index)
print (str(total) + " secuencias en "+fileName)
cont = 0
output_handle = open(output, 'w')

for record_id in wanted_ids:
    cont = cont + 1
    record = fasta_index[record_id]
    SeqIO.write(record, output_handle, 'fasta')

output_handle.close()
print(str(cont)+" secuencias selecionadas de "+str(total))

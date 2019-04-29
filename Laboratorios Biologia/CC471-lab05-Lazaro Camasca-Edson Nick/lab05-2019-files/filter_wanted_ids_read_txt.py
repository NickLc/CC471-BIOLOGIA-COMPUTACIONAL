# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:46:35 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO
wanted_ids = []
 
#ingresar "wanted_ids.txt"
#f = open(input(str("Ingrese la direccion del fichero wanted_ids: ")))
f = open("wanted_ids.txt")
for linea in f:
   wanted_ids.append(f.readline())
   
f.close()

input_filename = "PGSC_DM_v3.4_pep_representative.fasta"
output_filename = "wanted_potato_proteins.fasta"
count = 0
total = 0
output_handle = open(output_filename, "w")
# ...
# Your code here
for record in SeqIO.parse(input_filename, "fasta"):
    total +=1
    # Encuentra la sequencias queridas 
    if (record.id in wanted_ids):
        SeqIO.write(record, output_handle, "fasta")
        count +=1
# ...       
output_handle.close()
print(str(count) + " records selected out of " + str(total))
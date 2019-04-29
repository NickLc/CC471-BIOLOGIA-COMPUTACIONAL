# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:35:37 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO
wanted_ids = ["PGSC0003DMP400019313", "PGSC0003DMP400020381", "PGSC0003DMP400020972"]
input_filename = "PGSC_DM_v3.4_pep_representative.fasta"
output_filename = "wanted_potato_proteins.fasta"
count = 0
total = 0
output_handle = open(output_filename, "w")
# ...
# Your code here
for record in SeqIO.parse(input_filename, "fasta"):
    total +=1
    print(record.id)
    # Encuentra la sequencias queridas 
    if (record.id in wanted_ids):
        SeqIO.write(record, output_handle, "fasta")
        count +=1
# ...       
output_handle.close()
print(str(count) + " records selected out of " + str(total))
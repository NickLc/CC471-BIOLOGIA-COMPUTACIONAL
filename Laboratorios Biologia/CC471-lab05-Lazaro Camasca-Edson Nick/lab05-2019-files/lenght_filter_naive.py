# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:03:56 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO
input_filename = "NC_000913.faa"
output_filename = "NC_000913_long_only.faa"
count = 0
total = 0
for record in SeqIO.parse(input_filename, "fasta"):
    total = total + 1
    if 100 <= len(record):
        count = count + 1
        SeqIO.write(record, output_filename, "fasta")
#print("Contenido\n"+record+"\n")
print(str(count) + " records selected out of " + str(total))


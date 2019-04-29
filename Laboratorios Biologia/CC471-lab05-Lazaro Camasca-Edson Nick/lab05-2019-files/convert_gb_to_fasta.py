# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:50:23 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO

print("Convert gt to fasta")
#input_filename = "NC_000913.gbk"
input_filename = input(str("input_filename: "))
#output_filename = "NC_000913_converted.fasta"
output_filename = input(str("output_filename: "))
count = SeqIO.convert(input_filename, "gb", output_filename, "fasta")
print(str(count) + " records converted")
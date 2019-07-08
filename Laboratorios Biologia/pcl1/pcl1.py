# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 09:12:50 2019
@title: pcl1.py
@author: Edson Lazaro Camasca
"""

from Bio import SeqIO
import os

records = []
for filename in os.listdir("secfasta"):
	handle = open("secfasta" + "/" + filename)
	record = SeqIO.read( handle, "fasta")
	records.append(record)

print ("Numero de registros:", len(records))

SeqIO.write(records, "totalfasta.fasta", "fasta")
print("totalfasta es:", records)

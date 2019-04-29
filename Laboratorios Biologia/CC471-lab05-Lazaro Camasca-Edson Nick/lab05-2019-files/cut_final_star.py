# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:19:27 2019
@course: Biologia Computacional
@author: Lazaro Camasca Edson Nick
"""
from Bio import SeqIO
input_filename = "PGSC_DM_v3.4_pep_representative.fasta"
output_filename = "PGSC_DM_v3.4_pep_rep_no_stars.fasta"
output_handle = open(output_filename, "w")
for record in SeqIO.parse(input_filename, "fasta"):
    if (record[-1] == '*'):     #Solo los que tengan un * al final
        cut_record = record[:-1] # remove last letter
        SeqIO.write(cut_record, output_handle, "fasta")
output_handle.close()
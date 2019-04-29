from Bio import SeqIO
import os

records = []
for filename in os.listdir("file/LHBs"):
	handle = open("file/LHBs" + "/" + filename)
	record = SeqIO.read( handle, "swiss")
	records.append(record)

print ("Numero de registros:", len(records))

SeqIO.write(records, "LHBs_variants.fasta", "fasta")


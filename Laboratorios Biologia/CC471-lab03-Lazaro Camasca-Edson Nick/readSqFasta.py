from Bio import SeqIO

record = SeqIO.read("file/P17102.fasta", "fasta")
print ("Descripcion:", record.description)
print ("Secuencia:", record.seq)

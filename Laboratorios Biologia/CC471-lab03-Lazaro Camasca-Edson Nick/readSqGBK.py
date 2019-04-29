from Bio import SeqIO

record = SeqIO.read("file/P17102.gbk", "genbank")

print ("Codigo: ",record.id)
print ("Descricion:", record.description)
print ("Nombre: ",record.name)
print ("Secuencia:", record.seq)

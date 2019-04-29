from Bio import SeqIO
#P7
filename = "PGSC_DM_v3.4_pep_representative.fasta"
bad = 0
final_stop = 0
for record in SeqIO.parse(filename, "fasta"):
	if record.seq[-1] == '*':
		final_stop +=1 
		
	for i in record.seq:
		if i == '*':
			bad +=1 
			#print(record.description+" starts " + record.seq[0])
			break
			
print("Secuencias que tienen un *: "+str(bad))
print("Secuencias que tienen un * al final: "+str(final_stop))

from Bio import SeqIO
import pylab   # lib para realizar histograma

filename = "NC_000913.faa"
count_lenght_100 = 0
lenght_all_seq = 0
sizes = []
for record in SeqIO.parse(filename, "fasta"):
	#print("Record " + record.id + ", length " + str(len(record.seq)))
# Pregunta_2
	if len(record.seq)<100:
		count_lenght_100 +=1
# Pregunta_3
	lenght_all_seq += len(record.seq)
# Pregunta_4
	sizes.append(len(record.seq))
	
print("P2.Cantidad de secuencias que tienen una longitud menor a 100: "+ str(count_lenght_100))

print("P3.Longitud total de todas las secuencias:"+str(lenght_all_seq))

print("P4.El histograma es:")
pylab.hist(sizes, bins=50)
pylab.title("%i orchid sequences\nLengths %i to %i" \
            % (len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Sequence length (bp)")
pylab.ylabel("Count")
pylab.show()


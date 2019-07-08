#======================================================================================
# Count_fasta
from Bio import SeqIO
import sys

#filename = sys.argv[1]

def count_fasta_registro(filename):
    """Cuenta los registros de un archivo fasta
    filename es la direccion del archivo fasta"""
    count = 0
    for record in SeqIO.parse(filename, "fasta"):
        count = count + 1
    print("Existen " + str(count) + " registros en el archivo " + filename)

#======================================================================================
# record_lenght

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

#======================================================================================
# Check_start_met.py

from Bio import SeqIO
#P5
filename = "PGSC_DM_v3.4_pep_representative.fasta"
bad = 0
for record in SeqIO.parse(filename, "fasta"):
	if not record.seq.startswith("M"):
		bad = bad + 1
		#P6 
		print(record.description+" starts " + record.seq[0])
		#print(record.id + " starts " + record.seq[0])
	
print("Found " + str(bad) + " records in " + filename + " which did not start with M")

#======================================================================================
# Check_stop.py

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

#======================================================================================
# P8
from Bio import SeqIO

#FORMATO FASTA
fasta_record = SeqIO.read("NC_000913.fna", "fasta")
print("Formato fasta: "+fasta_record.id + " length " + str(len(fasta_record)))

#FORMATO GBANK
genbank_record = SeqIO.read("NC_000913.gbk", "genbank")
print("formato gebank: "+genbank_record.id + " length " + str(len(genbank_record)))

#======================================================================================

# convert_gb_to_fasta
#P9
from Bio import SeqIO

print("Convert gt to fasta")
#input_filename = "NC_000913.gbk"
input_filename = input(str("input_filename: "))
#output_filename = "NC_000913_converted.fasta"
output_filename = input(str("output_filename: "))
count = SeqIO.convert(input_filename, "gb", output_filename, "fasta")
print(str(count) + " records converted")

#======================================================================================

#lenght_filter_naive
# P10
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

#======================================================================================

#P10_2
#lenght_filter_corregido
from Bio import SeqIO
input_filename = "NC_000913.faa"
output_filename = "NC_000913_long_only.faa"
count = 0
total = 0
output_handle = open(output_filename, "w")
for record in SeqIO.parse(input_filename, "fasta"):
    total = total + 1
    if 100 <= len(record):
        count = count + 1
        SeqIO.write(record, output_handle, "fasta")
output_handle.close()
print(str(count) + " records selected out of " + str(total))

#======================================================================================
# P12
# Cut_final_start
from Bio import SeqIO
input_filename = "PGSC_DM_v3.4_pep_representative.fasta"
output_filename = "PGSC_DM_v3.4_pep_rep_no_stars.fasta"
output_handle = open(output_filename, "w")
for record in SeqIO.parse(input_filename, "fasta"):
    if (record[-1] == '*'):     #Solo los que tengan un * al final
        cut_record = record[:-1] # remove last letter
        SeqIO.write(cut_record, output_handle, "fasta")
output_handle.close()
#======================================================================================
# P13
#filter_wanted_ids
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

#======================================================================================
# P14
#filter_wanted_ids_read_txt

from Bio import SeqIO
wanted_ids = []
 
#ingresar "wanted_ids.txt"
#f = open(input(str("Ingrese la direccion del fichero wanted_ids: ")))
f = open("wanted_ids.txt")
for linea in f:
   wanted_ids.append(f.readline())
   
f.close()

input_filename = "PGSC_DM_v3.4_pep_representative.fasta"
output_filename = "wanted_potato_proteins.fasta"
count = 0
total = 0
output_handle = open(output_filename, "w")
# ...
# Your code here
for record in SeqIO.parse(input_filename, "fasta"):
    total +=1
    # Encuentra la sequencias queridas 
    if (record.id in wanted_ids):
        SeqIO.write(record, output_handle, "fasta")
        count +=1
# ...       
output_handle.close()
print(str(count) + " records selected out of " + str(total))

#======================================================================================
# P15
#filter_wanted_id_in_order

from Bio import SeqIO
wanted_ids = ["PGSC0003DMP400019313", "PGSC0003DMP400020381", "PGSC0003DMP400020972"]
fileName = 'PGSC_DM_v3.4_pep_representative.fasta'
output = 'wanted_potato_proteins_in_order.fasta'
fasta_index = SeqIO.index(fileName, 'fasta')
total = len(fasta_index)
print (str(total) + " secuencias en "+fileName)
cont = 0
output_handle = open(output, 'w')

for record_id in wanted_ids:
    cont = cont + 1
    record = fasta_index[record_id]
    SeqIO.write(record, output_handle, 'fasta')

output_handle.close()
print(str(cont)+" secuencias selecionadas de "+str(total))

#======================================================================================




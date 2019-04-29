from Bio import SeqIO

#FORMATO FASTA
fasta_record = SeqIO.read("NC_000913.fna", "fasta")
print("Formato fasta: "+fasta_record.id + " length " + str(len(fasta_record)))

#FORMATO GBANK
genbank_record = SeqIO.read("NC_000913.gbk", "genbank")
print("formato gebank: "+genbank_record.id + " length " + str(len(genbank_record)))
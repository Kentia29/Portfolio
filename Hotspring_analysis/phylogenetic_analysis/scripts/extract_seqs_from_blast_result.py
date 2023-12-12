import sys
from Bio import SeqIO

# read a BLAST output file, line by line, and save the name of the
# target queries. 
homologs = set()
for line in open(sys.argv[1]):
    # a list of values for each hit
    fields = list(map(str.strip, line.split('\t')))
    homologs.add(fields[1])

# Extract the sequence of each homolog and print it
for r in SeqIO.parse(sys.argv[2], format="fasta"):
    if r.id in homologs: 
        print ('>%s\n%s' %(r.id, r.seq))

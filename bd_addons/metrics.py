from Bio import SeqIO
import os

cur = os.getcwd()

swissprot = list(SeqIO.parse("./uniprot_sprot.fasta", "fasta"))
#!/bin/bash

# For this step you will need:
#- Blast+ in folder "binx"
#	NCBI ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 
#or 
#	wget -P ../binx/ ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
#	tar -xf ../binx/ncbi-blast-2.11.0+-x64-linux.tar.gz -C ../binx/

#-SwissProt Database in folder "data"
#	wget -P ../data/msa/ ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#	gunzip ../data/msa/uniprot_sprot.fasta.gz

echo "Choose MSA Method: (C)offee, Clustal-(O) or (M)uscle."
read msamethod
 
echo "Enter Try number"
read try

#create PSSM

../../../binx/ncbi-blast-2.11.0+/bin/psiblast -db ../../../data/uniprot_sprot.fasta -in_msa ../MSAs/MSA_${msamethod}_${try}.fa -out_ascii_pssm pssm_${msamethod}_${try}.txt -out out_psiblast_${msamethod}_${try}.txt
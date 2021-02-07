This folder contains the hmm models generated from our MSAs. For each MSA try we generate:
 - A HMM model;
 - The output of our HMM search on Swissprot, using such HMM as input. The output is in three different formats: align, tblout, domtblout.
 
 Steps to create HMMs and do the search (Works only for Linux I think):
# 1 Download HMMER and compile in folder of choice
wget -P ./<foldername>/ http://eddylab.org/software/hmmer/hmmer-3.3.1.tar.gz
tar -xf ./<foldername>/hmmer-3.3.1.tar.gz -C ../binx/
cd ./<foldername>/hmmer-3.3.1/ && ./configure && make

# Build a HMM from the MSA with hmmbuild
./<foldername>/hmmer-3.3.1/src/hmmbuild <output .hmm file path> <MSA .fa file path>
 
# Command line commands to search SwissProt, using the created HMM model as query. This commands will output all 3 different formats of output.
# NB: you will need to have the Swissprot database downloaded

../binx/hmmer-3.3.1/src/hmmsearch --tblout <tblout out file path> --domtblout <domtblout out file path> <query hmm model path> <path of database> > <align out file path>


To download Swissprot:

wget -P ../data/msa/ ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip ../data/msa/uniprot_sprot.fasta.gz


------ Naming Conventions -----

- HMM MODEL: hmm_model_X_Y;
- ALIGN OUTPUT: hmmsearch_out_X_Y.align;
- TBLOUT: hmm_search_out_X_Y.tblout;
- DOMTBLOUT: hmm_search_out_X_Y.domtblout;

For all it must be X = MSA method, Y = try number





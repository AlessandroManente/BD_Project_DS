# Here we take care of the mapping from the sifts table
from Bio import SeqIO
import pandas as pd
import os

def generate_pdb_df(sifts_path, model_output_path):
    model_prots_df = pd.read_csv(model_output_path)
    model_prots = list(model_prots_df.ids.values)

    # Import pdb - uniprot relation file
    sifts_tsv = pd.read_csv(sifts_path, sep = '\t', header = 1)
    sifts_tsv.columns = list(map(lambda x: x.lower(), sifts_tsv.columns.values))
    pdb_df = sifts_tsv.loc[sifts_tsv.sp_primary.isin(model_prots),['pdb','sp_primary','chain','pdb_beg', 'pdb_end','sp_beg','sp_end']].copy()
    pdb_df = pdb_df.reset_index()

    return pdb_df

def filter_pdb_db(pdb_db):
    denom = (pdb_db.sp_end - pdb_db.sp_beg+1).astype(int)
    num = (pdb_db.pdb_end.astype(int) - pdb_db.pdb_beg.astype(int)+1)
    coverage = num/denom
    pdb_db['coverage'] = coverage
    pdb_db = pdb_db[pdb_db.coverage > 0.8]

    return pdb_db
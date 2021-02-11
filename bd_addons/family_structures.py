from Bio import SeqIO
import pandas as pd
import os

def generate_pdb_dfs(sifts_path, model_output_path):
    model_prots_df = pd.read_csv(model_output_path)
    model_prots = list(model_prots_df.ids.values)

    # Import pdb - uniprot relation file
    sifts_tsv = pd.read_csv(sifts_path, sep = '\t', header = 1)
    sifts_tsv.columns = list(map(lambda x: x.lower(), sifts_tsv.columns.values))
    pdb_df = sifts_tsv.loc[sifts_tsv.sp_primary.isin(model_prots),['pdb','sp_primary','chain','sp_beg','sp_end']].copy()
    model_start = list(pdb_df.sp_primary.map(lambda x: model_prots_df[model_prots_df.ids == x].hit_start.values[0]))
    model_end = list(pdb_df.sp_primary.map(lambda x: model_prots_df[model_prots_df.ids == x].hit_end.values[0]))
    pdb_df['model_start'] = model_start
    pdb_df['model_end'] = model_end
    return pdb_df


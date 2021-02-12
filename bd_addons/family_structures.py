# Here we take care of the mapping from the sifts table
from Bio import SeqIO
import pandas as pd
import os
import itertools

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

def parse_tmalign_out(filename):
    """
    Takes as input the path of the .out file provided by TM-Align.
    Returns both the RMSD and the TMSCORE (the first one provided in the .out file) """
    lines = []
    temp_path = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\temp"
    with open(temp_path + '\\' + filename, 'r') as f:
        for line in f:
            for line in itertools.islice(f, 15,21):
                lines.append(line)
    RMSD = lines[0].split(',')[1].split('=')[1].strip()
    TMSCORE = lines[2].split('=')[1].split('(')[0].strip()

    return float(RMSD), float(TMSCORE)

def create_rmsd_matrix(best_model):
    """
    Takes all the files inside the temp folder, reads the rmsd values from them, 
    and stores them into a matrix"""

    model_path = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\pdbs_{}".format(best_model)

    if 'rmsds_{}.csv'.format(best_model) in os.listdir(model_path):
        rmsds_df = pd.read_csv(model_path + '\\' + 'rmsds_' + best_model + '.csv', index_col=0)
        return rmsds_df

    else:
        rmsd_dict = {}

        dir_to_parse = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\temp"
        files = os.listdir(dir_to_parse)
        for filename in files:
            ls = filename.split("_")
            o1 = ls[0].split(".")[0]
            o2 = ls[1].split(".")[0]
            rmsd_dict.setdefault(o1, {})
            rmsd_dict[o1][o2], _ = parse_tmalign_out(filename) 

        rmsd_df = pd.DataFrame.from_dict(rmsd_dict)
        return rmsd_df

def create_tmscores_matrix(best_model):
    """
    Takes all the files inside the temp folder, reads the tmscore values from them, 
    and stores them into a matrix"""

    model_path = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\pdbs_{}".format(best_model)

    if 'tmscores_{}.csv'.format(best_model) in os.listdir(model_path):
        tmscore_df = pd.read_csv(model_path + '\\' + 'tmscores_' + best_model + '.csv', index_col=0)
        return tmscore_df
    else:

        tmscore_dict = {}

        dir_to_parse = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\temp"
        files = os.listdir(dir_to_parse)
        for filename in files:
            ls = filename.split("_")
            o1 = ls[0].split(".")[0]
            o2 = ls[1].split(".")[0]
            tmscore_dict.setdefault(o1, {})
            _, tmscore_dict[o1][o2] = parse_tmalign_out(filename) 

        tmscore_df = pd.DataFrame.from_dict(tmscore_dict)
        tmscore_df.to_csv(model_path + '\\' + best_model + '.csv')
        return tmscore_df

def clear_temp_folder():
    temp_dir = ".\\data_team_1\\_part_2\\original_datasets\\family_structures\\temp"
    files = os.listdir(temp_dir)
    for filename in files:
        if filename.split('.')[-1] == 'out':
            os.remove(temp_dir + '\\' + filename)
    return
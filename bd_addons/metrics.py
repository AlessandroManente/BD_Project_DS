from Bio import SeqIO, SearchIO
import pandas as pd
import os
import numpy as np
from bd_addons.HmmPy import *

cur = os.getcwd()
path_swiss_prot = cur + '\\data_team_1\\swiss_prot\\uniprot_sprot.fasta'


def swiss_prot_parser():
    if 'uniprot_sprot.csv' not in os.listdir(cur + '\\data_team_1\\swiss_prot\\'):
        swissprot = list(SeqIO.parse(path_swiss_prot, "fasta"))
        swissprot = [[str(x.seq), x.id, x.description] for x in swissprot]

        columns_swissprot = ['seq', 'id', 'description']
        swissprot_df = pd.DataFrame(swissprot, columns=columns_swissprot)
        swissprot_df['ids'] = swissprot_df['id'].apply(lambda x: x.replace('|', '-').split('-')[1])
        swissprot_df.to_csv(cur + '\\data_team_1\\swiss_prot\\' + 'uniprot_sprot.csv')
    
    else:
        swissprot_df = pd.read_csv(cur + '\\data_team_1\\swiss_prot\\' + 'uniprot_sprot.csv')
    
    return swissprot_df


def parse_psiblast():
    dirs_to_parse = [cur + '\\data_team_1\\PSSMs\\PSSM_' + a + '\\to_parse' for a in ['C', 'M', 'O']]
    dirs_parsed = [cur + '\\data_team_1\\PSSMs\\PSSM_' + a + '\\parsed' for a in ['C', 'M', 'O']]
    
    parsed_dfs = {}

    for i, dir in enumerate(dirs_to_parse):
        files = os.listdir(dir)
        for filename in files:
            filename = filename.split('.')
            
            parsed_dfs['.'.join(filename)] = psiblast_parser(dir + '\\', filename[0], filename[1], dirs_parsed[i])
    
    return parsed_dfs


def psiblast_parser(dir_to_parse, filename, extension, dir_parsed):   
    if filename + '.csv' not in os.listdir(dir_parsed): 
        print("Parsing {} ...".format(filename))
        hits_dict = {}
        blast_records = SearchIO.parse(dir_to_parse + '\\' + filename + '.' + extension, 'blast-xml')
        for blast_record in blast_records:
            for rec in blast_record.hits:
                hits_dict.setdefault("ids",[]).append(rec.id.split('|')[1])
                hits_dict.setdefault("e_value",[]).append(rec.hsps[0].evalue)
                hits_dict.setdefault("hit_start",[]).append(rec.hsps[0].hit_range[0]+1)
                hits_dict.setdefault("hit_end",[]).append(rec.hsps[0].hit_range[1])
                hits_dict.setdefault("bitscore",[]).append(rec.hsps[0].bitscore)

        lengths = []
        with open(dir_to_parse + '\\' + filename + '.' + extension, 'r') as f:
            for line in f:
                if(line[2:11] == '<Hit_len>'):
                    lengths.append(line.split('>')[1].split('<')[0])
        hits_dict["protein_length"] = lengths
            
        hits_df = pd.DataFrame.from_dict(hits_dict)
        hits_df.to_csv(dir_parsed + '\\' + filename + '.csv')
    
    else:
        hits_df = pd.read_csv(dir_parsed + '\\' + filename + '.csv')
    
    return hits_df
    

def metrics_sequences(df, gt):
    gt_acc = gt.accession.to_list()
    df = df.drop_duplicates(subset=['ids'])
    df_ids = df.ids.to_list()
    swissprot_df = swiss_prot_parser()
    

    true_positives = df['ids'].apply(lambda x: 1 if x in gt_acc else 0).sum() # intersection between gt.ids and df.ids
    true_negatives = swissprot_df['ids'].apply(lambda x: 1 if (x not in gt_acc and x not in df_ids) else 0).sum() # rows in swissprot_df but not in gt.ids and df.ids
    false_positives = len(df) - true_positives
    false_negatives = len(swissprot_df) - len(df) - true_negatives#swissprot_df['ids'].apply(lambda x: 1 if x not in df_ids else 0).sum() - true_negatives
    
    
    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
    
    precision = true_positives / (true_positives + false_positives)

    recall = true_positives / (true_positives + false_negatives)

    specificity = true_negatives / (true_negatives + false_positives)

    balanced_accuracy = (recall + specificity) / 2

    mcc = ((true_positives * true_negatives) - (false_positives * false_negatives)) / np.sqrt((true_positives + false_positives) * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))

    f1_score = (2 * true_positives) / (2 * true_positives + false_positives + false_negatives)

    return [accuracy, precision, recall, specificity, balanced_accuracy, mcc, f1_score]


def metrics_8(gt):
    parsed_tblouts, parsed_domtblouts = parse_hmms()
    parsed_psiblast = parse_psiblast()

    metrics = []
    index_metrics = list(parsed_domtblouts.keys()) + list(parsed_psiblast.keys())
    columns_metrics = ['accuracy', 'precision', 'recall', 'specificity', 'balanced_accuracy', 'mcc', 'f1_score']

    for df in parsed_domtblouts.keys():
        metrics.append(metrics_sequences(parsed_domtblouts[df], gt))
    
    for df in parsed_psiblast.keys():
        metrics.append(metrics_sequences(parsed_psiblast[df], gt))
    
    metrics_df = pd.DataFrame(metrics, index_metrics, columns_metrics)
    metrics_df.to_csv(cur + '\\data_team_1\\metrics\\metrics_8.csv')

    return metrics_df, parsed_tblouts, parsed_domtblouts, parsed_psiblast

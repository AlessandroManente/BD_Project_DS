from Bio import SeqIO, SearchIO
import pandas as pd
import os
import numpy as np

cur = os.getcwd()
path_swiss_prot = cur + '\\data_team_1\\swiss_prot\\uniprot_sprot.fasta'


def swiss_prot_parser():
    if 'uniprot_sprot.csv' not in os.listdir(cur + '\\data_team_1\\swiss_prot\\'):
        swissprot = list(SeqIO.parse(path_swiss_prot, "fasta"))
        swissprot = [[str(x.seq), x.id, x.description] for x in swissprot]

        columns_swissprot = ['seq', 'id', 'description']
        swissprot_df = pd.DataFrame(swissprot, columns=columns_swissprot)
        
        swissprot_df.to_csv(cur + '\\data_team_1\\swiss_prot\\' + 'uniprot_sprot.csv')
    
    else:
        swissprot_df = pd.read_csv(cur + '\\data_team_1\\swiss_prot\\' + 'uniprot_sprot.csv')
    
    return swissprot_df


def psiblast_parser():
    iterations = 2
    xml_path = ".\data_team_1\PSSMs\PSSM_C\out_psiblast_C_1_swissprot_{}iterations.xml".format(iterations)
    
    hits_dict = {}
    blast_records = SearchIO.parse(xml_path, 'blast-xml')
    for blast_record in blast_records:
        for rec in blast_record.hits:
            hits_dict.setdefault("ids",[]).append(rec.id.split('|')[1])
            hits_dict.setdefault("e_value",[]).append(rec.hsps[0].evalue)
            hits_dict.setdefault("hit_start",[]).append(rec.hsps[0].hit_range[0])
            hits_dict.setdefault("hit_end",[]).append(rec.hsps[0].hit_range[1])
            hits_dict.setdefault("bitscore",[]).append(rec.hsps[0].bitscore)
        
    hits_df = pd.DataFrame.from_dict(hits_dict)
    
    return hits_df


def metrics_sequences(df, gt):
    gt_acc = gt.accession.to_list()
    swissprot_df = swiss_prot_parser()

    true_positives = df['ids'].apply(lambda x: split('-', x.replace('|', '-'))[1]).apply(lambda x: 1 if x in gt_acc else 0).sum()
    true_negatives = len(swissprot_df) - len(df)
    false_positives = len(df) - true_positives
    false_negatives = df['id'].apply(lambda x: split('-', x.replace('|', '-'))[1]).apply(lambda x: 1 if x in gt_acc else 0).sum() - true_positives
    
    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
    precision = true_positives / (true_positives + false_positives)
    recall = true_positives / (true_positives + false_negatives)
    balanced_accuracy = (recall + (true_negatives / (true_negatives + false_positives))) / 2
    mcc = ((true_positives * true_negatives) - (false_positives * false_negatives)) / np.sqrt((true_positives + false_positives) * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))
    f1_score = (2 * true_positives) / (2 * true_positives + false_positives + false_negatives)

    return [accuracy, precision, recall, balanced_accuracy, mcc, f1_score]


def metrics_8(dfs_list, gt):
    metrics = []
    index_metrics = ['df_psi_C', 'df_psi_M', 'df_psi_O', 'df_hmm_C', 'df_hmm_M', 'df_hmm_O']
    columns_metrics = ['accuracy', 'precision', 'recall', 'balanced_accuracy', 'mcc', 'f1_score']

    for df in dfs_list:
        metrics.append(df, gt)
    
    metrics_df = pd.DataFrame(metrics, index_metrics, columns_metrics)

    return metrics_df

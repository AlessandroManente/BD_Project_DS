from Bio import SeqIO, SearchIO
import pandas as pd
import os
from bd_addons.HmmPy import *
import time
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import random
sns.set_theme()
sns.set_style("whitegrid")
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


def metrics_8(gt, smart_update = True, force_recompute = False):
    parsed_tblouts, parsed_domtblouts = parse_hmms()
    parsed_psiblast = parse_psiblast()

    metrics = []
    index_metrics = list(parsed_domtblouts.keys()) + list(parsed_psiblast.keys())
    columns_metrics = ['accuracy', 'precision', 'recall', 'specificity', 'balanced_accuracy', 'mcc', 'f1_score']
    
    
    if 'metrics_8.csv' in os.listdir(cur + '\\data_team_1\\metrics'):
        old_metrics_df = pd.read_csv(cur + '\\data_team_1\\metrics\\metrics_8.csv', index_col=0)

        # Check if we have already all the statistics that we need in the the metrics8.csv
        if (old_metrics_df.index.to_list() == index_metrics) and not force_recompute:
            return old_metrics_df, parsed_tblouts, parsed_domtblouts, parsed_psiblast
        
        else:
            for df in parsed_domtblouts.keys():
                # recompute metrics only if the key wasnt already in the old csv!
                if smart_update:
                    if (df not in old_metrics_df.index.to_list()):
                        print("Computing metrics for: {}".format(df))
                        metrics.append(metrics_sequences(parsed_domtblouts[df], gt))
                else:
                    print("Computing metrics for: {}".format(df))
                    metrics.append(metrics_sequences(parsed_domtblouts[df], gt))

            
            for df in parsed_psiblast.keys():
                # recompute metrics only if the key wasnt already in the old csv!
                if (df not in old_metrics_df.index.to_list()) and  smart_update:
                    metrics.append(metrics_sequences(parsed_psiblast[df], gt))
            
            metrics_df = pd.DataFrame(metrics, index_metrics, columns_metrics)
            metrics_df.to_csv(cur + '\\data_team_1\\metrics\\metrics_8.csv')

            return metrics_df, parsed_tblouts, parsed_domtblouts, parsed_psiblast
    
    else:
        for df in parsed_domtblouts.keys():
            print("Computing metrics for: {}".format(df))
            metrics.append(metrics_sequences(parsed_domtblouts[df], gt))
        
        # qua c'era l'indentazione sbagliata?
        for df in parsed_psiblast.keys():
            print("Computing metrics for: {}".format(df))
            metrics.append(metrics_sequences(parsed_psiblast[df], gt))
            
        metrics_df = pd.DataFrame(metrics, index_metrics, columns_metrics)
        metrics_df.to_csv(cur + '\\data_team_1\\metrics\\metrics_8.csv')

        return metrics_df, parsed_tblouts, parsed_domtblouts, parsed_psiblast


def dfs_to_dicts(parsed_domtblouts, parsed_psiblast):
    # from dfs to dicts with no repeated columns -> fast af
    list_dfs = list(parsed_domtblouts.keys()) + list(parsed_psiblast.keys())
    list_dfs_not_repeated = []

    for df in parsed_domtblouts.keys():
        df_not_repeated = {}
        for i, row in parsed_domtblouts[df].iterrows():
            if 'domtblout' in df:
                df_not_repeated.setdefault(row['ids'],[]).append([row['from_ali_coord'], row['to_ali_coord'], row['target_len']])
            else:
                df_not_repeated.setdefault(row['ids'],[]).append([row['hit_start'], row['hit_end'], row['protein_length']])
        
        list_dfs_not_repeated.append(df_not_repeated)

    for df in parsed_psiblast.keys():
        df_not_repeated = {}
        for i, row in parsed_psiblast[df].iterrows():
            if 'domtblout' in df:
                df_not_repeated.setdefault(row['ids'],[]).append([row['from_ali_coord'], row['to_ali_coord'], row['target_len']])
            else:
                df_not_repeated.setdefault(row['ids'],[]).append([row['hit_start'], row['hit_end'], row['protein_length']])
        
        list_dfs_not_repeated.append(df_not_repeated)

    dict_dfs = dict(zip(list_dfs, list_dfs_not_repeated))

    return dict_dfs


def new_logic_case_1(a, b):#, tp, tn, fp, fn):
    global tp
    global tn
    global fp
    global fn

    if a == 0 and b == 1:
        fp += 1
    
    elif a == 1 and b == 1:
        tp += 1

    elif a == 1 and b == 0:
        fn += 1
    
    else:
        tn += 1
    

def new_logic_case_2(a, b):#, tp_, tn_, fp_, fn_):
    global tp_
    global tn_
    global fp_
    global fn_

    if a == 0 and b == 1:
        fp_ += 1
    
    elif a == 1 and b == 1:
        tp_ += 1

    elif a == 1 and b == 0:
        fn_ += 1
    
    else:
        tn_ += 1


def compute_con_matrix_9(parsed_domtblouts, parsed_psiblast, gt):
    new_logic_case_1_vectorized = np.vectorize(new_logic_case_1)
    new_logic_case_2_vectorized = np.vectorize(new_logic_case_2)

    gt_acc = gt.accession.to_list()

    dict_dfs = dfs_to_dicts(parsed_domtblouts, parsed_psiblast)

    conf_matrix = []

    #conf_df = pd.DataFrame(columns=['true_positives', 'true_negatives', 'false_positives', 'false_negatives'])
    
    for k, (df_name, df_dict) in enumerate(dict_dfs.items()):
        gt_int_df = [x for x in gt_acc if x in list(df_dict.keys())]
        df_not_gt = [x for x in list(df_dict.keys()) if x not in gt_acc]
        gt_not_df = [x for x in gt_acc if x not in list(df_dict.keys())]

        global tp
        global tn
        global fp
        global fn
        global tp_
        global tn_
        global fp_
        global fn_
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        tp_ = 0
        tn_ = 0
        fp_ = 0
        fn_ = 0

        for j, (ids, lists) in enumerate(df_dict.items()):
            if ids in gt_int_df:
                row_df = np.zeros((int(lists[0][2]),), dtype=int)

                for i, el in enumerate(lists):
                    row_df[int(el[0]) : int(el[1])] = 1
                
                row_gt = np.zeros((int(lists[0][2]),), dtype=int)
                row_gt[gt[gt.accession == ids].iloc[0,1] : gt[gt.accession == ids].iloc[0,2]] = 1
                
                new_logic_case_1_vectorized(row_gt, row_df)#, tp, tn, fp, fn)

            elif ids in df_not_gt:
                row_df = np.zeros((int(lists[0][2]),), dtype=int)

                for i, el in enumerate(lists):
                    row_df[int(el[0]) : int(el[1])] = 1
                
                row_gt = np.zeros((int(lists[0][2]),), dtype=int)
                
                new_logic_case_2_vectorized(row_gt, row_df)#, tp_, tn_, fp_, fn_)
        
        if gt_not_df != []:
            for ids in gt_not_df:
                row_gt = np.zeros((gt[gt.accession == ids].iloc[0,3],), dtype=int)
                row_gt[gt[gt.accession == ids].iloc[0,1] : gt[gt.accession == ids].iloc[0,2]] = 1

                row_df = np.zeros((gt[gt.accession == ids].iloc[0,3],), dtype=int)

                new_logic_case_1_vectorized(row_gt, row_df)#, tp, tn, fp, fn)
                
        conf_matrix.append([tp + tp_, tn + tn_, fp + fp_, fn + fn_])
    
    #conf_df = pd.DataFrame(conf_matrix, columns=['true_positives', 'true_negatives', 'false_positives', 'false_negatives'])

    conf_df = pd.DataFrame(conf_matrix, index=list(dict_dfs.keys()), columns=['true_positives', 'true_negatives', 'false_positives', 'false_negatives'])

    return conf_df


def metrics_9(parsed_domtblouts, parsed_psiblast, gt):
    conf_df = compute_con_matrix_9(parsed_domtblouts, parsed_psiblast, gt)

    index_metrics = list(parsed_domtblouts.keys()) + list(parsed_psiblast.keys())
    columns_metrics = ['accuracy', 'precision', 'recall', 'specificity', 'balanced_accuracy', 'mcc', 'f1_score']

    metrics_list = []

    if 'metrics_9.csv' in os.listdir(cur + '\\data_team_1\\metrics'):
        old_metrics_df = pd.read_csv(cur + '\\data_team_1\\metrics\\metrics_9.csv', index_col=0)

        if old_metrics_df.index.to_list() == index_metrics:
            return old_metrics_df, conf_df
        
        else:
            for i, row in conf_df.iterrows():
                true_positives = row['true_positives']
                true_negatives = row['true_negatives']
                false_positives = row['false_positives']
                false_negatives = row['false_negatives']

                accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
            
                precision = true_positives / (true_positives + false_positives)

                recall = true_positives / (true_positives + false_negatives)

                specificity = true_negatives / (true_negatives + false_positives)

                balanced_accuracy = (recall + specificity) / 2

                mcc = ((true_positives * true_negatives) - (false_positives * false_negatives)) / np.sqrt((true_positives + false_positives) * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))

                f1_score = (2 * true_positives) / (2 * true_positives + false_positives + false_negatives)

                metrics_list.append([accuracy, precision, recall, specificity, balanced_accuracy, mcc, f1_score])
            
            metrics_df = pd.DataFrame(metrics_list, index=conf_df.index, columns=columns_metrics)
            
            metrics_df.to_csv(cur + '\\data_team_1\\metrics\\metrics_9.csv')

            return metrics_df, conf_df
    
    else:
        for i, row in conf_df.iterrows():
            true_positives = row['true_positives']
            true_negatives = row['true_negatives']
            false_positives = row['false_positives']
            false_negatives = row['false_negatives']

            accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
        
            precision = true_positives / (true_positives + false_positives)

            recall = true_positives / (true_positives + false_negatives)

            specificity = true_negatives / (true_negatives + false_positives)

            balanced_accuracy = (recall + specificity) / 2

            mcc = ((true_positives * true_negatives) - (false_positives * false_negatives)) / np.sqrt((true_positives + false_positives) * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))

            f1_score = (2 * true_positives) / (2 * true_positives + false_positives + false_negatives)

            metrics_list.append([accuracy, precision, recall, specificity, balanced_accuracy, mcc, f1_score])
        
        metrics_df = pd.DataFrame(metrics_list, index=conf_df.index, columns=columns_metrics)
        
        metrics_df.to_csv(cur + '\\data_team_1\\metrics\\metrics_9.csv')

        return metrics_df, conf_df


def plot_metrics_8(metrics_df):
    col= list(metrics_df.columns)
    x = np.arange(metrics_df.shape[0])
    c = np.random.rand(metrics_df.shape[0],3)
    for i in col:
        fig, ax = plt.subplots(figsize=(20,10))
        ax.bar(metrics_df.index,list(metrics_df[i]), color = c)
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_df.index, rotation=65)
        ax.set_ylim((0,1))
        plt.title(i)
        #fig.savefig(i+".png")

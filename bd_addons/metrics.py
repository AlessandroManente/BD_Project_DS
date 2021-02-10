from Bio import SeqIO, SearchIO
import pandas as pd
import os

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
    pass

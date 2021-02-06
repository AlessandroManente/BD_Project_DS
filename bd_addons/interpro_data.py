import pandas as pd
import urllib
import json
import time
import os

cur = os.getcwd()

metadata_columns = ['key', 'accession', 'name', 'source_database', 'length', 'source_organism_taxID', 'source_organism_scientificName', 'source_organism_fullName']
entries_columns = ['key', 'accession', 'entry_protein_locations_fragments_start', 'entry_protein_locations_fragments_end', 'entry_protein_locations_fragments_dc-status', 'entry_protein_locations_model', 'entry_protein_locations_score', 'protein_length', 'source_database', 'entry_type', 'entry_integrated']


def json_to_df(df_metadata, df_entries, json, key):
    i = key
    for el in json['results']:
        meta = []
        
        for key in el['metadata'].keys():
            if key != 'source_organism':
                meta.append(el['metadata'][key])
                
            else:
                for key_2 in el['metadata'][key].keys():
                    meta.append(el['metadata'][key][key_2])
        
        entry = []
        
        for key in el['entries'][0]:
            if key != 'entry_protein_locations':
                entry.append(el['entries'][0][key])
                
            else:
                for key_2 in el['entries'][0][key][0].keys():
                    if key_2 != 'fragments':
                        entry.append(el['entries'][0][key][0][key_2])

                    else:
                        for key_3 in el['entries'][0][key][0][key_2][0].keys():
                            entry.append(el['entries'][0][key][0][key_2][0][key_3])
                        
        df_metadata = df_metadata.append(dict(zip(metadata_columns, [i] + meta)), ignore_index=True)
        df_entries = df_entries.append(dict(zip(entries_columns, [i] + entry)), ignore_index=True)

        i += 1

    return df_metadata, df_entries


def get_next_json(url):
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    return data


def get_df(url, qty, flag, team):
    data = get_next_json(url)

    metadata = pd.DataFrame(columns=metadata_columns)
    entries = pd.DataFrame(columns=entries_columns)
    

    while data['next'] is not None:
        try:
            if len(metadata) > 0:
                i = metadata.key.to_list()[-1] + 1
            else:
                i = 0

            remaining = qty - i
            print('Remaining:', remaining, 'proteins')

            metadata, entries = json_to_df(metadata, entries, data, i)
            
            data = get_next_json(data['next'])

            if data['next'] is None:
                if len(metadata) > 0:
                    i = metadata.key.to_list()[-1] + 1
                else:
                    i = 0

                metadata, entries = json_to_df(metadata, entries, data, i)
        
        except:
            checkpoint(metadata, entries, flag, team, remaining)
            time.sleep(900)

            if len(metadata) > 0:
                i = metadata.key.to_list()[-1] + 1
            else:
                i = 0

            remaining = qty - i
            print('Remaining:', remaining, 'proteins')

            metadata, entries = json_to_df(metadata, entries, data, i)
            
            data = get_next_json(data['next'])

            if data['next'] is None:
                if len(metadata) > 0:
                    i = metadata.key.to_list()[-1] + 1
                else:
                    i = 0

                metadata, entries = json_to_df(metadata, entries, data, i)
    
    checkpoint(metadata, entries, flag, team, remaining)

    return metadata, entries


def get_data(reviewed, unreviewed, general, team):
    general = get_next_json(general)
    general_reviewed = general['proteins']['reviewed']
    general_unreviewed = general['proteins']['unreviewed']

    metadata_reviewed, entries_reviewed = get_df(reviewed, general_reviewed, 'reviewed', team)
    print('Got reviewed data')
    metadata_unreviewed, entries_unreviewed = get_df(unreviewed, general_unreviewed, 'unreviewed', team)
    print('Got unreviewed data')

    return metadata_reviewed, entries_reviewed, metadata_unreviewed, entries_unreviewed

def checkpoint(df_metadata, df_entries, flag, team, remaining):
    df_metadata.to_csv(cur+'\\data_team_'+str(team)+'\\'+flag+'\\metadata\\metadata_'+flag+'_'+str(remaining)+'.csv')
    df_entries.to_csv(cur+'\\data_team_'+str(team)+'\\'+flag+'\\entries\\entries_'+flag+'_'+str(remaining)+'.csv')
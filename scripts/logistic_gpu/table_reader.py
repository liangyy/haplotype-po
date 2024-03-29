# to load phenotype and covariates from YAML:
# path: 'path_to_file'
# param:   # kwargs in pandas.read_csv
#     header: 0
#     sep: ','
# col:  # column name of phenotype of interest
#     - 'alzheimer_disease'
#     - 'heart_disease'
# indiv_col: 'eid'  # column name of individual id

import sys
import pandas as pd
import numpy as np
sys.path.insert(0, '../prs')
import gwas_reader

def load_table_from_yaml(yaml, rename_cols=False):
    '''
    Load the table specified in YAML.
    Rename individual ID column (specified in `indiv_col`)
    as `individual_id`
    '''
    if yaml is None:
        return pd.DataFrame({})
    
    file_dict = gwas_reader.read_yaml(yaml)
    file_read_kw = gwas_reader.clean_tab(file_dict['param'])
    df = pd.read_csv(
        file_dict['path'],
        **file_read_kw
    )
    if rename_cols is False:
        df = df[ file_dict['col'] + [file_dict['indiv_col']] ]
    else:
        df = df[ list(file_dict['col'].keys()) + [file_dict['indiv_col']] ]
        df = df.rename(columns=file_dict['col'])
    df = df.rename(columns={ file_dict['indiv_col']: 'individual_id' })
    df['individual_id'] = df['individual_id'].astype(str)
    return df

def load_file_as_lines(filename, skip=0):
    '''
    Read filename line by line and append to a list.
    Skip the first [skip] rows.
    '''
    o = []
    with open(filename, 'r') as f:
        for i in f:
            if skip > 0:
                skip -= 1
                continue
            o.append(i.strip())
    return o

def standardize_columns(df, except_cols=[]):
    '''
    Standardize columns except for the ones in except_cols
    '''        
    cols = df.columns.tolist()
    for cc in cols:
        if cc in except_cols:
            continue
        else:
            df[cc] = _standardize(df[cc])    
    return df

def _standardize(dat):
    # breakpoint()
    dat = (dat - np.mean(dat)) / np.std(dat)
    return dat


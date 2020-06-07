import argparse
parser = argparse.ArgumentParser(prog='prs_based_transcriptome_wide_imputation.py', description='''
    Script to run PRS-based imputation using Framingham transcriptome data. 
''')

parser.add_argument('--obs_expr', help='''
    Observed expression (indiv x gene).
    Add sample ID column name by: Filename:SampleCol
''')
parser.add_argument('--pred_expr', type=str, help='''
    Predicted expression (haplo-indiv x gene).
    Add sample ID column name and add excluding ID columns by: Filename:SampleCol:ExcludeCol1,ExcludeCol2 
''')
parser.add_argument('--pedigree', help='''
    Pedigree (one family per row).
    Added child ID, father ID, and mother ID by:
    Filename:ChildCol,FatherCol,MotherCol
''')
parser.add_argument('--covar', help='''
    Covariates (indiv x gene)
    Add sample ID column name by: Filename:SampleCol
''')
parser.add_argument('--output', help='''
    Output in TSV.GZ format.
''')
parser.add_argument('--nthread', default=None, type=int, help='''
    Number of threads to use.
''')

args = parser.parse_args()




import logging, sys
import torch
import pandas as pd
import sys, os
sys.path.insert(0, '../../haplotype_imputation')
import haplotype_imputer

def read_table(fname):
    if os.path.splitext(fname) == 'gz':
        kw = {'compression' : 'gzip'}
    else:
        kw = {}
    return pd.read_csv(fname, sep='\t', **kw)

def load_table_from_str(filestr, keep_only=False):
    filename = filestr.split(':')[0]
    df = read_table(filename)
    if keep_only is False:
        indiv_col = filestr.split(':')[1]
        df = df.rename({indiv_col : 'individual_id'})
    elif keep_only is True:
        keep_cols = filestr.split(':')[1].split(',')
        df = df[keep_cols]
    return df

def remove_dot(str_):
    return str_.split('.')[0]

def load_pred_expr_from_str(filestr):
    filestr_to_read = ':'.join(filestr.split(':')[:2])
    df = load_table_from_str(filestr_to_read)
    df = df.rename(columns={'hap_indiv_id' : 'individual_id'})
    
    to_drop_cols = filestr.split(':')[-1].split(',')
    for i in to_drop_cols:
        del df[i]
    
    df['suffix'] = [ i.split('_')[-1] for i in df['hap_indiv_id'].tolist() ]  # temporary column
    df['individual_id'] = [ i.split('_')[0] for i in df['hap_indiv_id'].tolist() ]
    df_h1 = df[ df['suffix'] == 'h1' ].reset_index(drop=True).copy()
    df_h2 = df[ df['suffix'] == 'h2' ].reset_index(drop=True).copy()
    del df_h1['suffix']  # remove the temporary column
    del df_h2['suffix']  # remove the temporary column
    
    # remove suffix in gene name 
    df_h1.columns = [ remove_dot(i) for i in df_h1.columns.tolist() ]
    df_h2.columns = [ remove_dot(i) for i in df_h2.columns.tolist() ]

    return df_h1, df_h2

def get_elements_in_common(alist, blist):
    return list(set(alist).intersection(set(blist)))

def extract_by_cols(df_list, cols_to_extract):
    out = []
    for df in df_list:
        out.append(out[cols_to_extract])
    return out

def get_rows(df, rows, by_col):
    df_ref = pd.DataFrame({'row': rows})
    df_ref = pd.merge(df_ref, df, left_on='row', right_on=by_col, how='inner')
    del df_ref['row']
    return df_ref

def extract_by_rows(df_, rows_to_extract, by_col, is_list=True):
    if is_list is True:
        return get_rows(df_, rows_to_extract, by_col)
    elif is_list is False:
        out = []
        for df in df_:
            out.append(get_rows(df, rows_to_extract, by_col))
        return out
    else:
        raise ValueError('Unexpected: is_list = ', is_list)

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

logging.info('Loading observed expression')
df_obs_expr = load_table_from_str(args.obs_expr)

logging.info('Loading covariates')
df_covar = load_table_from_str(args.covar)

logging.info('Loading pedigree')
df_ped = load_table_from_str(args.pedigree, keep_only=True)
df_ped.columns = ['individual_id', 'father', 'mother']

logging.info('Loading predicted expression')
df_h1, df_h2 = load_pred_expr_from_str(args.pedigree)

logging.info('Getting genes in common')
genes_in_common = get_elements_in_common(
    df_obs_expr.columns.tolist()[1:],
    df_h1.columns.tolist()[1:]
)
logging.info('--> There are {} genes in common'.format(len(genes_in_common)))
df_obs_expr, df_h1, df_h2 = extract_by_cols([ df_obs_expr, df_h1, df_h2 ], [ 'individual_id' ] + genes_in_common)

logging.info('Extracting individuals in pedigree')
df_obs_expr_father = extract_by_rows(df_obs_expr, df_ped['father'].tolist(), is_list=False)
df_obs_expr_mother = extract_by_rows(df_obs_expr, df_ped['mother'].tolist(), is_list=False)
df_h1, df_h2, df_covar = extract_by_rows([ df_h1, df_h2, df_covar ], df_ped['individual_id'].tolist())

mode = 'basic_em_py'
logging.info(f'Run imputation: mode = {mode}')
imputer = haplotype_imputer.HaploImputer()
beta, sigma2, out, lld = imputer.impute(
    df_father, df_mother, 
    h1, h2, hap_indiv_df, hap_pos_df,
    mode=mode,
    df_covar=df_covar,
    kwargs={'debug_cache': args.debug_cache_prefix}
)

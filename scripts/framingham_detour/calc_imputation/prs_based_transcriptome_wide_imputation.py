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
parser.add_argument('--imputer-output', type=str, help='''
    Pickle GZ imputer output
''')
parser.add_argument('--nthread', default=None, type=int, help='''
    Number of threads to use.
''')
parser.add_argument('--downsample', default=None, type=float, help='''
    Down-sample the genes. 
    Take a fraction which will retain only the fraction of genes.
''')

args = parser.parse_args()




import logging, sys
import torch
import pandas as pd
import numpy as np
import sys, os
import gzip, pickle
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
        df = df.rename(columns={indiv_col: 'individual_id'})
        df['individual_id'] = df['individual_id'].astype(str)
    elif keep_only is True:
        keep_cols = filestr.split(':')[1].split(',')
        df = df[keep_cols]
    return df

def remove_dot(str_):
    return str_.split('.')[0]

def load_pred_expr_from_str(filestr):
    filestr_to_read = ':'.join(filestr.split(':')[:2])
    df = load_table_from_str(filestr_to_read)
    df = df.rename(columns={'individual_id' : 'hap_indiv_id'})
    # breakpoint()
    
    to_drop_cols = filestr.split(':')[-1].split(',')
    for i in to_drop_cols:
        del df[i]
    
    df['suffix'] = [ i.split('_')[-1] for i in df['hap_indiv_id'].tolist() ]  # temporary column
    df['individual_id'] = [ i.split('_')[0] for i in df['hap_indiv_id'].tolist() ]
    df_h1 = df[ df['suffix'] == 'h1' ].reset_index(drop=True).copy()
    df_h2 = df[ df['suffix'] == 'h2' ].reset_index(drop=True).copy()
    del df_h1['suffix']  # remove the temporary column
    del df_h2['suffix']  # remove the temporary column
    del df_h1['hap_indiv_id']  # remove the column since we don't need them anymore
    del df_h2['hap_indiv_id']  # remove the column since we don't need them anymore
    
    # remove suffix in gene name 
    df_h1.columns = [ remove_dot(i) for i in df_h1.columns.tolist() ]
    df_h2.columns = [ remove_dot(i) for i in df_h2.columns.tolist() ]
    
    # deprecated since we want to filter gene std = 0 after subset on individuals
    # # remove genes with std = 0
    # # breakpoint()
    # df_h1 = df_h1.loc[:, [True] + (df_h1.std() != 0).tolist() + [True] ]
    # df_h2 = df_h2.loc[:, [True] + (df_h2.std() != 0).tolist() + [True] ]
    
    return df_h1, df_h2

def get_elements_in_common(alist, blist):
    return list(set(alist).intersection(set(blist)))

def extract_by_cols(df_list, cols_to_extract):
    out = []
    for df in df_list:
        out.append(df[cols_to_extract])
    return out

def get_rows(df, rows, by_col):
    df_ref = pd.DataFrame({'row': rows})
    # breakpoint()
    df_ref = pd.merge(df_ref, df, left_on='row', right_on=by_col, how='inner')
    del df_ref['row']
    return df_ref

def extract_by_rows(df_, rows_to_extract, by_col, is_list=True):
    if is_list is False:
        return get_rows(df_, rows_to_extract, by_col)
    elif is_list is True:
        out = []
        for df in df_:
            out.append(get_rows(df, rows_to_extract, by_col))
        return out
    else:
        raise ValueError('Unexpected: is_list = ', is_list)

def downsample(clist, fraction):
    n = int(len(clist) * fraction)
    selected_idx = np.random.choice(len(clist), n, replace=False)
    return list(np.array(clist)[selected_idx])

# set number of threads to use
if args.nthread is not None:
    torch.set_num_threads(args.nthread)

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
df_h1, df_h2 = load_pred_expr_from_str(args.pred_expr)

logging.info('Extracting individuals in pedigree')
# breakpoint()
df_obs_expr_father = extract_by_rows(df_obs_expr, df_ped['father'].astype(str).tolist(), 'individual_id', is_list=False)
df_obs_expr_mother = extract_by_rows(df_obs_expr, df_ped['mother'].astype(str).tolist(), 'individual_id', is_list=False)
df_obs_expr_father['individual_id'] = df_ped['individual_id'].astype(str).tolist()
df_obs_expr_mother['individual_id'] = df_ped['individual_id'].astype(str).tolist()
df_h1, df_h2, df_covar = extract_by_rows([ df_h1, df_h2, df_covar ], df_ped['individual_id'].astype(str).tolist(), 'individual_id')


logging.info('Getting genes in common')
genes_in_common = get_elements_in_common(
    df_obs_expr.columns.tolist()[1:],
    df_h1.columns.tolist()[1:]
)
logging.info('--> There are {} genes in common'.format(len(genes_in_common)))
df_obs_expr, df_h1, df_h2 = extract_by_cols([ df_obs_expr, df_h1, df_h2 ], [ 'individual_id' ] + genes_in_common)

if args.downsample is not None:
    logging.info('Downsample: fraction = {}'.format(args.downsample))
    gene_subset = downsample(genes_in_common, args.downsample)
    df_obs_expr, df_h1, df_h2 = extract_by_cols([ df_obs_expr, df_h1, df_h2 ], [ 'individual_id' ] + gene_subset)
    genes_in_common = gene_subset
    
# filter gene with std = 0
df_h1 = df_h1.iloc[:, [True] + (df_h1.std() != 0).tolist() ]
df_h2 = df_h2.iloc[:, [True] + (df_h2.std() != 0).tolist() ]
logging.info('SUMMARY (after removing gene with std = 0): {} genes are used'.format()

mode = 'basic_em_py'
logging.info(f'Run imputation: mode = {mode}')
imputer = haplotype_imputer.HaploImputer()
# breakpoint()
beta, sigma2, out, lld = imputer.impute(
    df_obs_expr_father, df_obs_expr_mother, 
    df_h1, df_h2,
    mode=mode,
    df_covar=df_covar,
    # kwargs={'debug_cache': args.debug_cache_prefix}
)

logging.info('Save imputer output')
with gzip.open(args.imputer_output, 'w') as f:
    pickle.dump((beta, sigma2, lld), f)


logging.info('Output')
out.to_csv(args.output, compression='gzip', sep='\t', index=False)


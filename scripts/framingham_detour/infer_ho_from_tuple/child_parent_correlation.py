import argparse
parser = argparse.ArgumentParser(prog='child_parent_correlation.py', description='''
    Compute correlation between child haplotypes and parent genotypes. 
    Also, father genotype and mother genotype.
''')

parser.add_argument('--h1', help='''
    haplotype 1
''')
parser.add_argument('--h2', default=None, help='''
    haplotype 2
''')
parser.add_argument('--pedigree', help='''
    pedigree file
''')
parser.add_argument('--maf-filter', default=0, type=float, help='''
    MAF filter to remove rare variants
''')
parser.add_argument('--child_col', help='''
    column name in pedigree file that indicates child
''')
parser.add_argument('--father_col', help='''
    column name in pedigree file that indicates father
''')
parser.add_argument('--mother_col', help='''
    column name in pedigree file that indicates mother
''')
parser.add_argument('--output', help='''
    a table in TSV.GZ
''')
args = parser.parse_args()

import logging, time, sys
import pandas as pd
import numpy as np
from tqdm import tqdm 

sys.path.insert(0, '../')
import misc_helper

def cor_dist(vec1, vec2):
    '''
    vec1, vec2 are pandas Series
    '''
    return np.corrcoef(vec1, vec2)[0, 1]

def inner_prod_dist(vec1, vec2):
    return (vec1 * vec2).mean()

def calc_pairwise_product(vec_list1, vec_list2, label_list1, label_list2):
    out_dict = {}
    for vec1, name1 in zip(vec_list1, label_list1):
        for vec2, name2 in zip(vec_list2, label_list2):
            # breakpoint()
            out_dict[f'{name1}_x_{name2}'] = inner_prod_dist(vec1, vec2)
    return out_dict

def _standardize(mat):
    mean = mat.mean()
    std = mat.std()
    return (mat - mean) / std

def standardize_rows(df):
    '''
    df is pandas DataFrame
    '''
    mat = df.values
    mat = np.apply_along_axis(_standardize, 1, mat)
    mat = pd.DataFrame(mat)
    mat.columns = df.columns
    return mat
    

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# load haplotypes
logging.info('Loading haplotypes')
df_h1 = pd.read_parquet(args.h1).astype(float)
df_h2 = pd.read_parquet(args.h2).astype(float)

# apply maf filter
logging.info('Applying maf filters: maf = {}'.format(args.maf_filter))
maf = (df_h1 + df_h2).apply(lambda x: min(x.mean() / 2, (1 - x).mean() / 2), axis=1)
df_h1 = df_h1[ maf > args.maf_filter ].reset_index(drop=True)
df_h2 = df_h2[ maf > args.maf_filter ].reset_index(drop=True)

# standardizing genotype
logging.info('Standardizing haplotypes')
# breakpoint()
df_h1 = standardize_rows(df_h1)
df_h2 = standardize_rows(df_h2)
# breakpoint()

# load pedigree data
logging.info('Loading pedigree data')
df_ped = misc_helper.read_table(args.pedigree)
df_ped = df_ped[ [ args.child_col, args.father_col, args.mother_col ] ]
for i in df_ped.columns:
    df_ped[i] = df_ped[i].astype(str)
df_ped.columns = [ 'child', 'father', 'mother' ]

# naive infer haplotype origin
clabel = [ f'child{i}' for i in range(1, 3) ]
plabel = [ 'father', 'mother' ]
df_list = []  # collect the resulting tables
for i in tqdm(range(df_ped.shape[0])):
    cid = df_ped.child[i]
    fid = df_ped.father[i]
    mid = df_ped.mother[i]
    pw_cor_child_x_parents = calc_pairwise_product(
        [ df_h1[cid], df_h2[cid] ], 
        [ (df_h1[fid] + df_h2[fid]) / np.sqrt(2), (df_h1[mid] + df_h2[mid]) / np.sqrt(2) ],
        clabel, plabel
    )
    pw_cor_father_x_mother = calc_pairwise_product(
        [ (df_h1[fid] + df_h2[fid]) / np.sqrt(2), (df_h1[mid] + df_h2[mid]) / np.sqrt(2) ],
        [ (df_h1[fid] + df_h2[fid]) / np.sqrt(2), (df_h1[mid] + df_h2[mid]) / np.sqrt(2) ],
        [ plabel[0], plabel[1] ], [ plabel[0], plabel[1] ]
    )
    df_ = pd.concat( 
        [ pd.DataFrame(pw_cor_child_x_parents, index=[0]),
          pd.DataFrame(pw_cor_father_x_mother, index=[0]) ],
        axis=1
    )
    df_['individual_id'] = cid
    df_list.append(df_)
df_to_save = pd.concat(df_list, axis=0)
df_to_save.to_csv(args.output, compression='gzip', sep='\t', index=False)


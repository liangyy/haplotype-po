import argparse
parser = argparse.ArgumentParser(prog='infer_ho_from_haplotype_in_family.py', description='''
    Infer haplotype origin from haplotypes from child, father, and mother.
    Output the most probable assignment in three main entries. 
    1. if haplotype 1 in child comes from father.
    2. the origin of haplotype 1 in child.
    3. the origin of haplotype 2 in child.
    And some other information about the matching:
    pairwise distance between child 1/2 and father 1/2;
    and pairwise distance between child 1/2 and mother 1/2
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
from tqdm import tqdm 

sys.path.insert(0, '../')
import misc_helper

def manhattan_dist(vec1, vec2):
    '''
    vec1, vec2 are pandas Series
    '''
    return (vec1 != vec2).sum()

def calc_pairwise_dist(vec_list1, vec_list2, label_list1, label_list2):
    out_dict = {}
    for vec1, name1 in zip(vec_list1, label_list1):
        for vec2, name2 in zip(vec_list2, label_list2):
            out_dict[f'{name1}_x_{name2}'] = manhattan_dist(vec1, vec2)
    return out_dict

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

# load haplotypes
logging.info('Loading haplotypes')
df_h1 = pd.read_parquet(args.h1)
df_h2 = pd.read_parquet(args.h2)

# load pedigree data
logging.info('Loading pedigree data')
df_ped = misc_helper.read_table(args.pedigree)
df_ped = df_ped[ [ args.child_col, args.father_col, args.mother_col ] ]
for i in df_ped.columns:
    df_ped[i] = df_ped[i].astype(str)
df_ped.columns = [ 'child', 'father', 'mother' ]

# naive infer haplotype origin
clabel = [ f'child{i}' for i in range(1, 3) ]
flabel = [ f'father{i}' for i in range(1, 3) ]
mlabel = [ f'mother{i}' for i in range(1, 3) ]
all_configs_h1_from_father = [ 
    f'child1_x_father{i}:child2_x_mother{j}' for j in range(1, 3) for i in range(1, 3) 
]
all_configs_h1_from_mother = [ 
    f'child2_x_father{i}:child1_x_mother{j}' for j in range(1, 3) for i in range(1, 3) 
]
all_configs = all_configs_h1_from_father + all_configs_h1_from_mother
df_list = []  # collect the resulting tables
for i in tqdm(range(df_ped.shape[0])):
    cid = df_ped.child[i]
    fid = df_ped.father[i]
    mid = df_ped.mother[i]
    pw_child_x_father = calc_pairwise_dist(
        [ df_h1[cid], df_h2[cid] ], [ df_h1[fid], df_h2[fid] ],
        clabel, flabel
    )
    pw_child_x_mother = calc_pairwise_dist(
        [ df_h1[cid], df_h2[cid] ], [ df_h1[mid], df_h2[mid] ],
        clabel, mlabel
    )
    config_dist = []
    config_name = []
    config_prob_z = []
    for conf_str in all_configs:
        conf_f, conf_m = conf_str.split(':')
        dist_ = pw_child_x_father[conf_f] + pw_child_x_mother[conf_m]
        config_dist.append(dist_)
        config_name.append(f'{conf_f}:{conf_m}')
        if conf_str in all_configs_h1_from_father:
            config_prob_z.append(1)
        elif conf_str in all_configs_h1_from_mother:
            config_prob_z.append(0)
        else:
            raise ValueError(f'conf_str = {conf_str} is not expected')
    df_ = pd.DataFrame({'dist': config_dist, 'name': config_name, 'prob_z': config_prob_z})
    min_dist = df_['dist'].min()
    # breakpoint()
    df_ = df_[ df_['dist'] == min_dist ].reset_index(drop=True).copy()
    # breakpoint()
    df_report = pd.concat([ pd.concat((pd.DataFrame(pw_child_x_father, index=[0]), pd.DataFrame(pw_child_x_mother, index=[0])), axis=1) ] * df_.shape[0])
    df_report = pd.concat((df_report, df_), axis=1)
    df_report['individual_id'] = cid
    df_list.append(df_report)
df_to_save = pd.concat(df_list, axis=0)
df_to_save.to_csv(args.output, compression='gzip', sep='\t', index=False)

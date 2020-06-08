import argparse
parser = argparse.ArgumentParser(prog='phased_genotype_to_parquet.py', description='''
    Convert phased genotype (one row is one variant) to 
    haplotype table (one row is one variant) in parquet format
''')

parser.add_argument('--genotype', help='''
    genotype file in tsv.gz format
''')
parser.add_argument('--samples', default=None, help='''
    if the header is missing in the genotype file, specify the sample list as header
''')
parser.add_argument('--output_prefix', help='''
    prefix of output
''')
args = parser.parse_args()

import logging, time, sys
import gzip
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, '../')
import misc_helper

def parse_phased_geno(txt):
    return txt.split('|')

# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

h1_mat = []
h2_mat = []
with gzip.open(args.genotype, 'rt') as f:
    for i in tqdm(f):
        h1_vec = []
        h2_vec = []
        i = i.strip().split('\t')
        for k in i:
            h1, h2 = parse_phased_geno(k)
            h1_vec.append(h1)
            h2_vec.append(h2)
        h1_mat.append(h1_vec)   
        h2_mat.append(h2_vec)

sample_list = []
with open(args.samples, 'r') as f:
    for i in f:
        sample_list.append(i.strip())

df_h1 = pd.DataFrame(h1_mat)
df_h1.columns = sample_list
df_h2 = pd.DataFrame(h2_mat)
df_h2.columns = sample_list

df_h1.to_parquet(args.output_prefix + '.h1.parquet')
df_h2.to_parquet(args.output_prefix + '.h2.parquet')


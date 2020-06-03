import argparse
parser = argparse.ArgumentParser(prog='merge_pred_expr.py', description='''
    merge predicted expression across all chromosomes
''')

parser.add_argument('--input-pattern', help='''
    input file pattern (with wildcards {chr_num})
''')
parser.add_argument('--output', help='''
    output txt.gz
''')
args = parser.parse_args()

import logging, time, sys
# configing util
logging.basicConfig(
    level = logging.INFO, 
    stream = sys.stderr, 
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

import pandas as pd


logging.info('Merge starts')
df = pd.read_csv(args.input_pattern.format(chr_num=1), sep='\t')
logging.info('Chromosome 1 finishes (initialize)')
for i in range(2, 23):
    df_ = pd.read_csv(args.input_pattern.format(chr_num=i), sep='\t')
    df = pd.merge(df, df_, left_on=['FID', 'IID'], right_on=['FID', 'IID'], how='inner')
    logging.info(f'Chromosome {i} finishes')
df.to_csv(args.output, compression='gzip', sep='\t', index=False)

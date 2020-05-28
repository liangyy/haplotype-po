import argparse
parser = argparse.ArgumentParser(prog='annot_rsid.py', description='''
    Add rsID to Neale lab UKBB GWAS
''')

parser.add_argument('--gwas', help='''
    Neale lab UKBB GWAS summary statistics
''')
parser.add_argument('--variant-meta', help='''
    Meta-information of the variants in Neale lab UKBB GWAS
''')
parser.add_argument('--output', help='''
    Output
''')
args = parser.parse_args()

import logging, time, sys, os
# configing util
logging.basicConfig(
    level = logging.INFO,
    stream = sys.stderr,
    format = '%(asctime)s  %(message)s',
    datefmt = '%Y-%m-%d %I:%M:%S %p'
)

import pandas as pd


logging.info('Loading GWAS')
gwas = pd.read_csv(args.gwas, sep='\t', compression='gzip')
logging.info('Loading meta information of variants')
var_meta = pd.read_csv(args.variant_meta, sep='\t', compression='gzip', usecols=['variant', 'rsid', 'ref', 'alt'])
logging.info('Merging')
gwas = pd.merge(gwas, var_meta[['variant', 'rsid', 'ref', 'alt']], left_on='variant', right_on='variant', how='inner')
gwas.rename(
    columns = {
        'ref':'non_effect_allele',
        'alt': 'effect_allele'
    }, 
    inplace = True
)
logging.info('Output')
gwas.to_csv(args.output, sep='\t', compression='gzip', index=False)

